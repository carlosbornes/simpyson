import numpy as np
from simpyson.converter import hz2ppm, ppm2hz

class Simpy:
    """
    A unified container Simpson NMR data that automatically handles conversions between formats.
    
    This class stores and manages Simpson's various formasts (FID, Spe in Hz and ppm) and 
    conversion between them.
    
    Attributes:
        b0: Magnetic field (e.g., '9.4T', '400MHz')
        nucleus: Nucleus type (e.g., '1H', '13C')
        
    Properties:
        fid: Time-domain data
        spe: Frequency-domain data (Hz)
        ppm: Chemical shift data (ppm)
    """
    def __init__(self, b0=None, nucleus=None):
        self._b0 = b0
        self._nucleus = nucleus
        self._fid_data = None  # Time-domain data
        self._spe_data = None  # Frequency-domain data
        self._xreim_data = None # xreim data
        self._metadata = {}    
    
    @property
    def b0(self):
        return self._b0
    
    @b0.setter
    def b0(self, value):
        self._b0 = value
        # Invalidate cached ppm data when b0 change
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']
    
    @property
    def nucleus(self):
        return self._nucleus
    
    @nucleus.setter
    def nucleus(self, value):
        self._nucleus = value
        # Invalidate cached ppm data if nucleus changs
        if self._spe_data and 'ppm' in self._spe_data:
            del self._spe_data['ppm']
    
    @property
    def fid(self):
        """Access time-domain data, converting from spectrum if needed."""
        if self._fid_data is None and self._spe_data is not None:
            self._compute_fid()
        return self._fid_data
    
    @property
    def spe(self):
        """Access frequency-domain data, converting from FID if needed."""
        if self._spe_data is None and self._fid_data is not None:
            self._compute_spectrum()
        return self._spe_data
    
    @property
    def ppm(self):
        """Access chemical shift data, computing from Hz if needed."""
        if self.spe and 'ppm' not in self.spe and self._b0 and self._nucleus:
            self._compute_ppm()
        
        if self.spe and 'ppm' in self.spe:
            return {
                'ppm': self.spe['ppm'],
                'real': self.spe['real'],
                'imag': self.spe['imag']
            }
        return None
    
    @property
    def xreim(self):
        """Access xreim data."""
        return self._xreim_data

    def _compute_spectrum(self):
        """Convert FID to spectrum."""
        if not self._fid_data:
            return
        
        npoints = self._fid_data['np']
        sw = self._fid_data['sw']
        signal = self._fid_data['real'] + 1j * self._fid_data['imag']
        spectrum = np.fft.fftshift(np.fft.fft(signal))
        
        self._spe_data = {
            'real': np.real(spectrum),
            'imag': np.imag(spectrum),
            'np': npoints,
            'sw': sw,
            'hz': np.linspace(-sw/2, sw/2, int(npoints))
        }
    
    def _compute_fid(self):
        """Convert spectrum to FID."""
        if not self._spe_data:
            return
            
        npoints = self._spe_data['np']
        sw = self._spe_data['sw']
        signal = self._spe_data['real'] + 1j * self._spe_data['imag']
        time_signal = np.fft.ifft(np.fft.ifftshift(signal))
        dt = 1.0 / sw
        
        self._fid_data = {
            'real': np.real(time_signal),
            'imag': np.imag(time_signal),
            'np': npoints,
            'sw': sw,
            'time': np.linspace(0, npoints*dt, int(npoints)) * 10e3
        }
    
    def _compute_ppm(self):
        """Calculate ppm scale from Hz."""
        if not (self._spe_data and 'hz' in self._spe_data):
            return
        
        if not (self._b0 and self._nucleus):
            return
            
        try:
            self._spe_data['ppm'] = hz2ppm(self._spe_data['hz'], self._b0, self._nucleus)
        except ValueError as e:
            print(f"Error converting to ppm: {e}")
    
    def from_fid(self, real, imag, np_value, sw, time=None):
        """Set data from FID values."""
        if time is None:
            dt = 1.0 / sw
            time = np.linspace(0, np_value*dt, int(np_value)) * 10e3
            
        self._fid_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'time': np.array(time)
        }
        
        # Clear cached spectrum data
        self._spe_data = None
        return self
    
    def from_spe(self, real, imag, np_value, sw, hz=None):
        """Set data from spectrum values."""
        if hz is None:
            hz = np.linspace(-int(sw) / 2, int(sw) / 2, int(np_value))
            
        self._spe_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'hz': np.array(hz)
        }
        
        # Calculate ppm if possible
        if self._b0 and self._nucleus:
            self._compute_ppm()
            
        # Clear cached FID data
        self._fid_data = None
        return self
    
    def from_xreim(self, time, real, imag):
        """Set data from xreim values."""

        self._xreim_data = {
            'time': np.array(time),
            'real': np.array(real),
            'imag': np.array(imag)
        }

        # Clear cached FID and spectrum data
        self._fid_data = None
        self._spe_data = None
        return self
    
    def from_csdf(self, real, imag, hz, np_value, sw):
        """Set data from CSDF values."""

        self._spe_data = {
            'real': np.array(real),
            'imag': np.array(imag),
            'np': np_value,
            'sw': sw,
            'hz': np.array(hz)
        }

        if self._b0 and self._nucleus:
            self._compute_ppm()

        self._fid_data = None
        return self
    
    def copy(self):
        """Create a copy of a Simpy object."""
        import copy as cp
    
        new_obj = Simpy(b0=self._b0, nucleus=self._nucleus)

        if self._fid_data is not None:
            new_obj._fid_data = cp.deepcopy(self._fid_data)

        if self._spe_data is not None:
            new_obj._spe_data = cp.deepcopy(self._spe_data)

        if self._xreim_data is not None:
            new_obj._xreim_data = cp.deepcopy(self._xreim_data)

        new_obj._metadata = cp.deepcopy(self._metadata)

        return new_obj
    
    def write(self, filename, format='csv'):
        """
        Write data to file in specified format.
        
        Args:
            filename: Output filename
            format: Format to save as ('csv', 'fid', 'spe')
            
        Returns:
            self for method chaining
        """
        # CSV format
        if format == 'csv':
            if self.spe:
                if 'ppm' in self.spe and self._b0 and self._nucleus:
                    x_data, x_label = self.spe['ppm'], 'ppm'
                else:
                    x_data, x_label = self.spe['hz'], 'Hz'
                y_data = self.spe['real']
            elif self.fid:
                x_data, x_label = self.fid['time'], 'Time (ms)'
                y_data = self.fid['real']
            elif self.xreim:
                x_data, x_label = self.xreim['time'], 'Time (ms)'
                y_data = self.xreim['real']
            else:
                raise ValueError("No data to save")
                
            np.savetxt(
                filename,
                np.column_stack((x_data, y_data)),
                delimiter=",",
                header=f"{x_label},Real",
                comments=""
            )
            return self

        # SIMPSON formats
        data_dict = None
        if format == 'spe':
            data_dict = self.spe
            data_type = 'SPE'
        elif format == 'fid':
            data_dict = self.fid
            data_type = 'FID'
        elif format == 'xreim':
            data_dict = self.xreim
            data_type = 'XREIM'
        else:
            raise ValueError(f"Unsupported save format: {format}")
            
        if not data_dict:
            raise ValueError(f"No data available to save in {format} format")
            
        # Write SIMPSON format file
        with open(filename, 'w') as f:
            f.write('SIMP\n')
            if 'np' in data_dict:
                f.write(f'NP={data_dict["np"]}\n')
            if 'sw' in data_dict:
                f.write(f'SW={data_dict["sw"]}\n')
            f.write(f'TYPE={data_type}\n')
            f.write('DATA\n')
            
            for re, im in zip(data_dict['real'], data_dict['imag']):
                f.write(f'{re} {im}\n')
                
            f.write('END')
            
        return self
