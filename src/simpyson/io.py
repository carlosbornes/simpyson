# This module contains functions for analyzing NMR data from SIMPSON simulations.
#
# It provides functions for reading NMR data from SPE and FID files, performing
# Fourier transforms, and plotting NMR spectra.

import numpy as np
import json
import os
import sys
from soprano.calculate.nmr.simpson import write_spinsys
import copy

class SimpReader:
    """
    A class to read and process NMR data from SIMPSON files.

    Attributes:
        filename (str): The name of the file to read.
        format (str): The format of the file (spe, fid, xreim).
        b0 (str, optional): The magnetic field strength in MHz or T.
        nucleus (str, optional): The nucleus type.

    Example:
        reader = SimpReader('spe_file', 'spe', b0='9.4T', nucleus='13C')
    """
    def __init__(self, filename, format, b0=None, nucleus=None):
        self.filename = filename
        self.format = format
        self.b0 = b0
        self.nucleus = nucleus
        self._read_file()

    def _read_file(self):
        if self.format == 'spe':
            self._read_spe()
        elif self.format == 'fid':
            self._read_fid()
        elif self.format == 'xreim':
            self._read_xreim()
        else:
            raise ValueError('Invalid format. Supported formats are spe, fid, and xreim.')
        
    def _read_spe(self):
        """
        This method reads NMR data from a SIMPSON SPE file.
        """
        if self.b0 is None and self.nucleus is None:      
            with open(self.filename) as f:
                data_sec = False
                real = []
                imag = []
                for line in f:
                    if line.startswith('NP'):
                        np_value = float(line.split('=')[1])
                    elif line.startswith('SW'):
                        sw = float(line.split('=')[1])
                    elif line.startswith('DATA'):
                        data_sec = True
                    elif data_sec and line.startswith('END'):
                        break
                    elif data_sec:
                        a, b = map(float, line.split())
                        real.append(a)
                        imag.append(b)
                hz = np.linspace(-int(sw) / 2, int(sw) / 2, int(np_value))
                real = np.array(real)
                imag = np.array(imag)
                hz = np.array(hz)
                self.data = {'real': real, 'imag': imag, 'np': np_value, 'sw': sw, 'hz': hz}
        elif self.b0 is not None and self.nucleus is None or self.b0 is None and self.nucleus is not None:
            raise ValueError('Both B0 and nucleus must be specified.')
        else:
            dir = os.path.dirname(os.path.realpath(__file__))
            isotope_file = os.path.join(dir, 'isotope_data.json')
            with open(self.filename) as f:
                data_sec = False
                real = []
                imag = []
                for line in f:
                    if line.startswith('NP'):
                        np_value = float(line.split('=')[1])
                    elif line.startswith('SW'):
                        sw = float(line.split('=')[1])
                    elif line.startswith('DATA'):
                        data_sec = True
                    elif data_sec and line.startswith('END'):
                        break
                    elif data_sec:
                        a, b = map(float, line.split())
                        real.append(a)
                        imag.append(b)
                hz = np.linspace(-int(sw) / 2, int(sw) / 2, int(np_value))
                real = np.array(real)
                imag = np.array(imag)
                hz = np.array(hz)
                
                isotope = int(''.join(filter(str.isdigit, self.nucleus)))
                element = ''.join(filter(str.isalpha, self.nucleus)).upper()

                b0_unit = ''.join(filter(str.isalpha, self.b0)).lower()

                with open(isotope_file) as f:
                    data = json.load(f)
                    gamma = data[element][str(isotope)]['Gamma']

                if b0_unit == 't':
                    b0 = int(''.join(filter(str.isdigit, self.b0)))
                    ppm = hz / (b0 * gamma)
                elif b0_unit == 'mhz':  
                    b0 = int(''.join(filter(str.isdigit, self.b0)))
                    gamma_h = data['H']['1']['Gamma']
                    ppm = hz / (b0/(gamma_h/gamma))
                else:
                    raise ValueError('B0 unit must be T or MHz.')
                self.data = {'real': real, 'imag': imag, 'np': np_value, 'sw': sw, 'hz': hz, 'ppm': ppm}

    def _read_fid(self):
        """
        This method reads NMR data from a SIMPSON FID file.
        """
        with open(self.filename) as f:
            data_sec = False
            real = []
            imag = []
            for line in f:
                if line.startswith('NP'):
                    np_value = float(line.split('=')[1])
                elif line.startswith('SW'):
                    sw = float(line.split('=')[1])
                elif line.startswith('DATA'):
                    data_sec = True
                elif data_sec and line.startswith('END'):
                    break
                elif data_sec:
                    a, b = map(float, line.split())
                    real.append(a)
                    imag.append(b)
            dt = 1.0 / sw
            time = np.linspace(0, np_value*dt, int(np_value))
            real = np.array(real)
            imag = np.array(imag)
            time = np.array(time)*10e3
            self.data = {'real': real, 'imag': imag, 'np': np_value, 'sw': sw, 'time': time}

    def _read_xreim(self):
        """
        This method reads NMR data from a SIMPSON saved with -xreim option.
        """
        with open(self.filename) as f:
            time = []
            real = []
            imag = []
            for line in f:
                time.append(float(line.split()[0]))
                real.append(float(line.split()[1]))
                imag.append(float(line.split()[2]))
            time = np.array(time)
            real = np.array(real)
            imag = np.array(imag)
            self.data = {'time': time, 'real': real, 'imag': imag}

    def to_spe(self):
        """
        Converts FID data to spectrum (SPE).

        Raises:
            ValueError: If the format is not FID.
        
        Returns:
            SimpReader: A new SimpReader instance with SPE format data.

        Example:
            spectrum = reader.to_spe()
        """
        if self.format != 'fid':
            raise ValueError('Only FID format can be converted to SPE.')

        spectrum = copy.deepcopy(self)

        npoints = spectrum.data['np']
        sw = spectrum.data['sw']
        raw_signal = spectrum.data['real'] + 1j * spectrum.data['imag']
        signal = np.fft.fftshift(np.fft.fft(raw_signal))
        real = np.real(signal)
        imag = np.imag(signal)
        hz = np.linspace(-sw/2, sw/2, int(npoints))
        spectrum.data = {'real': real, 'imag': imag, 'np': npoints, 'sw': sw, 'hz': hz}

        if spectrum.b0 is not None and spectrum.nucleus is not None:
            dir = os.path.dirname(os.path.realpath(__file__))
            isotope_file = os.path.join(dir, 'isotope_data.json')

            isotope = int(''.join(filter(str.isdigit, spectrum.nucleus)))
            element = ''.join(filter(str.isalpha, spectrum.nucleus)).upper()
            b0_unit = ''.join(filter(str.isalpha, spectrum.b0)).lower()

            with open(isotope_file) as f:
                data = json.load(f)
                if element in data:
                    gamma = data[element][str(isotope)]['Gamma']
                else:
                    raise ValueError('Nucleus not found.')

            if b0_unit == 't':
                b0 = int(''.join(filter(str.isdigit, spectrum.b0)))
                ppm = hz / (b0 * gamma)
            elif b0_unit == 'mhz':  
                # Convert from 1H MHz to MHz of the nucleus
                b0 = int(''.join(filter(str.isdigit, spectrum.b0)))
                gamma_h = data['H']['1']['Gamma']
                ppm = hz / (b0/(gamma_h/gamma))
            else:
                raise ValueError('B0 unit must be T or MHz.')

            spectrum.data['ppm'] = ppm

        spectrum.format = 'spe'

        return spectrum
    
    def to_fid(self):
        """
        Converts spectrum (SPE) data to FID.
    
        Raises:
            ValueError: If the format is not SPE.
        
        Returns:
            SimpReader: A new SimpReader instance with FID format data.
        """
        if self.format != 'spe':
            raise ValueError('Only SPE format can be converted to FID.')
    
        fid = copy.deepcopy(self)
    
        npoints = fid.data['np']
        sw = fid.data['sw']
        hz = fid.data['hz']
        signal = fid.data['real'] + 1j * fid.data['imag']
        signal = np.fft.ifft(np.fft.ifftshift(signal))
        real = np.real(signal)
        imag = np.imag(signal)
        dt = 1.0 / sw
        time = np.linspace(0, npoints*dt, int(npoints)) * 10e3  # Match _read_fid scaling
        fid.data = {'real': real, 'imag': imag, 'np': npoints, 'sw': sw, 'time': time}
    
        fid.format = 'fid'
    
        return fid
