# This module contains functions for analyzing NMR data from SIMPSON simulations.
#
# It provides functions for reading NMR data from SPE and FID files, performing
# Fourier transforms, and plotting NMR spectra.

import numpy as np
import json
import os
import sys
from soprano.calculate.nmr.simpson import write_spinsys
from soprano.data.nmr import _get_nmr_data
from soprano.properties.nmr import (
    MSIsotropy,
    MSReducedAnisotropy,
    MSAsymmetry,
    MSQuaternion,
    EFGQuadrupolarConstant,
    EFGAsymmetry,
    EFGQuaternion,
    DipolarCoupling,
)
import warnings
import copy

class SimpReader:
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
                element = ''.join(filter(str.isalpha, self.nucleus))

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
            element = ''.join(filter(str.isalpha, spectrum.nucleus))
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



# Write Simpson functions

_spinsys_template = """spinsys {{
{header}
{ms}
{efg}
{dipolar}
}}
"""

_header_template = """
channels {channels}
nuclei {nuclei}
"""
    
def write_spinsys(
    s,
    isotope_list=None,
    use_ms=False,
    ms_iso=False,
    q_order=0,
    dip_sel=None,
    self_coupling=False,
    path=None,
    ref={},
    grad=-1.0,
    avg_nmr=False,
):
    """
    Write a .spinsys input file for use with SIMPSON, given the details of a
    system. This is meant to be a low-level function, used by other higher
    level interfaces in NMRCalculator.

    | Args:
    |   s   (ase.Atoms):  atomic structure containing the desired spins. 
    |                     All atoms will be included - if that is not the 
    |                     desired result, this should be accomplished by 
    |                     making this a subset of the full structure.
    |   isotope_list ([int]): list of isotopes for each element in the system.
    |                         If left to None, default NMR-active isotopes
    |                         will be used.
    |   use_ms (bool): if True, include shift interactions from magnetic
    |                  shieldings.
    |   ms_iso (bool): if True, all magnetic shieldings will be made
    |                  isotropic.
    |   q_order(int): if greater than 0, include quadrupolar interactions from
    |                   Electric Field Gradients at the given order (1 or 2).
    |   dip_sel (AtomSelection): if not None, include dipolar couplings
    |                            between atoms belonging to this set.
    |   self_coupling (bool): if True, include dipolar couplings between
    |                         atoms and its image.
    |   path (str): path to save the newly created file to. If not provided,
    |               the contents will be simply returned as a string.
    |   ref (dict): dictionary of reference values for the calculation. This
    |               is used to convert from raw shielding values to chemical
    |               shifts. The dictionary should be in the form
    |               {element: value}, where value is the reference shielding
    |               for that element in ppm.
    |   grad (float|dict|list): gradient to use when converting from raw
    |                           shielding values to chemical shifts. If a
    |                           float is provided, it will be used for all
    |                           elements. If a dictionary is provided, it
    |                           should be in the form {element: value}, where
    |                           value is the gradient for that element. If a
    |                           list is provided, it should be have one value
    |                           per site. Defaults to a gradient of -1.0 for
    |                           all elements.
    |   avg_nmr (bool): if True and a trajectory is provided, NMR parameters
    |                   (shifts, cq, assymetries, dipolar couplings) will be 
    |                   averaged over the trajectory. This is useful for
    |                   MD simulations.
    | Returns:
    |   file_contents (str): spinsys file in string format. Only returned if
    |                        no save path is provided.

    """
    # Start by creating a proper isotope_list
    nmr_data = _get_nmr_data()

    nuclei = s.get_chemical_symbols()

    if isotope_list is None:
        isotope_list = [int(nmr_data[n]["iso"]) for n in nuclei]

    nuclei = [str(i) + n for i, n in zip(isotope_list, nuclei)]

    # Build header
    header = _header_template.format(
        channels=" ".join(set(nuclei)), nuclei=" ".join(nuclei)
    )

    # Build MS block
    ms_block = ""
    if use_ms:
        if not ref:
            warnings.warn(
                "No reference values provided for the calculation of "
                "chemical shifts. Assuming all zero."
                "To avoid this warning, provide a dictionary of the form "
                "{element: value}, where value is the reference shielding "
                "for that element in ppm."
            )
        if avg_nmr == False:
            msiso = MSIsotropy.get(s, ref=ref, grad=grad)
            if not ms_iso:
                msaniso = MSReducedAnisotropy.get(s)
                msasymm = MSAsymmetry.get(s)
                eulangs = (
                    np.array([q.euler_angles() for q in MSQuaternion.get(s)]) * 180 / np.pi
                )
            else:
                msaniso = np.zeros(len(s))
                msasymm = np.zeros(len(s))
                eulangs = np.zeros((len(s), 3))

        else:
            # Check if structure is a trajectory (list)
            if isinstance(s, list):
                s_avg = s[0].copy()
                all_ms_tensor = []
                for step in s:
                    all_ms_tensor.append(step.get_array('ms'))

                avg_ms_tensor = np.mean(all_ms_tensor, axis=0)
                s_avg.set_array('ms', avg_ms_tensor)
                msiso = MSIsotropy.get(s_avg, ref=ref, grad=grad)

                if not ms_iso:
                    msaniso = MSReducedAnisotropy.get(s_avg)
                    msasymm = MSAsymmetry.get(s_avg)
                    eulangs = (
                        np.array([q.euler_angles() for q in MSQuaternion.get(s_avg)]) * 180 / np.pi
                    )
                else:
                    msaniso = np.zeros(len(s_avg))
                    msasymm = np.zeros(len(s_avg))
                    eulangs = np.zeros((len(s_avg), 3))

            else:
                raise ValueError('The file is not a ASE trajectory')
            
    else:
        msiso = msaniso = msasymm = eulangs = None

    for i, ms in enumerate(msiso):
        ms_block += ("shift {0} {1}p {2}p " "{3} {4} {5} {6}\n").format(
            i + 1, ms, msaniso[i], msasymm[i], *eulangs[i]
            )

    # Build EFG block
    efg_block = ""
    if q_order > 0:
        if q_order > 2:
            raise ValueError("Invalid quadrupolar order")
        if avg_nmr == False:
            Cq = EFGQuadrupolarConstant(isotope_list=isotope_list)(s)
            eta_q = EFGAsymmetry.get(s)
            eulangs = (
                np.array([q.euler_angles() for q in EFGQuaternion.get(s)]) * 180 / np.pi
            )
        else:

            if isinstance(s, list):
                s_avg = s[0].copy()
                all_efg_tensor = []
                for step in s:
                    all_efg_tensor.append(step.get_array('efg'))

                avg_efg_tensor = np.mean(all_efg_tensor, axis=0)
                s_avg.set_array('efg', avg_efg_tensor)
                Cq = EFGQuadrupolarConstant(isotope_list=isotope_list)(s_avg)
                eta_q = EFGAsymmetry.get(s_avg)
                eulangs = (
                    np.array([q.euler_angles() for q in EFGQuaternion.get(s_avg)]) * 180 / np.pi
                )
            else:
                raise ValueError('The file is not a ASE trajectory')
            
        for i, cq in enumerate(Cq):
            if cq == 0:
                continue
            efg_block += ("quadrupole {0} {1} {2} {3}" " {4} {5} {6}\n").format(
                i + 1, q_order, cq, eta_q[i], *eulangs[i]
            )

    # Build dipolar block
    dip_block = ""
    if dip_sel is not None and len(dip_sel) > 1:
        if avg_nmr == False:
            dip = DipolarCoupling(sel_i=dip_sel, isotope_list=isotope_list, self_coupling=self_coupling)(s)
        else:
            if isinstance(s, list):
                s_avg = s[0].copy()
                all_dip = {}

                for step in s:
                    (i, j), (d, v) = DipolarCoupling(sel_i=dip_sel, isotope_list=isotope_list, self_coupling=self_coupling)(step)
                    all_dip[(i, j)] = (d, v)

                dip = {}
                for key in all_dip.keys():
                    dip[key] = (np.mean([all_dip[key][0]]), np.mean([all_dip[key][1]], axis=0))

            else:
                raise ValueError('The file is not a ASE trajectory')

        for (i, j), (d, v) in dip.items():
            a, b = (
                (np.array([np.arccos(-v[2]), np.arctan2(-v[1], -v[0])]) % (2 * np.pi))
                * 180
                / np.pi
            )
            dip_block += ("dipole {0} {1} {2} {3}" " {4} 0.0\n").format(
                i + 1, j + 1, d * 2 * np.pi, a, b
            )

    out_file = _spinsys_template.format(
        header=header, ms=ms_block, efg=efg_block, dipolar=dip_block
    )

    if path is None:
        return out_file
    else:
        with open(path, "w") as of:
            of.write(out_file)
