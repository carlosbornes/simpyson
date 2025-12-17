from __future__ import annotations

import json
import os

import numpy as np
from soprano.calculate.nmr.simpson import _header_template, _spinsys_template
from soprano.data.nmr import _get_isotope_list
from soprano.properties.nmr.dipolar import DipolarCoupling
from soprano.selection import AtomSelection


def get_gamma(nucleus, isotope_file=None):
    """
    Get gyromagnetic ratio for a given nucleus.
    
    Args:
        nucleus (str): Nucleus type (e.g., '1H' or '13C')
        isotope_file (str, optional): Path to isotope data file. If None, uses default.

    Returns:
        float: Gyromagnetic ratio in Hz/T
    """
    if isotope_file is None:
        dir = os.path.dirname(os.path.realpath(__file__))
        isotope_file = os.path.join(dir, 'isotope_data.json')

    isotope = int(''.join(filter(str.isdigit, nucleus)))
    element = ''.join(filter(str.isalpha, nucleus)).capitalize()

    with open(isotope_file) as f:
        data = json.load(f)
        if element in data and str(isotope) in data[element]:
            gamma = data[element][str(isotope)]['Gamma']
        else:
            raise ValueError(f'Nucleus {nucleus} not found in isotope data.')

    return gamma

def get_spin(nucleus, isotope_file=None):
    """
    Get spin quantum number for a given nucleus.
    
    Args:
        nucleus (str): Nucleus type (e.g., '1H' or '13C')
        isotope_file (str, optional): Path to isotope data file. If None, uses default.

    Returns:
        float: Spin quantum number (e.g. 0.5, 1.0, 1.5)
    """
    if isotope_file is None:
        dir = os.path.dirname(os.path.realpath(__file__))
        isotope_file = os.path.join(dir, 'isotope_data.json')

    isotope = int(''.join(filter(str.isdigit, nucleus)))
    element = ''.join(filter(str.isalpha, nucleus)).capitalize()

    with open(isotope_file) as f:
        data = json.load(f)
        if element in data and str(isotope) in data[element]:
            spin_str = data[element][str(isotope)]['Spin']
            # Spin is usually a string like "1/2" or "3/2" or integer "1"
            if '/' in str(spin_str):
                num, den = spin_str.split('/')
                return float(num) / float(den)
            else:
                return float(spin_str)
        else:
            raise ValueError(f'Nucleus {nucleus} not found in isotope data.')

def get_larmor_freq(b0, nucleus, isotope_file=None):
    """
    Convert Hz values to ppm values.
    
    Args:
        b0 (str): Magnetic field strength (e.g., '400MHz' or '9.4T')
        nucleus (str): Nucleus type (e.g., '1H' or '13C')
        isotope_file (str, optional): Path to isotope data file. If None, uses default.
    
    Returns:
        float: Larmor Frequency
        
    Raises:
        ValueError: If B0 unit is invalid or nucleus not found
    """

    if isotope_file is None:
        dir = os.path.dirname(os.path.realpath(__file__))
        isotope_file = os.path.join(dir, 'isotope_data.json')

    isotope = int(''.join(filter(str.isdigit, nucleus)))
    element = ''.join(filter(str.isalpha, nucleus)).capitalize()
    b0_unit = ''.join(filter(str.isalpha, b0)).lower()

    with open(isotope_file) as f:
        data = json.load(f)
        if element in data and str(isotope) in data[element]:
            gamma = get_gamma(nucleus, isotope_file=isotope_file)
        else:
            raise ValueError(f'Nucleus {nucleus} not found in isotope data.')

    if b0_unit == 't':
        b0_value = float(''.join(filter(lambda x: x.isdigit() or x == '.', b0)))
        larmor_freq = gamma * 1e7 * b0_value / (2 * np.pi * 1e6)
    elif b0_unit == 'mhz':
        b0_value = float(''.join(filter(lambda x: x.isdigit() or x == '.', b0)))
        gamma_h = get_gamma('1H', isotope_file=isotope_file)
        b0_value_T = 2 * np.pi * b0_value * 1e6 / (gamma_h * 1e7)
        larmor_freq = gamma * 1e7 * b0_value_T / (2 * np.pi * 1e6)
    else:
        raise ValueError('B0 unit must be T or MHz.')

    return larmor_freq

def add_spectra(spectra_list, b0=None, nucleus=None):
    """Combine multiple Simpy objects into a single spectrum."""
    if not spectra_list:
        return None

    result = spectra_list[0].copy()

    if b0:
        result.b0 = b0
    if nucleus:
        result.nucleus = nucleus

    for spectrum in spectra_list[1:]:
        result._spe_data['real'] += spectrum.spe['real']
        result._spe_data['imag'] += spectrum.spe['imag']

    if result.b0 and result.nucleus:
        result._compute_ppm()

    return result

# Code adapted from soprano to construct
# simple spin systems by hand
def simple_spinsys(
    atoms,
    isotopes,
    iso_ms=None,
    aniso_ms=None,
    eta_ms=None,
    euler_ms=None,
    cq=None,
    eta_q=None,
    q_order=None,
    euler_q=None,
    get_dipolar=False,
    dip_sel=None,
    obs_nuc=None,
):
    """
    Generates a SIMPSON .spinsys file string from directly provided NMR parameters.

    Args:
        atoms (ase.Atoms): Atoms to be considered.
        isotopes (dict): A dictionary mapping element to isotopes.
        iso_ms (list, optional): List of isotropic chemical shifts (ppm).
        aniso_ms (list, optional): List of shielding anisotropies (ppm).
        eta_ms (list, optional): List of shielding asymmetries.
        euler_ms (list, optional): List of Euler angles (alpha, beta, gamma) in degrees for shielding.
        cq (list, optional): List of quadrupolar coupling constants (Hz).
        eta_q (list, optional): List of quadrupolar asymmetry parameters.
        q_order (list, optional): List of quadrupolar orders (<=2).
        euler_q (list, optional): List of Euler angles (alpha, beta, gamma) in degrees for EFG.
        get_dipolar (bool, optional): If True, calculate and include dipolar couplings.
        dip_sel (AtomSelection, optional): Selection of atoms for dipolar couplings. Defaults to all.
        obs_nuc (str, optional): The nucleus to be observed.

    Returns:
        str: The contents of the .spinsys file as a string.
    """

    # Header Block
    symbols = atoms.get_chemical_symbols()
    isotope_list = _get_isotope_list(symbols, isotopes=isotopes)
    nuclei = [f"{iso}{el}" for el, iso in zip(symbols, isotope_list)]

    if obs_nuc and obs_nuc in nuclei:
        channels = [obs_nuc] + [n for n in sorted(set(nuclei)) if n != obs_nuc]
    else:
        channels = sorted(set(nuclei))

    header = _header_template.format(
        channels=" ".join(channels), nuclei=" ".join(nuclei)
    )

    # Magnetic shielding Block
    ms_block = ""
    if iso_ms is not None:
        num_atoms = len(atoms)
        aniso_ms_vals = aniso_ms if aniso_ms is not None else np.zeros(num_atoms)
        eta_ms_vals = eta_ms if eta_ms is not None else np.zeros(num_atoms)
        euler_ms_vals = euler_ms if euler_ms is not None else np.zeros((num_atoms, 3))

        for i in range(num_atoms):
            ms_block += "shift {0} {1}p {2}p {3} {4} {5} {6}\n".format(
                i + 1,
                iso_ms[i],
                aniso_ms_vals[i],
                eta_ms_vals[i],
                *euler_ms_vals[i],
            )

    # Quadrupolar Block
    efg_block = ""
    if cq is not None:
        num_atoms = len(atoms)
        eta_q_vals = eta_q if eta_q is not None else np.zeros(num_atoms)
        euler_q_vals = euler_q if euler_q is not None else np.zeros((num_atoms, 3))

        for i in range(num_atoms):
            if cq[i] != 0:
                if q_order is not None:
                    if q_order[i] > 2:
                        raise ValueError(
                            f"Quadrupolar order must be 2 or less, got {q_order[i]} for atom {i+1}"
                        )
                    efg_block += "quadrupole {0} {1} {2} {3} {4} {5} {6}\n".format(
                        i + 1, q_order[i], cq[i], eta_q_vals[i], *euler_q_vals[i]
                    )
                else:
                    efg_block += "quadrupole {0} 2 {1} {2} {3} {4} {5}\n".format(
                        i + 1, cq[i], eta_q_vals[i], *euler_q_vals[i]
                    )

    # Dipolar Block
    dip_block = ""
    if get_dipolar:
        if dip_sel is None:
            dip_sel = AtomSelection.all(atoms)

        if len(dip_sel) > 1:
            dip_couplings = DipolarCoupling.get(
                atoms, sel_i=dip_sel, isotope_list=isotope_list
            )
            for (i, j), (d, v) in dip_couplings.items():
                # Convert units
                d_rad_s = d * 2 * np.pi

                if np.allclose(v, [0, 0, 1]):
                    beta, alpha = 0.0, 0.0
                elif np.allclose(v, [0, 0, -1]):
                    beta, alpha = 180.0, 0.0
                else:
                    beta = np.arccos(v[2]) * 180 / np.pi
                    alpha = np.arctan2(v[1], v[0]) * 180 / np.pi

                dip_block += f"dipole {i + 1} {j + 1} {d_rad_s} {beta} {alpha} 0.0\n"

    return _spinsys_template.format(
        header=header, ms=ms_block, efg=efg_block, dipolar=dip_block
    )
