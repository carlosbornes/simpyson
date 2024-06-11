# This module contains functions for converting shielding and efg tensor
#
# It provides functions for reading NMR data from ASE ms and efg arrays,
# convert to distinct tensor convention and calculate isotropic shielding,
# anisotropy, cq, etc..

from ase.io import read
from ase.quaternions import Quaternion
import numpy as np
import fractions
import scipy.constants as const
import os
import json
import warnings
from ase.geometry import find_mic

# Magnetic shielding tensor functions

def get_symms(filename):
    ms = filename.get_array('ms')
    sym_ms = np.zeros((len(ms), 3, 3))
    for i in range(len(ms)):
        sym_ms[i] = 0.5 * (ms[i] + ms[i].T)
    return sym_ms
    
def get_diag(filename):
    sym_ms = get_symms(filename)
    eigvals, eigvecs = np.linalg.eigh(sym_ms, )
    return eigvals, eigvecs

def get_iso(filename, ref=None, grad=None):
    diag_ms = get_diag(filename)[0]
    iso = []
    if ref is not None and grad is not None:
        for i in diag_ms:
            iso.append(np.mean(i) * grad + ref)
    elif ref is not None and grad is None or ref is None and grad is not None:
        raise ValueError('Both ref and grad must be specified.')
    else:
        for i in diag_ms:
            iso.append(np.mean(i))
    return iso

def sort_diag_haeb(filename):
    diag = get_diag(filename)[0]
    iso = get_iso(filename)
    diff = np.abs(diag - iso[:, None])
    idx = np.argsort(diff, axis=1)
    sorted_diag = np.take_along_axis(diag, idx, axis=1)
    return sorted_diag

def sort_diag(filename):
    diag = get_diag(filename)[0]
    sorted_diag = np.sort(diag, axis=1)
    return sorted_diag
    
def get_aniso(filename):
    iso = get_iso(filename)[0]
    sorted_diag = sort_diag_haeb(filename)
    aniso = []
    for i in range(len(iso)):
        aniso.append((sorted_diag[i][-1]-(sorted_diag[i][1]+sorted_diag[i][0])/2))
    return aniso

def get_asym(filename):
    iso = get_iso(filename)[0]
    sorted_diag = sort_diag_haeb(filename)
    asym = []
    for i in range(len(iso)):
        asym.append((sorted_diag[i][0]-sorted_diag[i][1])/(sorted_diag[i][-1]-iso[i]))
    return asym

def get_span(filename):
    sorted_diag = sort_diag(filename)
    span = []
    for i in range(len(sorted_diag)):
        span.append(sorted_diag[i][-1]-sorted_diag[i][0])
    return span

def get_skew(filename):
    iso = get_iso(filename)
    span = get_span(filename)
    sorted_diag = sort_diag(filename)
    skew = []
    for i in range(len(iso)):
        skew.append((3*(sorted_diag[i][1]-iso[i])/(span[i])))
    return skew


# EFG tensor functions

def get_efg_diag(filename):
    efg = filename.get_array('efg')
    eigvals, eigvecs = np.linalg.eigh(efg)

    return eigvals, eigvecs

def sort_efg_diag(filename):
    diag = get_efg_diag(filename)[0]
    abs_eigval = np.abs(diag)
    idx = np.argsort(abs_eigval, axis=1)
    sorted_diag = np.take_along_axis(diag, idx, axis=1)
    return sorted_diag

def efg_vzz(filename):
    sorted_diag = sort_efg_diag(filename)
    vzz = []
    for i in range(len(sorted_diag)):
        vzz.append(sorted_diag[i][2])
    return vzz

def get_cq(filename, isotopes=None):
    elements = filename.get_chemical_symbols()
    dir = os.path.dirname(os.path.realpath(__file__))
    isotope_file = os.path.join(dir, 'isotope_data.json')
    with open(isotope_file, 'r') as f:
        isotope_data = json.load(f)

    Q = []
    if isotopes is None:
        isotopes = {}
    for element in elements:
        if element in isotopes:
            if element not in isotope_data or isotopes[element] not in isotope_data[element]:
                raise ValueError(f"The isotope {isotopes[element]} for element {element} does not exist in the isotope data.")
            Q.append(isotope_data[element][isotopes[element]]['QMoment'])
        else:
            if len(isotope_data[element]) > 1:
                valid_isotopes = {k: v for k, v in isotope_data[element].items() if float(fractions.Fraction(v['Spin'])) >= 1}
                if valid_isotopes:
                    max_abundance = max(valid_isotopes.items(), key=lambda x: x[1]['NatAbudance'])
                    isotopes[element] = max_abundance[0]
                    Q.append(max_abundance[1]['QMoment'])
                else:
                    raise ValueError(f"No isotopes with Spin >= 1 for element {element}.")
            else:
                isotopes[element] = list(isotope_data[element].items())[0][0]
                Q.append(list(isotope_data[element].items())[0][1]['QMoment'])

    vzz = efg_vzz(filename)
    cq = const.e * np.array(vzz) * np.array(Q) * const.physical_constants["atomic unit of electric field gradient"][0]*1e-30/ const.h
    return cq

def get_eta(filename):
    sorted_diag = sort_efg_diag(filename)
    eta = []
    for i in range(len(sorted_diag)):
        eta.append((sorted_diag[i][0]-sorted_diag[i][1])/(sorted_diag[i][2]))
    return eta

# Quaternion functions to express the orientation of tensors with respect to the cartesian axes

def ms_quaternion(filename):
    eigvecs = get_diag(filename)[1]
    eigvecs = np.array(eigvecs) * np.linalg.det(eigvecs)[:, None, None]
    ms_quat = []
    for vec in eigvecs:
        ms_quat.append(Quaternion.from_matrix(vec.T).q)
    return ms_quat

def efg_quaternion(filename):
    eigvecs = get_efg_diag(filename)[1] * np.linalg.det(get_efg_diag(filename)[1])[:, None, None]
    efg_quat = []
    for vec in eigvecs:
        efg_quat.append(Quaternion.from_matrix(vec.T).q)
    return efg_quat

# Dipolar coupling functions

def get_dipcoup(filename, atom_i=None, atom_j=None, isotopes=None, self_coupling=False):
    if atom_i is None and atom_j is None:
        atom_i = range(len(filename))
    elif atom_i is None and atom_j is not None:
        atom_i = atom_j

    elements = filename.get_chemical_symbols()
    dir = os.path.dirname(os.path.realpath(__file__))
    isotope_file = os.path.join(dir, 'isotope_data.json')
    with open(isotope_file, 'r') as f:
        isotope_data = json.load(f)

    gamma = []
    if isotopes is None:
        isotopes = {}
    for element in elements:
        if element in isotopes:
            if element not in isotope_data or isotopes[element] not in isotope_data[element]:
                raise ValueError(f"The isotope {isotopes[element]} for element {element} does not exist in the isotope data.")
            gamma.append(isotope_data[element][isotopes[element]]['Gamma'])
        else:
            if len(isotope_data[element]) > 1:
                valid_isotopes = {k: v for k, v in isotope_data[element].items()}
                if valid_isotopes:
                    max_abundance = max(valid_isotopes.items(), key=lambda x: x[1]['NatAbudance'])
                    isotopes[element] = max_abundance[0]
                    gamma.append(max_abundance[1]['Gamma'])
            else:
                isotopes[element] = list(isotope_data[element].items())[0][0]
                gamma.append(list(isotope_data[element].items())[0][1]['Gamma'])

    dipcoup = {}
    for i in atom_i:
        for j in atom_j:
            if i != j:
                vec, r = find_mic(filename.get_positions()[i] - filename.get_positions()[j], filename.get_cell(), pbc=True)
                coupling = -const.physical_constants["vacuum mag. permeability"][0]*const.hbar*gamma[i]*1e7*gamma[j]*1e7/(8*np.pi**2*(r*1e-10)**3)
                dipcoup.update({(i, j): {'dip coup': coupling, 'vector': vec}})
            elif i==j and self_coupling:
                supercell = filename.repeat((2, 2, 2))
                pos = filename.get_positions()[i]
                imgs = []
                for k in range(3):
                    img_pos = pos + filename.cell[k]
                    mics = find_mic(img_pos - pos, supercell.cell, pbc=True)
                    imgs.append(mics)

                vec, r = sorted(imgs, key=lambda x: x[-1])[0]
                coupling = -const.physical_constants["vacuum mag. permeability"][0]*const.hbar*gamma[i]*1e7*gamma[j]*1e7/(8*np.pi**2*(r*1e-10)**3)
                dipcoup.update({(i, i): {'dip coup': coupling, 'vector': vec}})
            elif i==j and not self_coupling:
                continue

    return dipcoup

# Convert FID to SPE

def fid_to_spe(fid: dict, b0=None, nucleus=None):
    """
    This function converts a FID to a SPE dictionary

    Args:
        FID (dict): dictionary returned by read_fid.

    Returns:
        dict: A dictionary containing the following keys:
        real (numpy.ndarray): Real part of the NMR data.
        img (numpy.ndarray): Imaginary part of the NMR data.
        np (int): Number of data points.
        sw (float): Spectral width (Hz).
        hz (numpy.ndarray): NMR frequencies (Hz).
        ppm (numpy.ndarray): NMR frequencies (ppm).
    """

    np = fid['np']
    sw = fid['sw']
    raw_signal = fid['real'] + 1j * fid['img']
    signal = np.fft.fftshift(np.fft.fft(raw_signal))
    real = np.real(signal)
    img = np.imag(signal)
    hz = np.linspace(-sw/2, sw/2, np)
    dict = {'real': real, 'img': img, 'np': np, 'sw': sw, 'hz': hz}

    if b0 is not None and nucleus is not None:
        dir = os.path.dirname(os.path.realpath(__file__))
        isotope_file = os.path.join(dir, 'isotope_data.json')

        isotope = int(''.join(filter(str.isdigit, nucleus)))
        element = ''.join(filter(str.isalpha, nucleus))
        b0_unit = ''.join(filter(str.isalpha, b0)).lower()

        # Get the gyromagnetic ratio from isotope_data.json
        with open(isotope_file) as f:
            data = json.load(f)
            if element in gamma:
                gamma = data[element][str(isotope)]['Gamma']
            else:
                raise ValueError('Nucleus not found.')

        if b0_unit == 't':
            b0 = int(''.join(filter(str.isdigit, b0)))
            ppm = hz / (b0 * gamma)
        elif b0_unit == 'mhz':  
            # Convert from 1H MHz to MHz of the nucleus
            b0 = int(''.join(filter(str.isdigit, b0)))
            gamma_h = data['H']['1']['Gamma']
            ppm = hz / (b0/(gamma_h/gamma))
        else:
            raise ValueError('B0 unit must be T or MHz.')
        
        dict['ppm'] = ppm

    return dict