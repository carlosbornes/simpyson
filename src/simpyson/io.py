from __future__ import annotations

import logging
import os
import warnings

import numpy as np

from simpyson.simpy import Simpy

logger = logging.getLogger("simpyson")


def read_simp(
    filename: str,
    format: str | None = None,
    b0: str | None = None,
    nucleus: str | None = None,
) -> Simpy:
    """
    Read SIMPSON NMR data from a file into a unified Simpy object.

    The file format is determined from the extension if ``format`` is not
    given explicitly.

    Parameters
    ----------
    filename : str
        Path to the SIMPSON output file.
    format : str or None
        File format (``'spe'``, ``'fid'``, ``'xreim'``, ``'csdf'``).
        If None, guessed from the file extension.
    b0 : str or None
        Magnetic field strength (e.g., ``'9.4T'``, ``'400MHz'``).
        Needed for ppm conversion.
    nucleus : str or None
        Nucleus type (e.g., ``'1H'``, ``'13C'``).
        Needed for ppm conversion.

    Returns
    -------
    Simpy
        Object containing the loaded data.

    Raises
    ------
    ValueError
        If the file format cannot be determined or is unsupported.
    OSError
        If the file cannot be read or parsed.
    """
    supported_formats = {'spe', 'fid', 'xreim', 'csdf'}

    if format is not None:
        format = format.lower()
        if format not in supported_formats:
            raise ValueError(f"Unsupported format {format}")
    else:
        ext = os.path.splitext(filename)[1].lower()
        if ext == '.spe':
            format = 'spe'
        elif ext == '.fid':
            format = 'fid'
        elif ext == '.xreim':
            format = 'xreim'
        elif ext == '.csdf':
            format = 'csdf'
        else:
            raise ValueError(f"Cannot determine file format of {filename}")

    simpy_data = Simpy(b0=b0, nucleus=nucleus)

    try:
        if format == 'spe':
            read_spe(filename, simpy_data)
        elif format == 'fid':
            read_fid(filename, simpy_data)
        elif format == 'xreim':
            read_xreim(filename, simpy_data)
        elif format == 'csdf':
            read_csdf(filename, simpy_data)
        else:
            raise ValueError(f"Unsupported format {format}")
    except (ValueError, KeyError, IndexError, OSError) as e:
        raise OSError(f"Error reading file {filename} as format {format}: {e!s}") from e
    return simpy_data


def read_spe(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON SPE file.

    Parameters
    ----------
    filename : str
        Path to the ``.spe`` file.
    simpy_data : Simpy
        Object to populate with spectrum data.

    Raises
    ------
    ValueError
        If required header fields (NP, SW) are missing.
    """
    with open(filename) as f:
        data_sec = False
        real: list[float] = []
        imag: list[float] = []
        ref = 0.0
        np_value: float | None = None
        sw: float | None = None
        for line in f:
            if line.startswith('NP'):
                np_value = float(line.split('=')[1])
            elif line.startswith('SW'):
                sw = float(line.split('=')[1])
            elif line.startswith('REF'):
                ref = float(line.split('=')[1])
            elif line.startswith('DATA'):
                data_sec = True
            elif data_sec and line.startswith('END'):
                break
            elif data_sec:
                a, b = map(float, line.split())
                real.append(a)
                imag.append(b)

        if np_value is None or sw is None:
            raise ValueError(
                f"Missing required header fields in {filename}: "
                f"{'NP' if np_value is None else ''}"
                f"{' and ' if np_value is None and sw is None else ''}"
                f"{'SW' if sw is None else ''} not found."
            )

        if ref != 0.0:
            warnings.warn(
                f"Using REF={ref} Hz from SPE header as frequency offset. "
                "This assumes ref equals the SIMPSON offset parameter.",
                stacklevel=2,
            )
        indices = np.arange(np_value)
        hz = sw * (indices / np_value - 0.5) - ref

        simpy_data.from_spe(real, imag, np_value, sw, hz)


def read_fid(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON FID file.

    Parameters
    ----------
    filename : str
        Path to the ``.fid`` file.
    simpy_data : Simpy
        Object to populate with FID data.

    Raises
    ------
    ValueError
        If required header fields (NP, SW) are missing.
    """
    with open(filename) as f:
        data_sec = False
        real: list[float] = []
        imag: list[float] = []
        np_value: float | None = None
        sw: float | None = None
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

        if np_value is None or sw is None:
            raise ValueError(
                f"Missing required header fields in {filename}: "
                f"{'NP' if np_value is None else ''}"
                f"{' and ' if np_value is None and sw is None else ''}"
                f"{'SW' if sw is None else ''} not found."
            )

        # Let from_fid() compute the time axis (avoids duplicating the calculation)
        simpy_data.from_fid(np.array(real), np.array(imag), np_value, sw)


def read_xreim(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON file saved with the ``-xreim`` option.

    Parameters
    ----------
    filename : str
        Path to the ``.xreim`` file.
    simpy_data : Simpy
        Object to populate with xreim data.
    """
    with open(filename) as f:
        time: list[float] = []
        real: list[float] = []
        imag: list[float] = []
        for line in f:
            parts = line.split()
            time.append(float(parts[0]))
            real.append(float(parts[1]))
            imag.append(float(parts[2]))

        simpy_data.from_xreim(np.array(time), np.array(real), np.array(imag))


def read_csdf(filename: str, simpy_data: Simpy) -> None:
    """
    Read NMR data from a SIMPSON CSDF file.

    Requires the ``csdmpy`` package to be installed.

    Parameters
    ----------
    filename : str
        Path to the ``.csdf`` file.
    simpy_data : Simpy
        Object to populate with spectrum data.
    """
    import csdmpy as csdm

    data = csdm.load(filename)
    hz = data.dimensions[0].coordinates.value
    real = data.dependent_variables[0].components[0].real
    imag = data.dependent_variables[0].components[0].imag
    np_value = np.array(len(hz))
    sw = np.abs(hz[-1] - hz[0])

    simpy_data.from_csdf(real, imag, hz, np_value, sw)


def write_simp(spinsys: str, out_name: str, **kwargs) -> object:
    """
    Create a SIMPSON input file with the specified parameters.

    This is a convenience wrapper around ``SimpCalc``. For full control, use
    ``SimpCalc`` directly from ``simpyson.calculator``.

    Parameters
    ----------
    spinsys : str
        Spin system definition (SIMPSON spinsys block or body).
    out_name : str
        Output file name (without extension).
    **kwargs
        Additional parameters forwarded to ``SimpCalc`` (e.g.,
        ``proton_frequency``, ``spin_rate``, ``sw``, ``np``, etc.).

    Returns
    -------
    SimpCalc
        Configured calculator object that can be saved with ``.save(path)``.
    """
    # Lazy import to avoid circular dependency (calculator.py -> io.py)
    from simpyson.calculator import SimpCalc

    return SimpCalc(spinsys=spinsys, out_name=out_name, **kwargs)
