# SimPYson: A Pythonic Interface for SIMPSON

[![DOI](https://zenodo.org/badge/813117518.svg)](https://doi.org/10.5281/zenodo.14041918)

SimPYson is a Python package designed to simplify the use of [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson), a powerful code to simulate solid-state NMR experiments. Born out of the need to streamline the process, SimPYson makes it easier to prepare input files from DFT calculations, run simulations, and analyze results from SIMPSON all within Python.

## Features

- **Convert DFT data to SIMPSON input files**: Prepare SIMPSON input files from DFT data (CASTEP, Quantum Espresso, VASP) using ASE and Soprano.

- **Run SIMPSON simulations from Python**: Use the `SimpCalc` calculator or the `simulate_spectrum()` convenience function to generate, run, and read SIMPSON simulations without leaving Python.

- **Read SIMPSON output files**: Load and manipulate NMR data from SIMPSON `.spe`, `.fid`, `.xreim`, and `.csdf` files directly in Python for further analysis and visualization.

- **Graphical User Interface**: Type `simpyson gui` in your terminal to launch the SimPYson GUI (PyQt5 + Plotly) and manipulate data without any coding.

- **Pulse sequence templates**: Ready-made templates for common experiments: no-pulse, 90-degree pulse, and CPMAS. Custom pulse sequences are also supported.

## Quick Start

```python
from simpyson.io import read_simp

# Read a SIMPSON spectrum file
data = read_simp("spectrum.spe", b0="400MHz", nucleus="13C")

# Access frequency-domain data
print(data.spe['hz'])   # Hz axis
print(data.ppm['ppm'])  # ppm axis (requires b0 and nucleus)
```

```python
from simpyson.calculator import simulate_spectrum

# Simulate a spectrum from a spin system string
spinsys = """
channels 13C
nuclei 13C
shift 1 10p 0 0 0 0 0
"""
result = simulate_spectrum(spinsys, proton_frequency=400e6)
```

## Documentation

Full documentation with tutorials is available at [carlosbornes.github.io/simpyson](https://carlosbornes.github.io/simpyson/).

## Installation

```bash
pip install git+https://github.com/nuts-org/simpyson.git
```

Requires Python >= 3.10 and a working [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) installation for running simulations.

## Planned Features

- Expand the number of pulse sequence templates for more complex NMR experiments.
- Improve support for additional DFT codes -- suggestions are welcome.
