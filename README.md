# SimPYson: A Pythonic Interface for SIMPSON

SimPYson tries to combine Simpson NMR simulation and the versatility of Python, providing a seamless and efficient way to interact with SIMPSON within your Python workflows. Currently, the code is in a embryonic stage, but we envision it will evolve into a comprehensive toolset that facilitates the following functionalities:

# Planned Features:

- DFT Data Parsing: Seamlessly parse and extract relevant data from DFT calculations into SIMPSON.

- Spinsys Preparation: Automate the process of preparing SIMPSON spinsys files from extracted DFT data.

- SIMPSON Input File Generation: Generate SIMPSON input files from prepared spinsys files, saving time and effort in setting up your simulations.

- Output File Analysis: Extract and analyze NMR data from SIMPSON output files, facilitating further data processing and exploration.

# Current Implementation:

- Read data from SIMPSON spectra: Read and manipulate NMR data from SIMPSON .spe and .fid files directly into Python, enabling further analysis and processing.

We are actively working on developing these features and plan to release the package soon. In the meantime, you can utilize the existing ``read_spe`` and ``read_fid`` functions to read and analyze NMR data from SIMPSON spectra. Here is an example code snippet:

```python
from simpyson.io import read_spe
import matplotlib.pyplot as plt

# Read data from the SPE file
spectrum = read_spe('ouput.spe')
fid = read_fid('output.fid')

# If B0 and nucleus are provided Simpyson calculates ppm scale automatically
spectrum_ppm = read_spe('ouput.spe', b0='800MHz', nucleus='13C')

# Plot the data
plt.plot(data['hz'], data['real'])
plt.xlabel('Frequency (Hz)')

plt.plot(fid['time'], data['real'])
plt.xlabel('Time (ms)')

plt.plot(data['ppm'], data['real'])
plt.xlabel('$^{13}$C (ppm)')
```
