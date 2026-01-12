from __future__ import annotations

from simpyson.templates import (
    CustomPulseSequence,
    PulseSequenceTemplate,
    get_template,
    pulseq_templates,
)


class SimpCalc:
    """
    Class to create SIMPSON simulation input files.
    
    This class handles all four main sections of a SIMPSON input file:
    - spinsys: Spin system definition
    - par: Simulation parameters 
    - pulseq: Pulse sequence
    - main: Processing section
    """

    def __init__(self, spinsys, pulse_sequence=None, **kwargs):
        if spinsys is None:
            raise ValueError("spinsys cannot be None")

        self.spinsys = spinsys
        self.parameters = kwargs
        self.output_config = {}

        output_keys = ['out_name', 'out_format', 'lb', 'zerofill']
        for key in output_keys:
            if key in self.parameters:
                self.output_config[key.replace('out_', '')] = self.parameters.pop(key)

        self.pulse_sequence = self._setup_pulse_sequence(pulse_sequence)

    def __str__(self):
        """Generate the complete SIMPSON input file"""
        sections = []
        sections.append(self.generate_spinsys())
        sections.append(self.generate_par())
        sections.append(self.generate_pulseq())
        sections.append(self.generate_main())
        return "\n".join(sections)


    def _setup_pulse_sequence(self, pulse_sequence):
        """Set up the pulse sequence based on user input."""
        if pulse_sequence is None:
             return None

        if isinstance(pulse_sequence, str):
            if pulse_sequence in pulseq_templates:
                template_class = pulseq_templates[pulse_sequence]
                template_instance = template_class()

                # Get variable parameters
                required_clean = {param.replace('variable_', '')
                                 for param in template_instance.get_required_parameters()}

                # Extract matching parameters from user input
                pulseq_params = {}
                for param in required_clean:
                    if param in self.parameters:
                        pulseq_params[param] = self.parameters[param]
                
                # Pass offset if it exists
                if 'variable_offset' in self.parameters:
                    pulseq_params['offset'] = self.parameters['variable_offset']

                # Create template with extracted parameters
                return get_template(pulse_sequence, **pulseq_params)
            else:
                # Custom string sequence
                return CustomPulseSequence(pulse_sequence, **self.parameters)

        elif isinstance(pulse_sequence, PulseSequenceTemplate):
            return pulse_sequence

        else:
            raise ValueError(f"Invalid pulse_sequence type: {type(pulse_sequence)}. "
                           "Must be a template name (str), a custom string, or " \
                           "a PulseSequenceTemplate object")

    def generate_spinsys(self):
        """
        Generates the spinsys section of the SIMPSON input file.
        It should be added either manually or using Soprano
        
        Returns:
            str: The spinsys section as a string.
        """
        if hasattr(self.spinsys, 'to_simpson'):
            self.spinsys = self.spinsys.to_simpson()

        elif isinstance(self.spinsys, str):
            if self.spinsys.strip().startswith("spinsys"):
                 # Already formatted block
                self.spinsys = self.spinsys
            elif self.spinsys.strip().startswith("channels") or self.spinsys.strip().startswith("nuclei"):
                 # Body of spinsys block
                spinsys_block = "spinsys {\n"
                spinsys_block += f'{self.spinsys}\n'
                spinsys_block += "}\n"
                self.spinsys = spinsys_block
            else:
                # Assume it's just the body if it doesn't start with known keywords but is a string
                # This is a bit risky but allows flexibility
                spinsys_block = "spinsys {\n"
                spinsys_block += f'{self.spinsys}\n'
                spinsys_block += "}\n"
                self.spinsys = spinsys_block
        else:
            raise ValueError(f"spinsys must be a string or a Soprano SpinSystem object. Got {type(self.spinsys)}")

        return str(self.spinsys)

    def generate_par(self):
        """
        Generates the par section of the SIMPSON input file.
        
        Returns:
            str: The par section as a string.
        """

        # Parameters required for every simulation
        required_params = {
            "proton_frequency", "spin_rate", "start_operator", "detect_operator",
            "np", "sw", "method", "crystal_file", "gamma_angles", "verbose"
        }


        missing_params = required_params - set(self.parameters.keys())
        if missing_params:
            raise ValueError(f"Missing required parameters: {', '.join(missing_params)}")

        par_block = "par {\n"


        for param in sorted(required_params):
            if param in self.parameters:
                value = self.parameters[param]
                par_block += f"   {param:<20} {value}\n"

        # Add pulse sequence parameters as variables
        pulse_sequence_vars = set()
        if self.pulse_sequence:
            for param, value in sorted(self.pulse_sequence.parameters.items()):
                if param.startswith('variable_'):
                    var_name = param.replace('variable_', '')
                    pulse_sequence_vars.add(var_name)
                    # Only add as variable if it's not already a standard parameter
                    if var_name not in required_params:
                        par_block += f"   variable {var_name:<15} {value}\n"

        # Add other variables from parameters
        for param, value in sorted(self.parameters.items()):
            if param.startswith('variable_'):
                var_name = param.replace('variable_', '')
                if var_name not in pulse_sequence_vars:
                    par_block += f"   variable {var_name:<15} {value}\n"

        # Add any remaining parameters
        processed_params = required_params | {k for k in self.parameters if k.startswith('variable_')} | pulse_sequence_vars
        remaining_params = {k: v for k, v in self.parameters.items()
                          if k not in processed_params}

        for param, value in sorted(remaining_params.items()):
            if not param.startswith('out_'):
                par_block += f"   {param:<20} {value}\n"

        par_block += "}\n"
        return par_block

    def generate_pulseq(self):
        """
        Generates the pulseq section of the SIMPSON input file.

        Returns:
            str: The pulseq section as a string.
        """
        if not self.pulse_sequence:
            # Return empty pulseq block if no sequence provided
            return "proc pulseq {} {}\n"

        return self.pulse_sequence.generate_code()

    def generate_main(self):
        """
        Generates the main section of the SIMPSON input file.
    
        Returns:
            str: The main section as a string.
        """

        out_format = self.parameters.get('out_format',
                     self.output_config.get('format', 'spe'))


        out_name = self.parameters.get('out_name',
                     self.output_config.get('name', '$par(name)'))

        lb = self.parameters.get('lb',
             self.output_config.get('lb', 0))

        zerofill = self.parameters.get('zerofill',
                  self.output_config.get('zerofill', 0))


        indent = "    "

        if out_format == "fid":
            return f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}faddlb $f {lb} 0
{indent}fzerofill $f {zerofill}
{indent}fsave $f {out_name}.fid
}}
"""
        elif out_format == "spe":
            main_block = f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}faddlb $f {lb} 0
{indent}fzerofill $f {zerofill}
{indent}fft $f
"""
            if 'variable_ref' in self.parameters or 'ref' in self.parameters:
                main_block += f"{indent}fset $f -ref $par(ref)\n"

            main_block += f"{indent}fsave $f {out_name}.spe\n}}\n"
            return main_block
        elif out_format == "xreim":
            return f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}fsave $f {out_name}.xreim -xreim
}}
"""
        else:
            raise ValueError(f"Unknown out_format '{out_format}'. Supported formats: 'fid', 'spe', 'xreim'")

    def save(self, filepath):
        """Save the SIMPSON input file"""
        with open(filepath, 'w') as file:
            file.write(str(self))

    def print(self):
        """Print the SIMPSON input file to console"""
        print(str(self))

    def run(self, filepath=None, timeout=None, read_output=False, delete_files=False, b0=None, nucleus=None, simpson_path=None, dry_run=False):
        """
        Run the SIMPSON simulation and optionally read the results.
        
        Args:
            filepath (str, optional): Path to save the input file. If None, a temporary file will be created.
            timeout (int, optional): Timeout in seconds for the SIMPSON process.
            read_output (bool): Whether to read the output file after running. Default is False.
            delete_files (bool): Whether to delete the input and output files after reading. Default is False.
            b0 (str, optional): Magnetic field strength (e.g., '400MHz', '9.4T'). Only used if read_output=True.
                If None, will be automatically derived from the 'proton_frequency' parameter.
            nucleus (str, optional): Nucleus type (e.g., '1H', '13C'). Only used if read_output=True.
                If None, will be automatically extracted from the spin system definition.
            simpson_path (str, optional): Custom path to the SIMPSON executable. If provided, this path will be used
                instead of searching for SIMPSON in the PATH or common installation locations.
            dry_run (bool): If True, only generates the input file and prints the command without running SIMPSON.
            
        Returns:
            str or Simpy or None: 
                - If read_output=True: A Simpy object containing the simulation results
                - If read_output=False: The command-line output from SIMPSON as a string
                - If dry_run=True: The command that would be executed
            
        Raises:
            FileNotFoundError: If SIMPSON is not found in the PATH, common locations, or the provided path.
            subprocess.TimeoutExpired: If the simulation doesn't complete within the timeout.
            subprocess.CalledProcessError: If SIMPSON returns a non-zero exit code.
        """
        import os
        import re
        import shutil
        import subprocess
        import tempfile

        # Find the SIMPSON executable
        simpson_executable = None

        if simpson_path:
            if os.path.exists(simpson_path):
                simpson_executable = simpson_path
            else:
                raise FileNotFoundError(f"SIMPSON executable not found at specified path: {simpson_path}")

        if simpson_executable is None:
            simpson_executable = shutil.which("simpson")

        if simpson_executable is None and not dry_run:
            raise FileNotFoundError(
                "SIMPSON executable not found in PATH or at the specified path.\n"
                "Please ensure SIMPSON is installed and available in your PATH, or provide the correct path."
            )

        temp_file = None
        if filepath is None:
            temp_fd, temp_file = tempfile.mkstemp(suffix='.in')
            os.close(temp_fd)
            filepath = temp_file

        base_filepath = os.path.splitext(filepath)[0]

        self.save(filepath)

        if dry_run:
            cmd = [simpson_executable if simpson_executable else "simpson", filepath]
            print(f"Dry run: Generated input file at {filepath}")
            print(f"Command: {' '.join(cmd)}")
            return ' '.join(cmd)

        # Determine expected output filename/locations
        out_format = self.parameters.get('out_format',
                        self.output_config.get('format', 'spe'))

        out_name = self.parameters.get('out_name',
                     self.output_config.get('name', base_filepath))
        if out_name == '$par(name)':
            out_name = base_filepath

        output_filename = f"{os.path.basename(out_name)}.{out_format}"

        possible_locations = [
            os.path.join(os.path.dirname(filepath), output_filename),
            os.path.join(os.getcwd(), output_filename),
            f"{out_name}.{out_format}"
        ]

        try:
            cmd = [simpson_executable, filepath]

            result = subprocess.run(cmd,
                                 check=True,
                                 capture_output=True,
                                 text=True,
                                 timeout=timeout)

            # Read ouput from run
            if read_output:
                from simpyson.io import read_simp

                output_file = None
                for location in possible_locations:
                    if os.path.exists(location):
                        output_file = location
                        break

                if output_file is None:
                    raise FileNotFoundError(
                        f"SIMPSON output file not found. Looked in: {possible_locations}\n"
                        f"SIMPSON stdout: {result.stdout}\n"
                        f"SIMPSON stderr: {result.stderr}"
                    )

                # Auto-extract magnetic field if not provided
                if b0 is None and 'proton_frequency' in self.parameters:
                    proton_freq = self.parameters['proton_frequency']
                    if isinstance(proton_freq, (int, float)):
                        # Convert to MHz if in Hz
                        if proton_freq > 1e6:
                            b0 = f"{proton_freq/1e6:.1f}MHz"
                        else:
                            b0 = f"{proton_freq}MHz"
                    elif isinstance(proton_freq, str):
                        match = re.match(r'(\d+(?:\.\d+)?)\s*([kMGT]?Hz|[kMGT]?hz)?', proton_freq)
                        if match:
                            value, unit = match.groups()
                            value = float(value)
                            if not unit:
                                b0 = f"{value}MHz"
                            elif unit.lower() in ['hz', 'khz', 'mhz', 'ghz', 'thz']:
                                b0 = proton_freq
                            else:
                                b0 = f"{value}MHz"

                # Auto-extract nucleus information if not provided
                if nucleus is None:
                    spinsys_str = self.generate_spinsys()
                    # First nucleus on the channels line is the observed nucleus
                    channels_match = re.search(r'channels\s+([\w\s]+)', spinsys_str)
                    if channels_match:
                        # Split by whitespace and take the first nucleus
                        nuclei_list = channels_match.group(1).split()
                        if nuclei_list:
                            nucleus = nuclei_list[0]
                    if nucleus is None:
                        nuclei_match = re.search(r'nuclei\s+(\w+)', spinsys_str)
                        if nuclei_match:
                            nucleus = nuclei_match.group(1)

                # Read the output file
                sim_result = read_simp(output_file, format=out_format, b0=b0, nucleus=nucleus)

                return sim_result
            else:
                return result.stdout

        finally:
            if delete_files:
                if temp_file and os.path.exists(temp_file):
                    os.remove(temp_file)
                elif filepath and os.path.exists(filepath):
                    os.remove(filepath)

                for location in possible_locations:
                    if os.path.exists(location):
                        try:
                            os.remove(location)
                        except OSError:
                            pass

def simulate_spectrum(spinsys, delete_files=True, filepath=None, **kwargs):
    """
    Easily simulate a spectrum from a spin system with smart defaults.
    
    This function automatically calculates the spectral width and center frequency
    based on the chemical shifts in the spin system, and runs a 'no_pulse' simulation.
    
    Args:
        spinsys: Soprano SpinSystem object or string definition
        delete_files (bool): Whether to delete input/output files after simulation. Default is True.
        filepath (str, optional): Path to save the SIMPSON input file. If None, uses a temporary file.
        **kwargs: Additional parameters to override defaults.
                  Common parameters: proton_frequency, spin_rate, lb, zerofill
                  
    Returns:
        Simpy: The simulated spectrum object
    """
    import re

    from simpyson.converter import ppm2hz
    from simpyson.utils import get_spin

    # Defaults
    defaults = {
        'proton_frequency': 800e6,
        'spin_rate': 30e3,
        'start_operator': 'Inx',
        'detect_operator': 'Inp',
        'crystal_file': 'rep168',
        'gamma_angles': 8,
        'np': 4096,
        'method': 'direct',
        'verbose': 0,
        'out_format': 'spe',
        'lb': 100,
        'zerofill': 4096,
        'pulse_sequence': 'no_pulse'
    }

    params = defaults.copy()
    # Update with user kwargs
    user_detect_op = 'detect_operator' in kwargs
    params.update(kwargs)

    # Extract shifts to calculate sw and offset
    spinsys_str = ""
    if hasattr(spinsys, 'to_simpson'):
        spinsys_str = spinsys.to_simpson()
    else:
        spinsys_str = str(spinsys)

    # Parse shifts
    shifts = []
    for line in spinsys_str.splitlines():
        if line.strip().startswith('shift'):
            parts = line.split()
            if len(parts) > 2:
                # shift index iso ...
                val_str = parts[2]
                if val_str.endswith('p'):
                    shifts.append(float(val_str[:-1]))
                else:
                    try:
                        shifts.append(float(val_str))
                    except ValueError:
                        pass

    if shifts:
        min_shift = min(shifts)
        max_shift = max(shifts)
        center_ppm = (min_shift + max_shift) / 2
        width_ppm = (max_shift - min_shift)

        # Ensure minimum width
        if width_ppm < 10: width_ppm = 10

        # Convert to Hz
        b0 = None
        proton_freq = params['proton_frequency']
        if isinstance(proton_freq, (int, float)):
            if proton_freq > 1e6:
                b0 = f"{proton_freq/1e6:.1f}MHz"
            else:
                b0 = f"{proton_freq}MHz"

        nucleus = '1H' # Default
        channels_match = re.search(r'channels\s+([\w\s]+)', spinsys_str)
        if channels_match:
            nuclei_list = channels_match.group(1).split()
            if nuclei_list:
                nucleus = nuclei_list[0]

        # Check for quadrupolar nucleus (spin > 0.5)
        if not user_detect_op:
            try:
                spin = get_spin(nucleus)
                if spin > 0.5:
                    params['detect_operator'] = 'Inc'
            except Exception:
                pass

        # Calculate center frequency for offset and ref
        # SIMPSON's offset command ADDS to peak positions: raw = true + offset
        # To center peaks at 0 in raw coordinates, use offset = -center_hz
        # This way: raw = true + (-center_hz) = true - center_hz ≈ 0 for peaks near center
        center_hz = ppm2hz(center_ppm, b0, nucleus)
        offset_value = -center_hz  # Negate to center peaks at 0
        
        # Calculate the spectral width needed to cover all peaks
        # With offset centering peaks at 0, SW only needs to cover the peak width
        min_hz = ppm2hz(min_shift, b0, nucleus)
        max_hz = ppm2hz(max_shift, b0, nucleus)
        width_hz = abs(max_hz - min_hz)
        required_sw = width_hz * 2
        
        # Round SW up to nearest multiple of spinning rate (n * spin_rate)
        spin_rate = params['spin_rate']
        n = 1
        while n * spin_rate < required_sw:
            n += 1
        sw_hz = n * spin_rate

        # Update params if not provided by user
        # Use variable_offset and variable_ref so they become variables in the par block
        # offset_value centers peaks at 0 in raw coordinates
        # ref should equal offset_value so that: hz = raw - ref restores true Hz
        if 'sw' not in kwargs:
            params['sw'] = sw_hz
        if 'variable_offset' not in kwargs:
            params['variable_offset'] = offset_value
        if 'variable_ref' not in kwargs:
            params['variable_ref'] = offset_value

    # Create calculator and run
    calc = SimpCalc(spinsys, **params)
    return calc.run(read_output=True, filepath=filepath, delete_files=delete_files)
