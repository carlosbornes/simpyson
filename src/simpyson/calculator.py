from simpyson.templates import get_template, pulseq_templates, CustomPulseSequence, PulseSequenceTemplate


# Simpson calculator
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
                
                # Create template with extracted parameters
                return get_template(pulse_sequence, **pulseq_params)
            else:
                # Custom string sequence
                return CustomPulseSequence(pulse_sequence)
                
        elif isinstance(pulse_sequence, PulseSequenceTemplate):
            return pulse_sequence
            
        else:
            raise ValueError("pulse_sequence must be a template name, a custom string, or " \
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
            if self.spinsys.startswith("spinsys"):
                self.spinsys = self.spinsys
            elif self.spinsys.startswith("nuclei"):
                spinsys_block = "spinsys {"
                spinsys_block += f'\n{self.spinsys}\n'
                spinsys_block += "}\n"
                self.spinsys = spinsys_block
            else:
                raise ValueError("Invalid spinsys format. It should start with 'spinsys' or 'nuclei'.")
        else:
            raise ValueError("spinsys must be a string or a Soprano SpinSystem object.")
    
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
        if self.pulse_sequence:
            for param, value in sorted(self.pulse_sequence.parameters.items()):
                if param.startswith('variable_'):
                    var_name = param.replace('variable_', '')
                    par_block += f"   variable {var_name:<15} {value}\n"
        
        # Add other variables from parameters
        for param, value in sorted(self.parameters.items()):
            if param.startswith('variable_'):
                var_name = param.replace('variable_', '')
                par_block += f"   variable {var_name:<15} {value}\n"
        
        # Add any remaining parameters
        processed_params = required_params | {k for k in self.parameters if k.startswith('variable_')}
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
            Warning("No pulse sequence provided was provided.")

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
            return f"""
proc main {{}} {{
{indent}global par
{indent}set f [fsimpson]
{indent}faddlb $f {lb} 0
{indent}fzerofill $f {zerofill}
{indent}fft $f
{indent}fsave $f {out_name}.spe
}}
"""
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
