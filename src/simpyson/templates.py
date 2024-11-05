class SimpSim:
    """
    Class to create a SIMPSON simulation input.

    Attributes:
        spinsys (str): Spin system from Soprano or a custom string.
        out_name (str): Output file name.
        out_format (str): Output format (fid, spe, xreim).
        spin_rate (float): Spin rate in Hz.
        np (int): Number of points.
        proton_freq (float): Proton frequency in Hz.
        start_op (str): Start operator.
        detect_op (str): Detect operator.
        crystal_file (str): Crystal file.
        gamma_angles (str): Gamma angles.
        sw (float): Spectral width.
        verbose (bool): Verbose output.
        lb (float): Line broadening.
        zerofill (int): Zero filling.
        method (str): Method of simulation (direct, indirect, ...).
        tsw (str, optional): Spectral width in time domain.
        pulse_sequence (str, optional): Pulse sequence from templates or a custom string.
        pH (float, optional): pulse for H in us.
        pX (float, optional): pulse for X in us.
        pY (float, optional): pulse for Y in us.
        plH (float, optional): power level for H in Hz.
        plX (float, optional): power level for X in Hz.
        plY (float, optional): power level for Y in Hz.
        phH (float, optional): phase for H pH.
        phX (float, optional): phase for X pH.
        phY (float, optional): phase for Y pH.

    Example:
        sim = SimpSim(spinsys=spinsys, out_name='output', out_format='spe', spin_rate=15e3, np=2048, proton_freq=400e6, start_op='Inx', detect_op='Inp', crystal_file='rep100', gamma_angles=4, sw=20e3, verbose=0, lb=20, zerofill=4096)
    """
    def __init__(self, spinsys, out_name, out_format, spin_rate, np, proton_freq, start_op, detect_op, crystal_file, gamma_angles, sw, verbose, lb, zerofill, method="direct", tsw=None, pulse_sequence=None, pH=None, pX=None, pY=None, plH=None, plX=None, plY=None, phH=None, phX=None, phY=None):
        self.spin_rate = spin_rate
        self.spinsys = spinsys
        self.out_name = out_name
        self.out_format = out_format
        self.np = np
        self.proton_freq = proton_freq
        self.start_op = start_op
        self.detect_op = detect_op
        self.crystal_file = crystal_file
        self.gamma_angles = gamma_angles
        self.sw = sw
        self.verbose = verbose
        self.tsw = tsw if tsw else f"1e6/{sw}"
        self.lb = lb
        self.zerofill = zerofill
        self.method = method
        self.pulse_sequence = pulse_sequence
        self.pH = pH
        self.pX = pX
        self.pY = pY
        self.plH = plH
        self.plX = plX
        self.plY = plY
        self.phH = phH
        self.phX = phX
        self.phY = phY

    def par_content(self):
        par_block = f"""
par {{
    spin_rate        {self.spin_rate}
    np               {self.np}
    proton_frequency {self.proton_freq}
    start_operator   {self.start_op}
    detect_operator  {self.detect_op}
    method           {self.method}
    crystal_file     {self.crystal_file}
    gamma_angles     {self.gamma_angles}
    variable sw      {self.sw}
    verbose          {self.verbose}
    variable tsw     {self.tsw}
"""
        if self.pH is not None:
            if self.plH is None or self.phH is None:
                raise ValueError("plH and phH must be defined if pH is defined")
            par_block += f"    pH               {self.pH}\n"
            par_block += f"    plH              {self.plH}\n"
            par_block += f"    phH              {self.phH}\n"

        if self.pX is not None:
            if self.plX is None or self.phX is None:
                raise ValueError("plX and phX must be defined if pX is defined")
            par_block += f"    pX               {self.pX}\n"
            par_block += f"    plX              {self.plX}\n"
            par_block += f"    phX              {self.phX}\n"

        if self.pY is not None:
            if self.plY is None or self.phY is None:
                raise ValueError("plY and phY must be defined if pY is defined")
            par_block += f"    pY               {self.pY}\n"
            par_block += f"    plY              {self.plY}\n"
            par_block += f"    phY              {self.phY}\n"

        par_block += "}\n"
        return par_block

    def pulseq_content(self):
        if self.pulse_sequence:
            return self.pulse_sequence
        else:
            return no_pulse

    def main_content(self):
        if self.out_format == "fid":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    faddlb $f {self.lb} 0
    fzerofill $f {self.zerofill}
    fsave $f {self.out_name}.fid
}}
"""
        elif self.out_format == "spe":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    faddlb $f {self.lb} 0
    fzerofill $f {self.zerofill}
    fft $f
    fsave $f {self.out_name}.spe
}}
"""
        elif self.out_format == "xreim":
            return f"""
proc main {{}} {{
    global par
    set f [fsimpson]
    fsave $f {self.out_name}.xreim -xreim
}}
"""
        else:
            raise ValueError(f"Unknown out_format {self.out_format}")

    def __str__(self):
        return f"{self.spinsys}\n{self.par_content()}{self.pulseq_content()}{self.main_content()}"
    
    def save(self, filepath):
        with open(filepath, 'w') as file:
            file.write(str(self))


# Predefined pulse sequences
# No pulse sequence
no_pulse = """
proc pulseq {} {
    global par
    acq_block {
    delay $par(tsw)
    }
}
"""

# 90 degree pulse
pulse_90 = """
proc pulseq {} {
    global par
    acq_block {
    pulse $par(pH) $par(plH) $par(phH)  
    delay $par(tsw)
    }
}
"""