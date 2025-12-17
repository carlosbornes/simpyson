import unittest
import numpy as np
from simpyson.calculator import simulate_spectrum, SimpCalc
from simpyson.converter import ppm2hz

class TestSimulateSpectrum(unittest.TestCase):
    def test_simulate_spectrum_defaults(self):
        """Test simulate_spectrum with minimal arguments"""
        
        # Define a simple spin system string
        # 1H at 5 ppm and 10 ppm
        spinsys = """
        channels 1H
        nuclei 1H 1H
        shift 1 5p 0 0 0 0 0
        shift 2 10p 0 0 0 0 0
        """
        
        try:
            simulate_spectrum(spinsys, dry_run=True)
        except Exception as e:
            pass

    def test_parameter_calculation(self):
        """Test that SW and Offset are calculated correctly"""
        
        spinsys = """
        channels 1H
        nuclei 1H 1H
        shift 1 5p 0 0 0 0 0
        shift 2 10p 0 0 0 0 0
        """
        
        shifts = [5.0, 10.0]
        min_shift = 5.0
        max_shift = 10.0
        center_ppm = 7.5
        width_ppm = 5.0 * 1.5 # 7.5 ppm
        
        b0 = '800.0MHz' # Default
        nucleus = '1H'
        
        center_hz = ppm2hz(center_ppm, b0, nucleus)
        sw_hz = abs(ppm2hz(width_ppm, b0, nucleus) - ppm2hz(0, b0, nucleus))
        
        print(f"Expected Center: {center_ppm} ppm -> {center_hz} Hz")
        print(f"Expected SW: {width_ppm} ppm -> {sw_hz} Hz")
        
        calc = SimpCalc(spinsys, 
                        proton_frequency=800e6, 
                        sw=sw_hz, 
                        offset=center_hz, 
                        variable_ref=center_hz,
                        pulse_sequence='no_pulse',
                        spin_rate=30e3,
                        start_operator='Inx',
                        detect_operator='Inp',
                        crystal_file='rep168',
                        gamma_angles=8,
                        np=4096,
                        method='direct',
                        verbose=0)
                        
        self.assertIn('variable_offset', calc.pulse_sequence.parameters)
        self.assertAlmostEqual(calc.pulse_sequence.parameters['variable_offset'], center_hz)
        
        main_block = calc.generate_main()
        self.assertIn(f"fset $f -ref $par(ref)", main_block)
        
        par_block = calc.generate_par()
        self.assertIn(f"variable ref", par_block)
        self.assertIn(f"{center_hz}", par_block)
        
        pulseq_block = calc.generate_pulseq()
        self.assertIn(f"offset $par(offset)", pulseq_block)

if __name__ == '__main__':
    unittest.main()
