import unittest
from simpyson.calculator import SimpCalc
from simpyson.templates import CustomPulseSequence

class TestCustomPulseSequence(unittest.TestCase):
    def test_custom_pulse_sequence_params(self):
        """Test that custom pulse sequence receives parameters from SimpCalc"""
        
        code = """
        pulse $par(my_param) 0 0 0 0
        """
        
        # Initialize SimpCalc with custom code and the parameter
        calc = SimpCalc(
            spinsys="spinsys { channels 1H }",
            pulse_sequence=code,
            proton_frequency=400e6,
            spin_rate=10000,
            start_operator="I1z",
            detect_operator="I1p",
            np=1024,
            sw=20000,
            method="direct",
            crystal_file="rep100",
            gamma_angles=10,
            verbose=0,
            my_param=10.0, # This must pass to CustomPulseSequence
            unused_param=20.0 # This should not
        )
        
        self.assertIn('variable_my_param', calc.pulse_sequence.parameters)
        self.assertEqual(calc.pulse_sequence.parameters['variable_my_param'], 10.0)
        
        self.assertNotIn('variable_unused_param', calc.pulse_sequence.parameters)
        
        # Verify generation
        par_block = calc.generate_par()
        import re
        self.assertRegex(par_block, r"variable\s+my_param\s+10\.0")
        
    def test_standard_params_in_custom_code(self):
        """Test using standard parameters in custom code"""
        code = """
        delay $par(np)
        """
        
        calc = SimpCalc(
            spinsys="spinsys { channels 1H }",
            pulse_sequence=code,
            proton_frequency=400e6,
            spin_rate=10000,
            start_operator="I1z",
            detect_operator="I1p",
            np=1024,
            sw=20000,
            method="direct",
            crystal_file="rep100",
            gamma_angles=10,
            verbose=0
        )
        
        # np is required by code
        self.assertIn('variable_np', calc.pulse_sequence.parameters)
        
        # It should appear in par block
        par_block = calc.generate_par()
        self.assertIn("np                   1024", par_block)
        
        # It should not appear as "variable np"
        self.assertNotIn("variable np", par_block)

if __name__ == '__main__':
    unittest.main()
