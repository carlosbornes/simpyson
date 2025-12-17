import unittest
from simpyson.calculator import SimpCalc

class TestSimpCalcImprovements(unittest.TestCase):
    def test_validation_spinsys(self):
        """Test that spinsys cannot be None"""
        with self.assertRaises(ValueError):
            SimpCalc(spinsys=None)

    def test_validation_pulse_sequence(self):
        """Test invalid pulse sequence types"""
        with self.assertRaises(ValueError):
            SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence=123)

    def test_missing_parameters(self):
        """Test missing required parameters"""
        calc = SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence="pulse_90")
        with self.assertRaises(ValueError) as cm:
            calc.generate_par()
        self.assertIn("Missing required parameters", str(cm.exception))

    def test_dry_run(self):
        """Test dry_run option"""
        calc = SimpCalc(
            spinsys="spinsys { channels 1H }",
            pulse_sequence="pulse_90",
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
        # Should not raise FileNotFoundError even if simpson is missing
        cmd = calc.run(dry_run=True)
        self.assertTrue(isinstance(cmd, str))
        self.assertIn("simpson", cmd)

if __name__ == '__main__':
    unittest.main()
