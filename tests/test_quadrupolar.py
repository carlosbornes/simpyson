import unittest
from simpyson.calculator import simulate_spectrum, SimpCalc
from simpyson.utils import get_spin

class TestQuadrupolarSupport(unittest.TestCase):
    def test_get_spin(self):
        """Test get_spin function"""
        self.assertEqual(get_spin('1H'), 0.5)
        self.assertEqual(get_spin('13C'), 0.5)
        self.assertEqual(get_spin('23Na'), 1.5) # 3/2
        self.assertEqual(get_spin('27Al'), 2.5) # 5/2
        self.assertEqual(get_spin('14N'), 1.0)
        
    def test_simulate_spectrum_quadrupolar(self):
        """Test that simulate_spectrum sets Inc for quadrupolar nuclei"""
        
        # 23Na 3/2
        spinsys = """
        channels 23Na
        nuclei 23Na
        shift 1 0 0 0 0 0 0
        """
        
        original_SimpCalc = SimpCalc
        
        captured_params = {}
        
        class MockSimpCalc:
            def __init__(self, spinsys, **kwargs):
                captured_params.update(kwargs)
                self.spinsys = spinsys
                self.parameters = kwargs
                
            def run(self, **kwargs):
                return "Mock Result"
                
        import simpyson.calculator
        simpyson.calculator.SimpCalc = MockSimpCalc
        
        try:
            simulate_spectrum(spinsys)
            self.assertEqual(captured_params.get('detect_operator'), 'Inc')
            
            # Test override
            captured_params.clear()
            simulate_spectrum(spinsys, detect_operator='Inp')
            self.assertEqual(captured_params.get('detect_operator'), 'Inp')
            
            # Test spin 1/2 (1H)
            captured_params.clear()
            spinsys_1H = """
            channels 1H
            nuclei 1H
            shift 1 0 0 0 0 0 0
            """
            simulate_spectrum(spinsys_1H)
            self.assertEqual(captured_params.get('detect_operator'), 'Inp')
            
        finally:
            # Restore
            simpyson.calculator.SimpCalc = original_SimpCalc

if __name__ == '__main__':
    unittest.main()
