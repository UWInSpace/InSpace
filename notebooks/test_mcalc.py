import unittest

from mcalc import get_prot_mass

class MolecularMassTestCase(unittest.TestCase):

    def test_lower(self):
        result_lower = get_prot_mass('mpfmvnniyvsfceikeivcaggsttkyadvlqenneqgrtvklq')
        self.assertEqual(results_lower, 5051.7509)
        
    def test_gaps(self):
        result_gaps = get_prot_mass('MPFMVNNIYVSF  CEIKEIV CAGGSTTKYADVLQEN NEQGRTVKLQ')
        self.assertEqual(result_gaps, 5051.7509)

    def test_numbers(self):
        # check that the function fails when fed values other than strings
        with self.assertRaises(AttributeError):
            get_prot_mass(123)

    def test_amino_acids(self):
        # check that the function fails when fed values other than strings containing real amino acid symbols
        with self.assertRaises(ValueError):
            get_prot_mass('MPFMVNNIYVSF528ZZZ')

if __name__ == '__main__':
    unittest.main()