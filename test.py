import unittest
from week1 import *

class TestFunctions(unittest.TestCase):
    def setUp(self):
        self.text = 'TATGGGGTGC'

    def test_kmer_composition_func(self):
        result1 = kmer_composition(self.text, 3)
        expected1 = ['ATG', 'GGG', 'GGG', 'GGT', 'GTG', 'TAT', 'TGC', 'TGG']
        self.assertEqual(result1, expected1)

        result2 = kmer_composition('CAATCCAAC', 5)
        expected2 = ['CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC']
        self.assertEqual(result2, result2)

    def test_reconstruct_genome_from_path_func(self):
        genome_path = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
        result = reconstruct_genome_from_path(genome_path)
        expected = 'ACCGAAGCT'
        self.assertEqual(result, expected)

class TestDBGClass(unittest.TestCase):
    def setUp(self):
        self.g = DBG('AAGATTCTCTAAGA', 4)

    def test_dbg_initialization(self):
        self.assertIsInstance(self.g, DBG)


if __name__ == '__main__':
    unittest.main()
