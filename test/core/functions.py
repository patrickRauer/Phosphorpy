import unittest
import numpy as np

from Phosphorpy.core import functions


class TestFunctions(unittest.TestCase):

    def test_gauss(self):
        self.assertEqual(functions.gaus(0, 0, 0, 1, 0, 0), 0)
        self.assertEqual(functions.gaus(0, 1, 0, 1, 0, 0), 1)
        self.assertEqual(functions.gaus(0, 0, 0, 1, 0, 1), 1)
        self.assertEqual(functions.gaus(0, 1, 0, 1, 0, 1), 2)
        self.assertEqual(functions.gaus(1, 1, 1, 1, 0, 0), 1)

        self.assertTrue(np.isnan(functions.gaus(0, 1, 0, 0, 0, 0)))

    def test_power_2_10(self):
        self.assertEqual(functions.power_2_10(0), 1)
        self.assertEqual(functions.power_2_10(2.5), 0.1)
        self.assertEqual(functions.power_2_10(-2.5), 10)

    def test_subtract(self):
        self.assertEqual(functions.subtract(1, 1), 0)

        np.testing.assert_array_equal(functions.subtract(np.arange(10),
                                                         np.arange(10)),
                                      np.zeros(10))
