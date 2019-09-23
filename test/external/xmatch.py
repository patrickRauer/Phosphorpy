import unittest
from astropy.table import Table, Column
import numpy as np
import pandas as pd
import os

from Phosphorpy.external import xmatch


class TestXMatch(unittest.TestCase):

    def setUp(self) -> None:
        self.data = Table({'ra': [12, 13., -4],
                           'dec': [0., 43., -6]})

    def test_check_row_id(self):
        xmatch.__check_row_id__(self.data, self.data.colnames)
        self.assertTrue('row_id' in self.data.colnames)
        np.testing.assert_array_equal(self.data['row_id'],
                                      Column(np.arange(len(self.data))+1))

    def test__write_temp_file__(self):
        tab = Table({'ra': [1, 2, 3],
                     'dec': [1, 2, 3]})
        xmatch.__write_temp_file__(tab.copy(), 'ra', 'dec')
        self.assertTrue(os.path.exists('temp.csv'))
        os.remove('temp.csv')
        xmatch.__write_temp_file__(tab.to_pandas(), 'ra', 'dec')
        self.assertTrue(os.path.exists('temp.csv'))
        os.remove('temp.csv')

        with self.assertRaises(TypeError):
            xmatch.__write_temp_file__('f', 'ra', 'dec')

    def test__output_columns__(self):

        pass

    def test_find_suffix(self):
        cols = (('ABmag', 'mag'),
                ('kmag', 'mag'),
                ('ABmag1', 'mag1'),
                ('kmag1', 'mag1'),
                ('ABAperMag', 'AperMag'),
                ('kAperMag', 'AperMag'),
                ('ABAperMag1', 'AperMag1'),
                ('kAperMag1', 'AperMag1'),
                ('ABrMag', 'Mag'),
                ('kMag', 'Mag'),
                ('ABMag1', 'Mag1'),
                ('kMag1', 'Mag1'))
        for c in cols:
            self.assertEqual(xmatch.find_suffix((c[0], )),
                             c[1])

    def test_compute_gaia_mags(self):
        pass

    def test_xmatch(self):
        pass
