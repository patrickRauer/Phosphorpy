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
        gaia_test = pd.DataFrame({
            'phot_g_mean_mag': np.array([5., 6., 7.]),
            'phot_bp_mean_mag': np.array([5., 6., 7.]),
            'phot_rp_mean_mag': np.array([5., 6., 7.]),
            'parallax': np.array([1, 0.5, 0.25]),
            'parallax_error': np.array([1, 0.5, 0.25]),
            'phot_g_mean_flux_error': np.array([1e-17, 2e-17, 3e-17]),
            'phot_bp_mean_flux_error': np.array([1e-17, 2e-17, 3e-17]),
            'phot_rp_mean_flux_error': np.array([1e-17, 2e-17, 3e-17]),
            'phot_g_mean_flux': np.power(10., -0.4*np.array([5., 6., 7.])),
            'phot_bp_mean_flux': np.power(10., -0.4*np.array([5., 6., 7.])),
            'phot_rp_mean_flux': np.power(10., -0.4*np.array([5., 6., 7.])),
        })
        xmatch._compute_gaia_mags(gaia_test)

    def test___output_columns__(self):
        for k in xmatch.SURVEY_DATA.keys():
            xmatch.__output_columns__(k)

    def test_xmatch(self):
        coords = pd.DataFrame({
            'ra': np.array([114.084986, 247.083831, 044.9961159]),
            'dec': np.array([25.144718, 40.933285, 00.005620003])})
        check = {}
        for k in xmatch.SURVEY_DATA.keys():
            # exclude surveys/catalogs which are not available via XMatch
            if 'Bailer' not in k and 'SkyMapper' not in k and 'GPS1' not in k:
                print(k)
                try:
                    xmatch.xmatch(coords.copy(), 'ra', 'dec', k)
                    check[k] = True
                except ValueError:
                    check[k] = False
        print(check)
        if not all(check.items()):
            print(check)
            raise ValueError()


def test_large_file():
    try:
        path = '/Volumes/UNTITLED/sed/sdss_2mass/330_60_15.fits'
        data = Table.read(path)[['RAJ2000', 'DEJ2000']]
        rs = xmatch.xmatch(data.to_pandas(), 'RAJ2000', 'DEJ2000', 'GAIA')

        print(len(data), len(rs))
        pass
    except:
        pass


if __name__ == '__main__':
    unittest.main()
