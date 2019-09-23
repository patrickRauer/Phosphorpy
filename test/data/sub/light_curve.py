from astropy.table import Table
import numpy as np
import pandas as pd
import os
import unittest

from Phosphorpy.data.sub import light_curve


class TestLightCurve(unittest.TestCase):

    def setUp(self) -> None:
        self.coord = pd.DataFrame(
            {
                'ra': [18.15629, 166.12397, 260.78809],
                'dec': [33.49162, 8.6418, 48.31078]
            }
        )
        self.lc = light_curve.LightCurves(self.coord)

    def test_str(self):
        print(self.lc)

    def test_stats(self):
        print(self.lc.stats())

    def test_average(self):
        print(self.lc.average().stats())

        with self.assertRaises(ValueError):
            self.lc.average(-1)

    def test_get_light_curve(self):
        l = self.lc.get_light_curve(0)
        l.plot.light_curve(0)

        self.lc.plot.light_curve([0, 1])


class TestLightCurveDataSet(unittest.TestCase):

    def setUp(self) -> None:
        from Phosphorpy import DataSet
        Table.read('/Users/patrickr/Documents/temp/comb_ra_lte_180.fits')[:10].write('temp.fits', overwrite=True)
        self.ds = DataSet.load_coordinates('temp.fits', 'fits', 'ra', 'dec')
        # self.ds.load_from_vizier(['SDSS', '2MASS'])
        os.remove('temp.fits')

    def test_get_light_curves(self):
        print(self.ds.light_curves)
