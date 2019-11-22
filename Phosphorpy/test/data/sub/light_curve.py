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

    def test_align(self):
        lc = self.lc.align_light_curves(inplace=False)
        stats = lc.stats()
        stats0 = self.lc.stats()
        for i in stats.index.get_level_values(0).values:
            st = stats.loc[i]
            self.assertEqual(
                np.round(np.std(st[('mag', 'median')].values), 6), 0
            )
            st0 = stats0.loc[i]
            self.assertIsNot(
                np.round(np.std(st0[('mag', 'median')].values), 6), 0
            )

        self.lc.align_light_curves()
        stats = self.lc.stats()
        for i in stats.index.get_level_values(0).values:
            st = stats.loc[i]
            self.assertEqual(
                np.round(np.std(st[('mag', 'median')].values), 6), 0
            )

    def test_get_light_curve(self):
        l = self.lc.get_light_curve(0)
        l.plot.light_curve(0)

        self.lc.plot.light_curve([0, 1])

    def test_plot_light_curve(self):
        self.lc.plot.light_curve(0)

        self.lc.plot.light_curve(0, min_mjd=58000, max_mjd=58500)

        self.lc.plot.light_curve([0, 1])

        self.lc.plot.light_curve([0, 1], min_mjd=58000, max_mjd=58500)


if __name__ == '__main__':
    unittest.main()
