import numpy as np
import pandas as pd

import unittest

from Phosphorpy.data.sub import astrometry
from Phosphorpy.data.sub.coordinates import CoordinateTable
from Phosphorpy.core.structure import Mask

from Phosphorpy import DataSet


def test_only_nearest():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078],
            'angDist': [0.1, 0.02, 1],
            'row_id': [1, 1, 2]
        }
    )
    d = astrometry._only_nearest(coord)
    assert len(d) == 2


def test_download_gaia_data():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    coord = CoordinateTable(coord, mask=Mask(len(coord)))

    bj = astrometry._download_gaia_data(coord)

    assert len(bj) <= len(coord)


def test_download_bailer_jones_data():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    coord = CoordinateTable(coord, mask=Mask(len(coord)))

    bj = astrometry._download_bailer_jones_data(coord.data)

    assert len(bj) <= len(coord)

    bj = astrometry._download_bailer_jones_data(coord)

    assert len(bj) <= len(coord)


def test_load_to_dataset():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    # coord = CoordinateTable(coord, mask=Mask(len(coord)))
    ds = DataSet(coordinates=coord)
    astrometry.AstrometryTable.load_to_dataset(ds)


def test_load_astrometry():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    coord = CoordinateTable(coord, mask=Mask(len(coord)))
    astro = astrometry.AstrometryTable.load_astrometry(coord)
    assert len(astro) <= len(coord)


class TestAstrometryTable(unittest.TestCase):

    def setUp(self) -> None:
        self.coord = pd.DataFrame(
            {
                'ra': [18.15629, 166.12397, 260.78809],
                'dec': [33.49162, 8.6418, 48.31078]
            }
        )
        self.coord = CoordinateTable(self.coord, mask=Mask(len(self.coord)))
        self.astro = astrometry.AstrometryTable.load_astrometry(self.coord)

    def test_proper_motion(self):
        pm = self.astro.proper_motion()
        self.assertLessEqual(len(pm), len(self.coord))
        self.assertEqual(
            list(pm.columns.values),
            ['pmra', 'pmdec', 'pmra_err', 'pmdec_err']
        )

        pm_cos = self.astro.proper_motion(True)
        self.assertTrue(
            (np.nan_to_num(np.abs(pm['pmra'].values)) >=
             np.nan_to_num(np.abs(pm_cos['pmra'].values))).all()
        )

    def test_total_proper_motion(self):
        pm = self.astro.total_proper_motion()
        self.assertLessEqual(
            len(pm), len(self.coord)
        )
        self.assertEqual(
            list(pm.columns.values),
            ['pm', 'err']
        )

    def test_distance(self):
        with self.assertRaises(ValueError):
            self.astro.distance('test')

        self.astro.distance()

        self.astro.distance('bj')

        self.astro.distance('simple')

        self.astro.distance('parallax')

    def test_set_distance_limit(self):
        self.astro.set_distance_limit(0.0, 0.5)

        with self.assertRaises(ValueError):
            self.astro.set_distance_limit(0.5, 0.0)

    def test_set_parallax_limit(self):
        self.astro.set_parallax_limit(0, 0.5)

        with self.assertRaises(ValueError):
            self.astro.set_parallax_limit(0.5, 0.)

        self.astro.set_parallax_limit(0, 0.5, True)

    def test_set_total_proper_motion_limit(self):
        self.astro.set_total_proper_motion_limit(0, 10)
        self.astro.set_total_proper_motion_limit(0, 10, True)

        with self.assertRaises(ValueError):
            self.astro.set_total_proper_motion_limit(10, 0)

    def test_set_proper_motion_limit(self):
        with self.assertRaises(ValueError):
            self.astro.set_proper_motion_limit('k', 0, 1)

        self.astro.set_proper_motion_limit('ra', -10, 10)

        with self.assertRaises(ValueError):
            self.astro.set_proper_motion_limit('ra', 10, -10)

        self.astro.set_proper_motion_limit('ra', -10, 10, True)

        self.astro.set_proper_motion_limit('dec', -10, 10)


if __name__ == '__main__':
    unittest.main()
