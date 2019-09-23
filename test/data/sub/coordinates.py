from astropy.table import Table
import numpy as np
import pandas as pd

import unittest

from Phosphorpy.data.sub.coordinates import CoordinateTable


class TestCoordinateTable(unittest.TestCase):

    def setUp(self) -> None:
        data = pd.DataFrame({
            'ra': 360*np.random.random(10),
            'dec': 180*np.random.random(10)-90
        })
        self.data = data
        self.coord = CoordinateTable(data.copy())

    def test_to_table(self):
        self.assertTrue(
            (self.coord.to_table() == self.data).all().all()
        )
        self.assertEqual(
            len(self.coord.to_table(True).columns), 4
        )
        self.assertEqual(
            list(self.coord.to_table(True).columns.values),
            ['ra', 'dec', 'l', 'b']
        )

    def test_to_sky_coord(self):
        self.assertEqual(
            len(self.coord.to_sky_coord()), len(self.data)
        )

        index = 1
        c = self.coord.to_sky_coord(index)
        self.assertEqual(
            c.ra.degree, self.data['ra'].values[index]
        )
        self.assertEqual(
            c.dec.degree, self.data['dec'].values[index]
        )

    def test_to_astropy_table(self):
        self.assertTrue(
            (self.coord.to_astropy_table()[['ra', 'dec']] ==
             Table.from_pandas(self.data)).all()
        )

        self.assertEqual(
            self.coord.to_astropy_table().meta['category'],
            'coordinates'
        )

    def test_write(self):
        raise NotImplementedError()
