from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import os
import unittest

from Phosphorpy.data.sub.coordinates import CoordinateTable


class TestCoordinateTable(unittest.TestCase):

    def setUp(self) -> None:
        data = pd.DataFrame({
            'ra': 360*np.random.random(10),
            'dec': 180*np.random.random(10)-90
        })
        self.data = data
        self.sk = SkyCoord(data['ra'].values*u.deg, data['dec'].values*u.deg)
        self.coord = CoordinateTable(data.copy())

    def test_initialization(self):
        coord = CoordinateTable(self.sk)
        coord = CoordinateTable(self.data.values)
        coord = CoordinateTable(Table.from_pandas(self.data))

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
        self.coord.write('./temp.fits', data_format='fits')
        self.coord.write('./temp.csv', data_format='csv')

        os.remove('./temp.fits')
        os.remove('./temp.csv')

    def test_get_item(self):
        for c in ['ra', 'dec']:
            self.assertTrue((self.coord[c] == self.data[c]).all())
        self.coord['l']
        self.coord['b']

        self.coord[1]

        with self.assertRaises(KeyError):
            self.coord['g']

    def test_len(self):
        self.assertEqual(
            len(self.coord), len(self.data)
        )

    def test_eq(self):
        with self.assertRaises(ValueError):
            self.coord == self.data

        self.assertTrue(
            self.coord == self.sk
        )
        self.assertFalse(
            self.coord == self.sk[1:]
        )

    def test_match(self):
        self.coord.match(self.data)

        self.coord.match(Table.from_pandas(self.data))

        self.coord.match(self.data.values)

    def test_combine_coordinates(self):
        self.coord._combine_coordinates(self.data)
        self.coord._combine_coordinates(self.data.values)
        self.coord._combine_coordinates(self.sk)

        with self.assertRaises(ValueError):
            self.coord._combine_coordinates(1)


if __name__ == '__main__':
    unittest.main()