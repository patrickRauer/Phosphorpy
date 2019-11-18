from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import unittest

from Phosphorpy.external import vizier
from Phosphorpy.config import names


class TestVizier(unittest.TestCase):

    def setUp(self) -> None:
        self.survey = vizier.Vizier(name=names.TWO_MASS)
        self.m92_coordinates = SkyCoord('17h17m7.390s +43d08m9.40s')
        self.m3_coordinates = SkyCoord('13h42m11.620s +28d22m38.2s')
        self.m92_m3_coordinates = SkyCoord(['17h17m7.390s +43d08m9.40s',
                                            '13h42m11.620s +28d22m38.2s'])

        self.table = Table()
        self.table['ra'] = self.m92_m3_coordinates.ra.degree
        self.table['dec'] = self.m92_m3_coordinates.dec.degree

    def test__adjust_coordinate_names__(self):
        survey = vizier.Vizier(name=names.TWO_MASS)
        ra_name = survey.ra_name
        dec_name = survey.dec_name

        survey.__adjust_coordinate_names__()

        self.assertEqual(survey.ra_name, f'{survey.name}_{ra_name}')
        self.assertEqual(survey.dec_name, f'{survey.name}_{dec_name}')

    def test___query__(self):
        rs = self.survey.__query__(self.m92_coordinates)

        self.assertEqual(len(rs), 1)

        rs = self.survey.__query__(self.m3_coordinates)

        self.assertEqual(len(rs), 1)

        rs = self.survey.__query__(self.m92_m3_coordinates)

    def test_query_coordinate(self):
        rs = self.survey.query_coordinate(self.m3_coordinates)

        rs = self.survey.query_coordinate([self.m3_coordinates.ra.degree,
                                           self.m3_coordinates.dec.degree])

        rs = self.survey.query_coordinate((self.m3_coordinates.ra.degree,
                                           self.m3_coordinates.dec.degree))

        with self.assertRaises(ValueError):
            self.survey.query_coordinate([1, 2, 3])

        with self.assertRaises(ValueError):
            self.survey.query_coordinate([1])

        with self.assertRaises(ValueError):
            self.survey.query_coordinate(self.m92_m3_coordinates)

    def test_query(self):
        with self.assertRaises(ValueError):
            self.survey.query('a')

        rs = self.survey.query(self.table)
        print(type(rs))
        rs = self.survey.query(self.table, blank=True)
        print(rs.columns)

        print(np.hypot(rs['ra_input']-self.m3_coordinates.ra.degree,
                       rs['dec_input']-self.m3_coordinates.dec.degree)*3600)

    def test_query_constrain(self):

        self.survey.query_constrain(Jmag='<1')
