import unittest
from astropy.table import Table, vstack
import numpy as np
import warnings
import random


from Phosphorpy.config import match

warnings.simplefilter('ignore')


class TestMatch(unittest.TestCase):

    def setUp(self) -> None:
        size = 10
        self.data = Table({'ra': np.arange(size)+1.0,
                           'dec': np.arange(size)+1.0,
                           'label': np.arange(size)+1.0})

    def test_to_pandas(self):
        self.assertTrue((match.to_pandas(self.data) ==
                         self.data.to_pandas()).all().all())
        self.assertTrue((match.to_pandas(self.data.to_pandas()) ==
                         self.data.to_pandas()).all().all())

        with self.assertRaises(TypeError):
            match.to_pandas('a')

    def test_fill_missing_labels(self):
        match.fill_missing_labels(self.data, 10)

        with self.assertRaises(ValueError):
            match.fill_missing_labels(self.data[['ra', 'dec']], 10)

        with self.assertRaises(ValueError):
            match.fill_missing_labels(self.data[['ra', 'dec']].to_pandas(), 10)

        with self.assertRaises(TypeError):
            match.fill_missing_labels('a', 10)

    def test_convert_input_data(self):
        data_pd = self.data.to_pandas()
        self.assertTrue((match.convert_input_data(data_pd) ==
                         self.data).all())
        self.assertTrue((match.convert_input_data(self.data) ==
                         self.data).all())

    def test_next_neighbour_id(self):
        rand_sample = random.sample(range(0, len(self.data)),
                                    len(self.data))
        d = self.data[[rand_sample]]
        ids = match.next_neighbour_id(self.data, d,
                                      'ra', 'dec', 'ra', 'dec')
        self.assertTrue((d == self.data[ids]).all())

    def test_match_catalogs(self):
        matched = match.match_catalogs(self.data, self.data,
                                       'ra', 'dec', 'ra', 'dec')
        self.assertEqual(len(matched), len(self.data))
        self.assertTrue((matched['ra'] == self.data['ra']).all())

    def test_group_by_coordinates(self):
        d = vstack([self.data, self.data])
        grouped = match.group_by_coordinates(d, 'ra', 'dec')

        self.assertEqual(len(grouped), len(self.data))
        print(grouped == self.data)
        self.assertTrue((grouped == self.data).all())
