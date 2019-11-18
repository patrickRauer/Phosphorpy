import pandas as pd
import unittest

from Phosphorpy.data.sub.tables.magnitude import Magnitude


class TestMagnitude(unittest.TestCase):

    def setUp(self) -> None:
        data = pd.DataFrame({'umag': [1, 2, 4., -9.],
                             'gmag': [-1.2, -4.3, -6.9, 19.9]})
        names = ['umag', 'gmag']
        self.magnitude = Magnitude(data, names, 'sdss')

    def test_str(self):
        self.assertEqual(str(self.magnitude),
                         'Magnitude of sdss with 4 entries\n')
