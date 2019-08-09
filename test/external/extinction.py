import unittest
from astropy import units as u
from Phosphorpy.external.extinction import get_extinctions


class TestExtinction(unittest.TestCase):

    def setUp(self) -> None:
        self.ra = [12, 14, 15.3]
        self.dec = [1, 12.0, 23.0]

    def test_get_extinctions(self):
        get_extinctions(self.ra[0], self.dec[0])

        get_extinctions(self.ra, self.dec)

        get_extinctions(self.ra[0], self.dec[0], 0.4717, 'a')

        with self.assertRaises(ValueError):

            get_extinctions(self.ra[0], self.dec[0], 4717*u.AA)


