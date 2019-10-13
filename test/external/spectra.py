from astropy.coordinates import SkyCoord
from astropy import units as u
import unittest
import numpy as np

from Phosphorpy.external import spectra


class TestGetSDSSSpectra(unittest.TestCase):

    def setUp(self) -> None:
        self.coord = SkyCoord(114.084986*u.deg, 25.144718*u.deg)
        self.coords = SkyCoord(
            np.array([114.084986, 247.083831])*u.deg,
            np.array([25.144718, 40.933285])*u.deg
        )

    def test_single_spectra(self):
        print('single coordinate')
        sp = spectra.get_sdss_spectra(self.coord)
        print(sp)

    def test_single_coords_with_ids(self):
        sp = spectra.get_sdss_spectra(self.coord, np.arange(1))

    def test_multiple_coords(self):
        print('multiple coordinate')
        sp = spectra.get_sdss_spectra(self.coords)
        print(sp)

    def test_multiple_coords_with_ids(self):
        sp = spectra.get_sdss_spectra(self.coords, np.arange(len(self.coords)))

        with self.assertRaises(ValueError):
            sp = spectra.get_sdss_spectra(self.coords, np.arange(len(self.coords)-1))

        with self.assertRaises(ValueError):
            sp = spectra.get_sdss_spectra(self.coords, np.arange(len(self.coords)+1))

        sp = spectra.get_sdss_spectra(self.coords, list(np.arange(len(self.coords))))


class TestGetGAMASpectra(unittest.TestCase):

    def setUp(self) -> None:
        self.coord = SkyCoord(114.084986*u.deg, 25.144718*u.deg)
        self.coord2 = SkyCoord(211.73487*u.deg, -1.59471*u.deg)
        self.coords = SkyCoord(
            np.array([114.084986, 247.083831, 211.73487])*u.deg,
            np.array([25.144718, 40.933285, -1.59471])*u.deg
        )

    def test_single_spectra(self):
        print('single coordinate')
        sp = spectra.get_gama_spectra(self.coord2)
        sp = spectra.get_gama_spectra(self.coord)
        print(sp)

    def test_multiple_coords(self):
        print('multiple coordinate')
        sp = spectra.get_gama_spectra(self.coords)
        print(sp)
