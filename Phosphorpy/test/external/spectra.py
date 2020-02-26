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


class TestGetLAMOSTSpectra(unittest.TestCase):

    def setUp(self) -> None:
        # the sources have the following amount of spectra in LAMOST DR5
        # 1
        # 2
        # 0
        # 1
        self.coords = SkyCoord(
            np.array([332.368745, 149.0749413, 211.73487, 78.823442])*u.deg,
            np.array([-01.955771, +33.4598012, -1.59471, 38.263469])*u.deg
        )

    def test_single_spectra(self):
        # test the source with a single spectra
        sp = spectra.get_lamost_spectra(self.coords[0])
        self.assertEqual(len(sp), 1)

        # test the source with two spectra
        sp = spectra.get_lamost_spectra(self.coords[1])
        self.assertEqual(len(sp), 2)

        # test a source with one spectra
        sp = spectra.get_lamost_spectra(self.coords[2])
        self.assertEqual(len(sp), 0)

    def test_multiple_coords(self):
        # test all together
        sp = spectra.get_lamost_spectra(self.coords)
        self.assertEqual(len(sp), 4)
