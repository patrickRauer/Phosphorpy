from astropy.coordinates import SkyCoord
from astropy import units as u
<<<<<<< HEAD
import unittest
import os
import shutil
import numpy as np

from Phosphorpy.external import image


class TestImage(unittest.TestCase):

    def setUp(self) -> None:
        self.image = image.Image()

    def test_check_size(self):

        self.assertEqual(self.image._check_size(1*u.arcmin), 1*u.arcmin)
        self.assertEqual(self.image._check_size(1*u.arcmin), 1/60*u.deg)
        self.assertEqual(self.image._check_size(1*u.arcmin), 60*u.arcsec)
        # size in wrong units
        with self.assertRaises(ValueError):
            self.image._check_size(2*u.m)
        # too small
        with self.assertRaises(ValueError):
            self.image._check_size(2*u.arcsec)
        # to big
        with self.assertRaises(ValueError):
            self.image._check_size(2*u.deg)

    def test_check_bands(self):
        self.image._available_bands = ['g', 'r', 'i', 'z', 'y']
        bands = [['g', 'r', 'i'], ['i', 'z', 'y']]
        for b in bands:
            self.assertEqual(self.image._check_bands_(b), b)

        # wrong band
        with self.assertRaises(ValueError):
            self.image._check_bands_(['u', 'g', 'r'])

        # not enough bands
        with self.assertRaises(ValueError):
            self.image._check_bands_(['u', 'g'])

        # too many bands
        with self.assertRaises(ValueError):
            self.image._check_bands_(['u', 'g', 'r', 'i'])


class TestPanstarrsImage(unittest.TestCase):
    ps_image = None
    coordinates_m92 = None
    coordinates_outside_ps = None

    def setUp(self) -> None:
        self.coordinates_m92 = SkyCoord('17h17m07.39s +43d08m09.4s')
        self.coordinates_outside_ps = SkyCoord(100*u.deg, -50*u.deg)
        self.ps_image = image.PanstarrsImage()

    def test_get_normalized_imaged(self):
        # default case with coordinates of M92
        rgb_image, hdus = self.ps_image.get_normalized_imaged(self.coordinates_m92, 0)

        self.assertFalse(np.isnan(rgb_image).any())

        # outside of the covered filed of PS1
        with self.assertRaises(ValueError):
            self.ps_image.get_normalized_imaged(self.coordinates_outside_ps, 0)

        # smoothing
        # below the lower limit
        with self.assertRaises(ValueError):
            self.ps_image.get_normalized_imaged(self.coordinates_m92, -1)

        # above the upper limit
        with self.assertRaises(ValueError):
            self.ps_image.get_normalized_imaged(self.coordinates_m92, 20)

    def test_get_color_image(self):
        test_path = './test/test.fits'

        # save as a bunch of fits files
        self.ps_image.get_color_image(self.coordinates_m92, test_path)
        self.assertTrue(os.path.exists('./test/test_g.fits') &
                        os.path.exists('./test/test_r.fits') &
                        os.path.exists('./test/test_z.fits'))
        shutil.rmtree('./test/')

        # save as a colored png
        self.ps_image.get_color_image(self.coordinates_m92, './test/img.png')

        # just show
        self.ps_image.get_color_image(self.coordinates_m92)

    def tearDown(self) -> None:
        for p in ['./temp/', './test/']:
            try:
                shutil.rmtree(p)
            except TypeError:
                pass
            except FileNotFoundError:
                pass


from Phosphorpy.external.image import PanstarrsImage, SDSSImage
import unittest


class TestPSImage(unittest.TestCase):

    def setUp(self) -> None:
        self.ps = PanstarrsImage()
        self.coord = SkyCoord(163.05867*u.deg, -24.79324*u.deg)

    def test_get_normalized_image(self):
        self.ps.get_normalized_imaged(self.coord, 2)

        with self.assertRaises(ValueError):
            self.ps.get_normalized_imaged(self.coord, -2)

        with self.assertRaises(ValueError):
            self.ps.get_normalized_imaged(self.coord, 20)

        with self.assertRaises(ValueError):
            self.ps.get_normalized_imaged(self.coord, 2, temp_path=-1)

    def test_get_color_image(self):
        self.ps.get_color_image(self.coord)
        self.ps.get_color_image(self.coord, smooth=6)


def test_smoothing():
    ps = PanstarrsImage()
    coord = SkyCoord(163.05867 * u.deg, -24.79324 * u.deg)
    ps.get_color_image(coord)
    # ps.get_color_image(coord, smooth=6)
<<<<<<< HEAD
>>>>>>> testing
=======
>>>>>>> master
