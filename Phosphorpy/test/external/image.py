from astropy.coordinates import SkyCoord
from astropy import units as u
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
