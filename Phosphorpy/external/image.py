"""
This script provides an interface to download and create color images from SDSS and Pan-STARRS.
These images can be overplotted with markers of the required source and a proper motion arrow.
"""

import os

import numpy as np
import pylab as pl
from astropy import units as u
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS
from astroquery.sdss import SDSS
import shutil

from Phosphorpy.external.panstarrs import download_all_bands
from Phosphorpy.core.functions import smooth2d


def write(path, c_image, bands):
    """
    Writes the different images to a fits file

    :param path: Path to the images
    :param c_image: the images
    :param bands: the name of the bands
    :return:
    """
    if '.fit' in path:
        path = path.split('.fit')[0]

    if not os.path.exists(os.path.split(path)[0]):
        os.makedirs(os.path.split(path)[0])

    for img, b in zip(c_image, bands):
        img.writeto(f'{path}_{b}.fits', overwrite=True)


def max_count(rgb, center, box_size=2):
    """
    Finds the maximal count in the box for every layer
    :param rgb: The rgb layers
    :type rgb: ndarray
    :param center: The central coordinates (pixel)
    :type center: tuple, list
    :param box_size: The half box size
    :type box_size: bool
    :return: The maximal counts
    :rtype: list
    """
    center_counts = []
    for i in range(3):
        im = rgb[:, :, i].copy()

        med = 2 * np.nanmedian(np.nanmean(im, axis=0))
        im[im < med] = med
        # im -= med

        # collect the center count median
        center_counts.append(np.nanmedian(im[center[0] - box_size: center[0] + box_size,
                                          center[1] - box_size: center[1] + box_size]))
    return center_counts


class Image:
    last_coordinate = None
    color_image_bands = None
    color_image_radius = None
    _available_bands = None

    def _check_bands_(self, bands):
        """
        Checks if the band is a valid PS1 band

        :param bands: The input bands
        :type bands: Union
        :return:

        :raises ValueError: If one of the bands is not included in the valid PS1 bands
        """
        if bands is None:
            return self.color_image_bands
        else:
            if len(bands) != 3:
                raise ValueError("Bands must include three bands (red, green, blue).")
            for b in bands:
                if b not in self._available_bands:
                    raise ValueError(f'{b} is not a valid band. Choose one of {self._available_bands}.')

            return bands

    def _check_size(self, size):
        """
        Checks the size and returns to size in arcmin.

        :param size: The size of the image
        :type size: float, Quantity
        :return: The size of the image in arcmin
        :raises ValueError: If the unit of the size is not equivalent to degree.
        """
        if size is None:
            return self.color_image_radius
        else:
            if type(size) != u.Quantity:
                return size*u.arcmin
            else:
                if size.unit.is_equivalent(u.deg):
                    size = size.to(u.arcmin)

                    if size < 30*u.arcsec:
                        raise ValueError('Size is too small. Minimal size is 30 arcsec.')

                    if size > 30*u.arcmin:
                        raise ValueError('Size is too big. Maximal size is 30 arcmin.')

                    return size
                else:
                    raise ValueError('The unit of \'size\' must be equivalent to degrees.')

    def get_color_image(self, s, path='', bands=None, size=None):
        """
        Interface for the return of a rgb-color image from a survey

        :param s: The central coordinates of the target
        :type s: astropy.coordinates.SkyCoord
        :param path: The path to the storage place.
        :type path: str
        :param bands: The used bands, if others than the default are required
        :type bands: Union
        :param size: The size of the image, if not the default size should be used.
        :type size: float
        :return: The three different color channels
        :rtype: numpy.ndarray
        """
        raise NotImplementedError()


class SDSSImage(Image):

    def __init__(self):
        Image.__init__(self)
        self.last_coordinate = None
        self.color_image_bands = ['z', 'r', 'u']
        self.color_image_radius = 2 * u.arcmin

    def get_color_image(self, s, path='', bands=None, size=None, smooth=2):
        """
        Download the SDSS images and create an RGB image out of them.

        :param s: Coordinates of the target
        :type s: astropy.coordinates.SkyCoord
        :param path:
            The path to the save place. Default is '' which means that the image will be shown only.
        :type path: str
        :param bands: Individual band combination or None. Default is None which means that the default bands are used.
        :type bands: None, tuple, list
        :param size: The wanted size of the image
        :type size: float, astropy.units.Quantity
        :param smooth: Number of smooths
        :type smooth: int
        :return:
        """
        if bands is None:
            bands = self.color_image_bands
        sd = SDSS.get_images(s, band=bands)
        if sd is None:
            return ''
        # if SDSS took more than one image the area
        if len(sd) > 5:
            run = 0
            r_min = 999999

            n_color_bands = len(self.color_image_bands)

            # check only one image per epoch
            for i, hdu in enumerate(sd[::n_color_bands]):
                wcs = WCS(hdu[0].header)
                x, y = wcs.all_world2pix(s.ra.degree, s.dec.degree, 0)
                r = np.hypot(hdu[0].data.shape[0] // 2 - y, hdu[0].data.shape[1] // 2 - x)

                # if the current distance to the center is smaller than the previous distance to the center
                if r < r_min:
                    run = i
                    r_min = r

            # take the images from the closest epoch
            sd = sd[run * n_color_bands: (run + 1) * n_color_bands]

        imgs = []
        head = None
        wcs = None
        for hdu in sd:
            if head is None:
                head = hdu[0].header
                wcs = WCS(head)
                wcs_o = wcs
            else:
                wcs_o = WCS(hdu[0].header)

            data = smooth2d(hdu[0].data, smooth)

            # make a cut around the target
            if size is not None:
                if type(size) != u.Quantity:
                    size = size * u.arcmin
            else:
                size = self.color_image_radius

            cut = Cutout2D(data, s, size, wcs=wcs_o)
            imgs.append(cut.data)

        if '.fit' in path:
            write(path, imgs, bands)
            return

        pl.clf()
        # create a subplot and use the wcs projection
        sp = pl.subplot(projection=wcs)
        rgb = make_lupton_rgb(imgs[2], imgs[1], imgs[0], Q=13,
                              stretch=0.9, minimum=0)
        sp.imshow(rgb, origin='lower')

        sp.set_xlabel('RA')
        sp.set_ylabel('Dec')

        if path == '':
            pl.show()
        else:
            if not os.path.exists(os.path.split(path)[0]):
                os.makedirs(os.path.split(path)[0])
            pl.savefig(path)


class PanstarrsImage(Image):

    def __init__(self):
        """
        Interface to `Pan-STARRS cutout service <https://ps1images.stsci.edu/cgi-bin/ps1cutouts>`_
        """
        Image.__init__(self)
        self.last_coordinate = None
        self.color_image_bands = ['z', 'r', 'g']
        self.color_image_radius = 2 * u.arcmin
        self._available_bands = ['g', 'r', 'i', 'z', 'y']

    def get_normalized_imaged(self, s, smooth, bands=None, size=None, temp_path='./temp/'):
        """
        Returns the normalized images and the HDU's

        :param s: The central coordinates of the image
        :type s: astropy.coordinates.SkyCoord
        :param smooth: The number of linear smooths
        :type smooth: int
        :param bands:
            Optional: Set of three bands for the color image (g, r, i, z, y). Order must be red, green and blue
        :type bands: None, list, tuple
        :param size: Optional: A customized size of the image in arcmin
        :type size: float
        :param temp_path: Path to the temporary storage place. Default is './temp/.
        :type temp_path: str
        :return: The normalized images and the original HDU's
        """
        if smooth < 0 or smooth > 10:
            raise ValueError('\'smooth\' must be between 0 and 10.')

        if type(temp_path) != str:
            raise ValueError('temp_path must be a string.')

        bands = self._check_bands_(bands)

        size = self._check_size(size)

        # download Pan-STARRS images
        paths = download_all_bands(s.ra.degree, s.dec.degree, size,
                                   temp_path)
        imgs = []
        for c in bands:
            imgs.append(fits.open(paths[c]))

        # create an RGB-array
        rgb = np.zeros((imgs[2][0].shape[0], imgs[2][0].shape[1], 3))
        for i in range(3):
            data = imgs[i][0].data
            # med = np.median(data)/2
            # data -= med
            # data[data < 0] = 0
            rgb[:, :, i] = data
        # rgb[:, :, 2] = imgs[2][0].data-
        # rgb[:, :, 1] = imgs[1][0].data
        # rgb[:, :, 0] = imgs[0][0].data
        # take the center
        center = [rgb.shape[0] // 2, rgb.shape[1] // 2]

        # make a lower cut to exclude the noise
        center_counts = max_count(rgb, center)

        if 0 in center_counts:
            center_counts = max_count(rgb, center, rgb.shape[0]//2)

        # normalize to the center median counts
        for i in range(3):
            img = rgb[:, :, i]
            img /= center_counts[i]
            img[img > 100] = 100
            img -= np.nanmin(img)
            img = np.log10(smooth2d(img, smooth))
            img -= np.nanmedian(img)
            img[img < 0] = 0
            img /= np.nanmax(img)

            rgb[:, :, i] = img

        return rgb, imgs

    def get_color_image(self, s, path='', smooth=2, mark_source=False, proper_motion=None, bands=None,
                        size=None):
        """
        Download the Pan-STARRS images and create an RGB image out of them.

        :param s: Coordinates of the target
        :type s: astropy.coordinates.SkyCoord
        :param path:
            The path to the save place. Default is '' which means that the image will be shown only.
        :type path: str
        :param smooth: Number of smooth
        :type smooth: int
        :param mark_source: True if a circle should be plotted around the source, else False. Default is False.
        :type mark_source: bool
        :param proper_motion:
            The proper motion of the source or None. If the proper motion is given, an arrow will indicate
            the proper motion, if None no proper motion arrows will be drawn. Default is None.
        :type proper_motion: None, list
        :param bands: Optional argument with bands other than the defaults (z, r, g)
        :type bands: list
        :param size: The wanted size of the image
        :type size: float, astropy.units.Quantity
        :return:
        """
        temp_path = f'./temp_{np.random.randint(0, 1000000)}/'
        rgb, imgs = self.get_normalized_imaged(s, smooth, bands=bands, size=size,
                                               temp_path=temp_path)

        if '.fit' in path:
            write(path, imgs, self._check_bands_(bands))
            shutil.rmtree(temp_path)
            return

        pl.clf()
        # make a chart with a WCS projection
        sp = pl.subplot(projection=WCS(imgs[0][0].header))
        sp.imshow(rgb, origin='lower')

        # draw the source marker
        if mark_source:
            sp.scatter(s.ra.degree, s.dec.degree, 20, marker='o', c='r')

        # draw the proper motion arrow
        if proper_motion is not None:
            sp.arrow(s.ra.degree, s.dec.degree, proper_motion[0], proper_motion[1], color='r')

        pl.subplots_adjust(top=0.969, bottom=0.118,
                           left=0.146, right=0.973)

        # set the axis labels
        sp.set_xlabel('RA')
        sp.set_ylabel('Dec')

        if path == '':
            pl.show()

        else:
            if not os.path.exists(os.path.split(path)[0]):
                os.makedirs(os.path.split(path)[0])
            pl.savefig(path)
        for i in imgs:
            i.close()
        shutil.rmtree(temp_path)
