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

from Phosphorpy.external.panstarrs import download_all_bands


def smooth2d(mat, c=5):
    """
    Smooths a 2d-array

    :param mat: The input data
    :type mat: np.ndarray
    :param c: The number of smooths
    :param c: int
    :return: The c-times smoothed input data
    """
    if c == 0:
        return mat
    out = np.zeros(mat.shape)

    out[0, 0] = (2*mat[0, 0]+mat[1, 0]+mat[0, 1])/4
    out[0, -1] = (2*mat[0, -1]+mat[1, -1]+mat[0, -2])/4
    out[-1, 0] = (2*mat[-1, 0]+mat[-2, 0]+mat[-1, 1])/4
    out[-1, -1] = (2*mat[-1, -1]+mat[-2, -1]+mat[-1, -2])/4

    out[1:-1, 0] = (3*mat[1:-1, 0]+mat[1:-1, 1]+mat[:-2, 0]+mat[2:, 0])/6
    out[1:-1, -1] = (3*mat[1:-1, -1]+mat[1:-1, -2]+mat[:-2, -1]+mat[2:, -1])/6
    out[0, 1:-1] = (3*mat[0, 1:-1]+mat[1, 1:-1]+mat[0, :-2]+mat[0, 2:])/6
    out[-1, 1:-1] = (3*mat[-1, 1:-1]+mat[-2, 1:-1]+mat[-1, :-2]+mat[-1, 2:])/6

    out[1:-1, 1:-1] = (4*mat[1:-1, 1:-1]+mat[:-2, 1:-1]+mat[2:, 1:-1]+mat[1:-1, :-2]+mat[1:-1, 2:])/8
    return smooth2d(out, c=c-1)


class SDSSImage:
    last_coordinate = None
    color_image_bands = ['z', 'r', 'u']
    color_image_radius = 2 * u.arcmin

    def __init__(self):
        pass

    def get_color_image(self, s, path='', bands=None):
        """
        Download the SDSS images and create an RGB image out of them.

        :param s: Coordinates of the target
        :type s: astropy.coordinates.SkyCoord
        :param path:
            The path to the save place. Default is '' which means that the image will be shown only.
        :type path: str
        :param bands: Individual band combination or None. Default is None which means that the default bands are used.
        :type bands: None, tuple, list
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
            # print(hdu[0].header)
            if head is None:
                head = hdu[0].header
                wcs = WCS(head)
                wcs_o = wcs
            else:
                wcs_o = WCS(hdu[0].header)

            data = smooth2d(hdu[0].data)

            # make a cut around the target
            cut = Cutout2D(data, s, self.color_image_radius, wcs=wcs_o)
            imgs.append(cut.data)

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


class PanstarrsImage:
    last_coordinate = None
    color_image_bands = ['z', 'r', 'g']
    color_image_radius = 2 * u.arcmin

    def __init__(self):
        pass

    def get_normalized_imaged(self, s, smooth, bands=None):
        """
        Returns the normalized images and the HDU's

        :param s: The central coordinates of the image
        :type s: astropy.coordinates.SkyCoord
        :param smooth: The number of linear smooths
        :type smooth: int
        :param bands:
            Optional: Set of three bands for the color image (g, r, i, z, y). Order must be red, green and blue
        :type bands: None, list, tuple
        :return: The normalized images and the original HDU's
        """
        if bands is None:
            bands = self.color_image_bands
        else:
            if len(bands) != 3:
                raise ValueError("Bands must include three bands")
        # download Pan-STARRS images
        paths = download_all_bands(s.ra.degree, s.dec.degree, self.color_image_radius,
                                   './temp/')
        imgs = []
        for c in bands:
            imgs.append(fits.open(paths[c]))

        # create an RGB-array
        rgb = np.zeros((imgs[2][0].shape[0], imgs[2][0].shape[1], 3))
        rgb[:, :, 2] = imgs[2][0].data
        rgb[:, :, 1] = imgs[1][0].data
        rgb[:, :, 0] = imgs[0][0].data
        rgb = np.log10(rgb)

        # take the center
        center = [rgb.shape[0] // 2, rgb.shape[1] // 2]
        center_counts = []

        rgb = np.nan_to_num(rgb)
        # make a lower cut to exclude the noise
        for i in range(3):
            im = rgb[:, :, i]

            med = 2*np.median(np.nanmean(im, axis=0))

            im[im < med] = med
            im -= med

            # collect the center count median
            center_counts.append(np.nanmedian(im[center[0] - 2: center[0] + 2,
                                              center[1] - 2: center[1] + 2]))

        # normalize to the center median counts
        for i in range(3):
            rgb[:, :, i] /= center_counts[i]
            rgb[:, :, i] = smooth2d(rgb[:, :, i], smooth)
        return rgb, imgs

    def get_color_image(self, s, path='', smooth=2, mark_source=False, proper_motion=None):
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
        :return:
        """
        rgb, imgs = self.get_normalized_imaged(s, smooth)

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
