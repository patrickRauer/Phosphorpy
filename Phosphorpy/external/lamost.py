from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u
from astropy import constants as cons
import numpy as np
import urllib
import numbers
import os

from Phosphorpy.data.sub.spectra import Spectra, SpectraList


class LamostSpectraList(SpectraList):

    def __init__(self, spectra=None):
        SpectraList.__init__(self, spectra)

    def as_dataframe(self, as_velocity=False, redo=False):
        """
        Returns the main information of the stored spectra as a pandas DataFrame

        :param as_velocity:
        :param redo:
        :return: The main information as a dataframe
        :rtype: pandas.DataFrame
        """
        out = []
        for s, d_id in zip(self._spectra, self._ids):
            properties = s.estimate_line_properties(as_velocity=as_velocity,
                                                    redo=redo)
            properties['obsID'] = s.obs_id
            properties['ID'] = d_id
            out.append(properties)
        out = vstack(out).to_pandas()
        return out.set_index('ID')


class LamostSpectra(Spectra):
    obs_id = None
    _wavelength = None
    _flux = None
    header = None
    _lines = None
    estimations = None

    def __init__(self, obsid):
        """
        LamostSpectra provides an interface to work with LAMOST spectra.

        :param obsid: The observation ID of the spectra
        :type obsid: int
        """
        Spectra.__init__(self)
        self.obs_id = obsid
        data, header = download_spectra(obsid)
        self._wavelength = data[2]
        self._flux = data[0]
        self.header = header

    @property
    def radial_velocity(self):
        z = float(self.header['z'])
        if z == -9999:
            z = np.nan
        return (1-z)*cons.c.to(u.km/u.s)

    @property
    def seeing(self):
        return self.header['SEEING']

    @staticmethod
    def from_coordinates(ra, dec, ids=None):
        """
        Queries Lamost spectra from the coordinate(s) and returns a LamostSpectra object in the case
        of a single result and a list of spectra in the case of multiple spectra.

        :param ra: The RA coordinate
        :type ra: float, Union
        :param dec: The Dec coordinate
        :type dec: float, Union
        :param ids: List with ids. ids will be ignored, if ra and dec are floats.
        :type ids: Union
        :return: The results or None, if no spectra is found
        :rtype: LamostSpectra, list
        """
        spectra_list = get_spectra_list_ra_dec(ra, dec)
        if spectra_list is None:
            return None
        elif len(spectra_list) == 1:
            return LamostSpectra(spectra_list['ObsID'][0])
        else:
            out = LamostSpectraList()
            for s in spectra_list:
                if ids is None:
                    out.append(LamostSpectra(s['ObsID']))
                else:
                    out.append(LamostSpectra(s['ObsID']), ids[s['_q']-1])
            return out


def _convert_to_array(x):
    """
    Checks if the input is an numpy array. If not convert it to it
    :param x:
    :return: The input data as numpy array
    :rtype: np.ndarray
    """
    if type(x) != np.ndarray:
        if isinstance(x, numbers.Number):
            x = (x, )
        x = np.array(x)
    return x


def get_spectra_list_ra_dec(ra, dec):
    """
    Starts a query of the LAMOST DR4 catalog and searches for spectra.

    :param ra: RA of the target(s)
    :type ra: float, list, np.ndarray
    :param dec: Dec of the target(s)
    :type dec: float, list, np.ndarray
    :return: A table with the results or None, if no data are found.
    :rtype: Table
    """
    ra = _convert_to_array(ra)
    dec = _convert_to_array(dec)
    return get_spectra_list(SkyCoord(ra*u.deg, dec*u.deg))


def get_spectra_list(coordinates):
    """
    Starts a query of the LAMOST DR4 catalog and searches for spectra.

    :param coordinates: Coordinates of the target(s)
    :type coordinates: SkyCoord
    :return: A table with the results or None, if no data are found.
    :rtype: Table
    """
    if type(coordinates) != SkyCoord:
        raise TypeError('Coordinates must be a astropy.coordinates.SkyCoord object.')

    lamost = Vizier(columns=['_q', 'RAJ2000', 'DEJ2000', 'snru', 'snrg', 'snrr', 'snri', 'snrz',
                             'ObsID', 'z', 'SubClass', 'class'])
    lamost.ROW_LIMIT = -1

    try:
        rs = lamost.query_region(coordinates, 1*u.arcsec, catalog='V/153/dr4')[0]

        return rs
    except IndexError:
        return None


def download_spectra(obsid):
    """
    Downloads a spectra from the LAMOST server

    :param obsid: The observation ID of the spectra.
    :type obsid: int
    :return: The data of the spectra and the header of it
    :rtype: np.ndarray, fits.Header
    """
    temp_path = './temp.fits'
    url = 'http://dr4.lamost.org/./spectrum/fits/{}?token='
    urllib.request.urlretrieve(url.format(obsid),
                               temp_path)
    with fits.open(temp_path) as fi:
        data, header = fi[0].data, fi[0].header

    os.remove(temp_path)
    return data, header
