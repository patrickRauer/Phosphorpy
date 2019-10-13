from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.sdss import SDSS as sdss
import numpy as np

from Phosphorpy.external.xmatch import xmatch
from Phosphorpy.data.sub.spectra import Spectra, SpectraList


SDSS = 'sdss'
LAMOST = 'lamost'
GAMA = 'gama'


def get_lamost_spectra(coord):
    download = 'http://dr4.lamost.org/./spectrum/fits/{}?token='


def get_sdss_spectra(coord, ids=None):
    """
    Downloads SDSS spectra for objects at the given coordinates, if spectra are available.

    :param coord: The coordinates of the objects for which spectra are wanted
    :type coord: SkyCoord
    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [0, len(coord)-1].
        Default is None.
    :type ids: Union
    :return: SpectraList with the found spectra or an empty SpectraList, if no spectra was found.
    :rtype: SpectraList
    """
    if coord.isscalar:
        coord = SkyCoord(
            np.array([coord.ra.degree])*coord.ra.unit,
            np.array([coord.dec.degree])*coord.dec.unit
        )

    if ids is None:
        ids = np.arange(len(coord))
    else:
        if len(ids) != len(coord):
            raise ValueError('If ID\'s are given, they must have the same length as the given coordinates')

    spec_list = SpectraList()
    rs = sdss.query_region(coord, spectro=True)

    # if no SDSS spectra was found
    if rs is None:
        return spec_list

    # convert the output coordinates to SkyCoord and get the order of the input coordinate in
    # respect to tine output coordinates
    sdss_coord = SkyCoord(rs['ra']*u.deg,
                          rs['dec']*u.deg)
    sort = sdss_coord.match_to_catalog_sky(coord)[0]
    ids = ids[sort]

    # download the SDSS spectra
    sp = sdss.get_spectra(matches=rs)

    # read the spectra and convert the wavelength to a linear scale
    # add the resulting Spectra object to the SpectraList
    for spec, index in zip(sp, ids):
        spec = Table(spec[1].data)
        spec['wavelength'] = np.power(10., spec['loglam'])
        # create a new Spectra object for the spectra with the wavelength and the flux.
        # Use angstrom as default wavelength units
        spec = Spectra(wavelength=spec['wavelength'],
                       flux=spec['flux'],
                       wavelength_unit=u.angstrom)
        spec_list.append(spec, index)

    return spec_list


def get_gama_spectra(coord):
    pass


def get_spectra(coord, source='SDSS'):
    """
    Search for spectra of sources at the coordinate(s) in the specific survey

    :param coord: The coordinate(s) of the required sources
    :type coord: SkyCoord
    :param source: The source of the spectra
    :type source: str
    :return: The spectra, if any spectra was found.
    """
    source = source.lower()

    if coord.isscalar:
        coord = SkyCoord(
            np.array([coord.ra.degree])*coord.ra.unit,
            np.array([coord.dec.degree])*coord.dec.unit
        )

    if source == SDSS:
        return get_sdss_spectra(coord)
    elif source == LAMOST:
        return get_lamost_spectra(coord)
    elif source == GAMA:
        return get_gama_spectra(coord)
    else:
        raise ValueError(f'{source} is unknown for spectra.')
