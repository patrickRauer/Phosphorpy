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


def _check_coordinates(coord):
    """
    Checks if the coordinates are scalar and if this is true, convert them to non scalar
    :param coord: The coordinates
    :type coord: SkyCoord
    :return: The input coordinates as non scalar
    :rtype: SkyCoord
    """
    if coord.isscalar:
        coord = SkyCoord(
            np.array([coord.ra.degree])*coord.ra.unit,
            np.array([coord.dec.degree])*coord.dec.unit
        )
    return coord


def _check_ids(ids, coord):
    """
    Checks, if the ids are None and if this is true, create a default ID list.
    Checks also, if the length of the ids is the same as the length coordinates.

    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [0, len(coord)-1].
        Default is None.
    :type ids: Union
    :param coord: The coordinates of the objects for which spectra are wanted
    :type coord: SkyCoord
    :return: The list of ID's
    :rtype: ndarray
    """
    if ids is None:
        ids = np.arange(len(coord))
    else:
        if len(ids) != len(coord):
            raise ValueError('If ID\'s are given, they must have the same length as the given coordinates')
        elif type(ids) != np.ndarray:
            ids = np.array(ids)
    return ids


def get_lamost_spectra(coord, ids=None):
    coord = _check_coordinates(coord)

    ids = _check_ids(ids, coord)

    spec_list = SpectraList()

    lamost_download_url = 'http://dr4.lamost.org/./spectrum/fits/{}?token='
    lamost = Vizier(
        columns=['_q', 'RAJ2000', 'DEJ2000', 'ObsID', 'snru', 'snrg', 'snrr', 'snri', 'snrz', 'z', 'SubClass'])
    lamost.ROW_LIMIT = -1
    rs = lamost.query_region(coord, 1 * u.arcsec, catalog='V/153/dr4')[0]

    if len(rs) == 0:
        return spec_list
    
    lamost_coord = SkyCoord(rs['ra']*u.deg,
                            rs['dec']*u.deg)
    sort = lamost_coord.match_to_catalog_sky(coord)[0]
    ids = ids[sort]

    # create a temporary path with a random number at the end to avoid potential overwriting
    temp_path = f'temp_lamost_spec_{np.random.randint(0, 10000)}.fits'

    for obs_id, index in zip(rs['ObsID'], ids):
        urllib.request.urlretrieve(lamost_download_url.format(obs_id), temp_path)
        with fits.open(temp_path) as fi:
            d = fi[0].data
            wave = d[2]
            fl = d[0]

            spec = Spectra(wavelength=wave,
                           flux=fl,
                           wavelength_unit=u.angstrom)
            spec_list.append(spec, index)
    return spec_list


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
    coord = _check_coordinates(coord)

    ids = _check_ids(ids, coord)

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


def get_gama_spectra(coord, ids=None):
    pass


def get_spectra(coord, ids=None, source='SDSS'):
    """
    Search for spectra of sources at the coordinate(s) in the specific survey

    :param coord: The coordinate(s) of the required sources
    :type coord: SkyCoord
    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [0, len(coord)-1].
        Default is None.
    :type ids: Union
    :param source: The source of the spectra
    :type source: str
    :return: The spectra, if any spectra was found.
    """
    source = source.lower()

    coord = _check_coordinates(coord)
    ids = _check_ids(ids, coord)

    if source == SDSS:
        return get_sdss_spectra(coord, ids=ids)
    elif source == LAMOST:
        return get_lamost_spectra(coord, ids=ids)
    elif source == GAMA:
        return get_gama_spectra(coord, ids=ids)
    else:
        raise ValueError(f'{source} is unknown for spectra.')
