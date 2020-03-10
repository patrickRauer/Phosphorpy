from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astroquery.sdss import SDSS as sdss
from astroquery.gama import GAMA as GAMA_INTERFACE
from astroquery.vizier import Vizier
import numpy as np
import urllib
import os
import warnings
from glob import glob

from Phosphorpy.config.survey_data import get_survey_data
from Phosphorpy.data.sub.spectra import Spectra, SpectraList


SDSS = 'sdss'
LAMOST = 'lamost'
GAMA = 'gama'
DFGRS = '2fDGR'


def _remove_temp_files():
    """
    Removes all temporary spectra
    :return:
    """
    for f in glob('temp_*_spec_*.fits'):
        os.remove(f)


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
        ids = np.arange(len(coord))+1
    else:
        if len(ids) != len(coord):
            raise ValueError('If ID\'s are given, they must have the same length as the given coordinates')
        elif type(ids) != np.ndarray:
            ids = np.array(ids)
    return ids


def get_lamost_spectra(coord, ids=None):
    """
    Downloads spectra from the LAMOST DR5.
    The cross-match radius is set to 1 arcsec.

    :param coord: The coordinates of the targets
    :type coord: SkyCoord
    :param ids:
        ID's corresponding to the coordinates. Default is None, which means that the positions in the coord
        is used as ID.
    :return: The downloaded spectra
    :rtype: SpectraList
    """
    coord = _check_coordinates(coord)

    ids = _check_ids(ids, coord)

    spec_list = SpectraList()

    lamost_data = get_survey_data('LAMOST')
    lamost = Vizier(columns=lamost_data['columns'])
    lamost.ROW_LIMIT = -1
    rs = lamost.query_region(coord, 1 * u.arcsec, catalog=lamost_data['vizier'])

    if len(rs) == 0:
        return spec_list

    rs = rs[0]

    lamost_coord = SkyCoord(rs['RAJ2000'],
                            rs['DEJ2000'])
    sort = lamost_coord.match_to_catalog_sky(coord)[0]
    ids = ids[sort]

    # create a temporary path with a random number at the end to avoid potential overwriting
    temp_path = f'temp_lamost_spec_{np.random.randint(0, 10000)}.fits'

    for obs_id, index in zip(rs['ObsID'], ids):
        urllib.request.urlretrieve(lamost_data['url'].format(obs_id), temp_path)
        with fits.open(temp_path) as fi:
            d = fi[0].data
            wave = d[2]
            fl = d[0]

            spec = Spectra(wavelength=wave,
                           flux=fl,
                           wavelength_unit=u.angstrom,
                           survey=LAMOST,
                           meta=fi[0].header)
            spec_list.append(spec, index)

    os.remove(temp_path)

    return spec_list


def get_sdss_spectra(coord, ids=None):
    """
    Downloads SDSS spectra for objects at the given coordinates, if spectra are available.

    :param coord: The coordinates of the objects for which spectra are wanted
    :type coord: SkyCoord
    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [1, len(coord)].
        Default is None.
    :type ids: Union
    :return: SpectraList with the found spectra or an empty SpectraList, if no spectra was found.
    :rtype: SpectraList
    """
    coord = _check_coordinates(coord)

    ids = _check_ids(ids, coord)

    spec_list = SpectraList()
    try:
        rs = sdss.query_region(coord, spectro=True)
    except urllib.error.HTTPError:
        warnings.warn('Connection problem to SDSS')
        rs = None

    # if no SDSS spectra was found
    if rs is None:
        return spec_list
    # convert the output coordinates to SkyCoord and get the order of the input coordinate in
    # respect to tine output coordinates
    sdss_coord = SkyCoord(rs['ra']*u.deg,
                          rs['dec']*u.deg)
    sort = sdss_coord.match_to_catalog_sky(coord)[0]
    ids = ids[sort]

    try:
        # download the SDSS spectra
        sp = sdss.get_spectra(matches=rs)
    except urllib.error.HTTPError:
        warnings.warn('Connection problem to SDSS')
        return spec_list

    # read the spectra and convert the wavelength to a linear scale
    # add the resulting Spectra object to the SpectraList
    for spec, index in zip(sp, ids):
        spec = Table(spec[1].data)
        spec['wavelength'] = np.power(10., spec['loglam'])
        # create a new Spectra object for the spectra with the wavelength and the flux.
        # Use angstrom as default wavelength units
        spec = Spectra(wavelength=spec['wavelength'],
                       flux=spec['flux'],
                       wavelength_unit=u.angstrom,
                       survey=SDSS,
                       meta=spec[1].header)
        spec_list.append(spec, index)

    return spec_list


def get_gama_spectra(coord, ids=None):
    """
    Downloads GAMA spectra for objects at the given coordinates, if spectra are available.

    :param coord: The coordinates of the objects for which spectra are wanted
    :type coord: SkyCoord
    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [0, len(coord)-1].
        Default is None.
    :type ids: Union
    :return: SpectraList with the found spectra or an empty SpectraList, if no spectra was found.
    :rtype: SpectraList
    """
    spec_list_gama = SpectraList()
    spec_list_2dfgrs = SpectraList()
    coord = _check_coordinates(coord)

    ids = _check_ids(ids, coord)

    radius = 2./3600
    sql = 'SELECT RA, DEC, URL FROM SpecAll WHERE '
    sql += ' OR '.join(
        [
            f'(RA BETWEEN {ra-radius} AND {ra+radius} AND DEC BETWEEN {dec-radius} AND {dec+radius})'
            for ra, dec in zip(coord.ra.degree, coord.dec.degree)
        ]
    )
    try:
        rs = GAMA_INTERFACE.query_sql(sql)
    except ValueError:
        return spec_list_gama, spec_list_2dfgrs

    gama_coord = SkyCoord(rs['RA']*u.deg, rs['DEC']*u.deg)
    sort = gama_coord.match_to_catalog_sky(coord)[0]
    ids = ids[sort]

    # create a temporary path with a random number at the end to avoid potential overwriting
    temp_path = f'temp_gama_spec_{np.random.randint(0, 10000)}.fits'
    for r, index in zip(rs, ids):
        urllib.request.urlretrieve(r['URL'], temp_path)
        with fits.open(temp_path) as fi:
            # if it is only a PrimaryHDU, then it is 2dFGRS spectra
            if len(fi) == 1:
                header = fi[0].header
                data = fi[0].data
                wavelength = header['CDELT1']*(np.arange(header['NAXIS1'])-header['CRPIX1'])+header['CRVAL1']
                spec = Spectra(wavelength=wavelength,
                               flux=data[0], wavelength_unit=u.angstrom,
                               survey=DFGRS,
                               meta=header)
                spec_list_2dfgrs.append(spec, index)
            else:
                wcs = WCS(fi[0].header)
                waves = wcs.all_pix2world(np.arange(len(fi[0].data[0])), 0, 0)[0]
                waves = np.power(10., waves)
                spec = Spectra(wavelength=waves,
                               flux=fi[0].data[0],
                               wavelength_unit=u.angstrom,
                               survey=GAMA,
                               meta=fi[0].header)
                spec_list_gama.append(spec, index)
    # remove the temporary file
    os.remove(temp_path)
    return spec_list_gama, spec_list_2dfgrs


def get_spectra(coord, ids=None, source='SDSS'):
    """
    Search for spectra of sources at the coordinate(s) in the specific survey

    :param coord: The coordinate(s) of the required sources
    :type coord: SkyCoord
    :param ids:
        ID's of the spectra. If no ID's are given, the ID's are set to [1, len(coord)].
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


def get_all_spectra(coordinates, index=None):
    """
    Downloads all available spectra from SDSS, GAMA and LAMOST DR5

    :param coordinates: The coordinates of the targets with potential spectra
    :type coordinates: SkyCoord
    :param index:
        Optional, individual indices of the targets. If no indices are given, the return spectra will have indices
        from 1 to len(coordinates) with the same order as the given coordinates.
        Default is None.
    :type index: Union
    :return: The found spectra from the three surveys
    :rtype: SpectraList
    """
    specs = None
    try:
        specs = get_spectra(coordinates, index)
    except:
        pass

    try:
        lamost = get_spectra(coordinates, index, LAMOST)
        if specs is None:
            specs = lamost
        else:
            specs.merge(lamost)
    except:
        pass

    try:
        gama = get_spectra(coordinates, index, GAMA)
        if specs is None:
            specs = gama[0]
            specs.merge(gama[1])
        else:
            specs.merge(gama[0])
            specs.merge(gama[1])
    except:
        pass

    _remove_temp_files()

    return specs
