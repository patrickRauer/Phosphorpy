from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.sdss import SDSS as sdss
import numpy as np
import math

from Phosphorpy.external.xmatch import xmatch
from Phosphorpy.data.sub.spectra import Spectra, SpectraList


SDSS = 'sdss'
LAMOST = 'lamost'
GAMA = 'gama'


def get_lamost_spectra(coord):
    download = 'http://dr4.lamost.org/./spectrum/fits/{}?token='


def get_sdss_spectra(coord, ids=None):

    if coord.isscalar:
        coord = SkyCoord(
            np.array([coord.ra.degree])*coord.ra.unit,
            np.array([coord.dec.degree])*coord.dec.unit
        )

    if ids is None:
        ids = np.arange(len(coord))

    spec_list = SpectraList()
    rs = sdss.query_region(coord, spectro=True)

    # if no SDSS spectra was found
    if rs is None:
        return spec_list

    sdss_coord = SkyCoord(rs['ra']*u.deg,
                          rs['dec']*u.deg)
    print(sdss_coord.match_to_catalog_sky(coord))
    sp = sdss.get_spectra(matches=rs)

    spec_count = 0
    for c, index in zip(coord, ids):
        distance = math.hypot(c.ra.degree - rs['ra'][spec_count],
                              c.dec.degree - rs['dec'][spec_count])*3600
        print(distance)
        if distance < 2.5:
            spec = Table(sp[spec_count][1].data)
            spec['wavelength'] = np.power(10, spec['loglam'])
            spec = Spectra(wavelength=spec['wavelength'],
                           flux=spec['flux'])
            spec_list.append(spec, index)
            spec_count += 1

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
        print(coord.isscalar)

    if source == SDSS:
        return get_sdss_spectra(coord)
    elif source == LAMOST:
        return get_lamost_spectra(coord)
    elif source == GAMA:
        return get_gama_spectra(coord)
    else:
        raise ValueError(f'{source} is unknown for spectra.')
