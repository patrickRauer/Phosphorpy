from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.irsa import Irsa
import pandas as pd
import numpy as np

ZWICKY_URL = 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?' \
             'POS=CIRCLE{:+09.4f}{:+09.4f}{:+07.4f}&BANDNAME=r&NOBS_MIN=3&BAD_CATFLAGS_MASK=32768&FORMAT=ipac_table'


def progress_bar(iteration, total, prefix='',
                 suffix='', length=100, fill='█', print_end="\r"):
    percent = f'{100*iteration/total:04.1f}'
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()


def download_light_curve(ra, dec, radius=None, index=None):
    """
    Downloads the light curve from the `Zwicky-survey <https://www.ztf.caltech.edu/>`_ DR2

    :param ra: RA coordinate
    :type ra: float, Union
    :param dec: Dec coordinate
    :type dec: float, Union
    :param radius:
        The search radius. If no units are given, the search radius unit is set to arcsec.
        Default is 2 arcsec.
    :type radius: float, Quantity
    :param index:
        The index of the coordinates, if more than one coordinate set is given.
        Otherwise index will be ignored.
    :return: The light curve(s) from PTF
    :rtype: DataFrame
    """
    if radius is None:
        radius = 0.0014

    if type(radius) == u.Quantity:
        radius = radius.to(u.deg)
        radius = radius.value

    if not isinstance(ra, np.float):
        if index is None:
            index = np.arange(len(ra))
        out = []
        counter = 0
        print('Download ZTF light curves')
        for r, d, i in zip(ra, dec, index):
            progress_bar(counter, len(ra), 'Download', '')
            try:
                t = download_light_curve(r, d, radius)
                t['row_id'] = i
                out.append(t)
            except:
                pass
            counter += 1
        progress_bar(counter, len(ra), 'Download', 'Complete')
        return pd.concat(out)

    t = Table.read(ZWICKY_URL.format(ra, dec, radius), format='ascii.ipac')

    return t[['mjd', 'oid', 'mag', 'magerr', 'ra', 'dec']].to_pandas()


def download_ptf(ra, dec, radius=None, index=None):
    """
    Downloads the `PTF <https://www.ptf.caltech.edu/>`_ light curve and reformat the columns
    to align the names to match the naming of the other light curves

    :param ra: RA coordinate
    :type ra: float, Union
    :param dec: Dec coordinate
    :type dec: float, Union
    :param radius:
        The search radius. If no units are given, the search radius unit is set to arcsec.
        Default is 2 arcsec.
    :type radius: float, Quantity
    :param index:
        The index of the coordinates, if more than one coordinate set is given.
        Otherwise index will be ignored.
    :return: The light curve(s) from PTF
    :rtype: DataFrame
    """
    # if no radius is given
    if radius is None:
        radius = 2*u.arcsec
    # if the radius is just a float
    elif type(radius) == float:
        radius = radius*u.arcsec

    if not isinstance(ra, np.float):
        if index is None:
            index = np.arange(len(ra))
        out = []
        counter = 0
        print('Download PTF light curves')
        for r, d, i in zip(ra, dec, index):
            progress_bar(counter, len(ra), 'Download', 'Complete')
            t = download_ptf(r, d, radius, i)
            t['row_id'] = i
            out.append(t)
            counter += 1
        progress_bar(counter, len(ra), 'Download', 'Complete')
        return pd.concat(out)
    s = SkyCoord(ra*u.deg, dec*u.deg)
    table = Irsa.query_region(s, catalog="ptf_lightcurves", spatial="Cone",
                              radius=radius)
    table = table[['obsmjd', 'mag_autocorr', 'magerr_auto', 'ra', 'dec', 'oid']].to_pandas()
    table = table.rename({'obsmjd': 'mjd', 'mag_autocorr': 'mag',
                          'magerr_auto': 'magerr'}, axis='columns')
    return table
