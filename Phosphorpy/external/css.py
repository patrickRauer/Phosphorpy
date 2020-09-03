"""
This script provides an interface to the CSS-server to download light curves from it.

(http://nesssi.cacr.caltech.edu/DataRelease/)

@author: Jean Patrick Rauer
"""
from bs4 import BeautifulSoup
from requests.exceptions import ConnectionError
import requests
import pylab as pl
import pandas as pd
import numpy as np
import os
import warnings


LIGHT_CURVES_URL = 'http://nesssi.cacr.caltech.edu/cgi-bin/getmulticonedb_release2.cgi'


def smooth(d, c=5):
    """
    Smooth a 1D data set linearly

    :param d: The input data
    :type d: np.ndarray
    :param c: The number of smooth
    :param c: int
    :return: The c-times smoothed input data
    :rtype: np.ndarray
    """
    if c == 0:
        return d
    x = np.zeros(len(d))
    x[0] = (d[0]+d[1])/2
    x[1:-1] = (2*d[1:-1]+d[2:]+d[:-2])/4
    x[-1] = (d[-2]+d[-1])/2
    return smooth(x, c=c-1)


def smooth_err(d, err, c=5):
    """
    Smooth a 1D data set with respect to the errors

    :param d: The input data
    :type d: np.ndarray
    :param err: The errors of the data
    :type err: np.ndarray
    :param c: The number of smooths
    :type c: int
    :return: The c-times smoothed input data
    :rtype: np.ndarray
    """
    if c == 0:
        return d
    e = 1/err
    x = np.zeros(len(d))
    x[0] = (d[0]*e[0]+d[1]*e[1])/(e[0]+e[1])
    x[1:-1] = (2*d[1:-1]*e[1:-1]+d[2:]*e[2:]+d[:-2]*e[:-2])/(2*e[1:-1]+e[2:]+e[:-2])
    x[-1] = (d[-2]*e[-2]+d[-1]*e[-1])/(e[-2]+e[-1])
    return smooth_err(x, err, c=c-1)


def vari_index(d, err):
    """
    Computes the variability index

    :param d: The input data
    :type d: np.ndarray
    :param err: The errors of the input data
    :type err: np.ndarray
    :return: The variability index
    :rtype: float
    """
    return 1/(len(d)-1)*np.sum(np.square(d-np.mean(d))/err**2)


def download_light_curve(ra, dec):
    """
    Downloads the light curve of the target from the Catalina Sky Survey
    
    :param ra: RA coordinate in degree
    :type ra: float
    :param dec: Dec coordinate in degree
    :type dec: float
    :returns: The light curve of the target
    :rtype: pandas.DataFrame
    """
    try:
        r = requests.post('http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi',
                          data={'RA': ra, 'Dec': dec, 'DB': 'photcat', 'OUT': 'csv', 'SHORT': 'short'})

        # because the CSS HTML code has 'html' and 'HTML' as tags, BeautifulSoup takes only
        # the first (html) section, therefore we split the receiving HTML text into two parts.

        soup = BeautifulSoup(r.text.split('</html>')[-1])
        lc = pd.read_csv(soup.find_all('a', href=True)[0]['href'])
        return lc
    except ValueError:
        raise ValueError('No light curve available.')


def _download_light_curve_parts(temp_path='temp.txt'):
    with open(temp_path) as f:
        r = requests.post(LIGHT_CURVES_URL,
                          data={'DB': 'photcat', 'OUT': 'csv',
                                'SHORT': 'short'},
                          files={'upload_file': f})
        soup = BeautifulSoup(r.text)
        download_url = soup.find('a', href=True)['href']

        d = pd.read_csv(download_url)

        return d


def download_light_curves(ra, dec):
    """
    Downloads a set of CSS light curves from the CSS server

    :param ra: The RA coordinates
    :type ra: list
    :param dec: The Dec coordinates
    :type dec: list
    :return: The downloaded light curves
    :rtype: pandas.DataFrame
    """
    results = []
    try:
        part_size = 80
        parts = len(ra) // part_size + 1
        coords = pd.DataFrame({'ra': np.round(ra, 5), 'dec': np.round(dec, 5)})
        for i in range(parts):

            coords[i * part_size: (i + 1) * part_size].to_csv('temp.txt', sep='\t', header=False)

            d = _download_light_curve_parts()

            results.append(d)
    except ConnectionError:
        warnings.warn('Connection problem')

    finally:
        os.remove('temp.txt')

    try:
        return pd.concat(results)
    except ValueError:
        return None


def daily_average(d):
    """
    Takes the daily average of the light curve to reduce the noise.
    """
    d = d.copy()
    d['MJD_day'] = np.int32(d['MJD'].values)
    d = d.groupby('MJD_day')
    d = d.aggregate(np.mean)
    return d


def plot_light_curve(ra, dec):
    """
    PLots the CSS light curve of the target
    
    :param ra: RA coordinate in degree
    :type ra: float
    :param dec: Dec coordinate in degree
    :type dec: float
    """
    lc = download_light_curve(ra, dec)
    
    pl.clf()
    sp = pl.subplot(211)
    avg = daily_average(lc)
    sp.errorbar(avg['MJD'],
                avg['Mag'],
                avg['Magerr'],
                fmt='.k',
                capsize=2)
    sp.invert_yaxis()
    sp = pl.subplot(212)
    avg_smooth = smooth_err(avg['Mag'].values, avg['Magerr'].values)
    sp.scatter(avg['MJD'],
               avg_smooth,
               marker='.',
               c='k')
    
    sp.set_xlabel('MJD')
    sp.set_ylabel('mag')
    sp.invert_yaxis()
    pl.show()
