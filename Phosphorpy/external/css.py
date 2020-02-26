"""
This script provides an interface to the CSS-server to download light curves from it.

(http://nesssi.cacr.caltech.edu/DataRelease/)

@author: Jean Patrick Rauer
"""

import requests
import pylab as pl
import urllib
import pandas
import numpy as np
import os
import warnings


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
        url = r.text.split('save as to')[-1]
        url = url.split('href=')[-1]
        url = url.split('>download')[0]
        urllib.request.urlretrieve(url, "temp.csv")
        pd = pandas.read_csv('temp.csv')
        os.remove('temp.csv')
        return pd
    except ValueError:
        raise ValueError('No light curve available.')


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
    url = 'http://nesssi.cacr.caltech.edu/cgi-bin/getmulticonedb_release2.cgi'
    url2 = 'http://nesssi.cacr.caltech.edu/cgi-bin/getdatamulti.cgi?ID={}&txtInput=0000'
    results = []
#    try:
    try:
        part_size = 80
        parts = len(ra)//part_size+1
        for i in range(parts):
            d_ra = ra[i*part_size: (i+1)*part_size]
            d_dec = dec[i*part_size: (i+1)*part_size]
            with open('temp.txt', 'w') as f:
                for j, (r, d) in enumerate(zip(d_ra, d_dec)):
                    row_id = i*part_size+j
                    f.write(f"{row_id+1}\t{round(r, 5)}\t{round(d, 5)}\n")

            with open('temp.txt') as f:
                r = requests.post(url,
                                  data={'DB': 'photcat', 'OUT': 'csv', 
                                        'SHORT': 'short'},
                                  files={'upload_file': f})

            r = r.text.split('location.href=\'')[-1]
            r = r.split('\'')[0]

            if 'Query service results' in r:
                r = url2.format(
                    r.split('name="ID" value="')[-1].split('"')[0]
                )
                r = requests.get(r).text

                r = r.split('href=')[-1].split('>download')[0]

            urllib.request.urlretrieve(r, "temp.csv")
            pd = pandas.read_csv('temp.csv')
            os.remove('temp.csv')
            os.remove('temp.txt')
            results.append(pd)
    except:
        warnings.warn('Connection problem')

#    except ValueError:
#        raise ValueError('No light curve available.')
    return pandas.concat(results)


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
