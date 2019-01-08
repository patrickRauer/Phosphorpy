#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 18:07:21 2018

@author: patrickr
"""

import requests
import pylab as pl
import urllib
import pandas
import numpy as np
import os


def smooth(d, c=5):
    if c == 0:
        return d
    x = np.zeros(len(d))
    x[0] = (d[0]+d[1])/2
    x[1:-1] = (2*d[1:-1]+d[2:]+d[:-2])/4
    x[-1] = (d[-2]+d[-1])/2
    return smooth(x, c=c-1)


def smooth_err(d, err, c=5):
    if c == 0:
        return d
    e = 1/err
    x = np.zeros(len(d))
    x[0] = (d[0]*e[0]+d[1]*e[1])/(e[0]+e[1])
    x[1:-1] = (2*d[1:-1]*e[1:-1]+d[2:]*e[2:]+d[:-2]*e[:-2])/(2*e[1:-1]+e[2:]+e[:-2])
    x[-1] = (d[-2]*e[-2]+d[-1]*e[-1])/(e[-2]+e[-1])
    return smooth_err(x, err, c=c-1)


def vari_index(d, err):
    return 1/(len(d)-1)*np.sum(np.square(d-np.mean(d))/err**2)


def download_light_curve(ra, dec):
    """
    Downloads the light curve of the target from the Catalina Sky Survey
    
    :param ra: RA coordinate in degree
    :type ra: float
    :param dec: Dec coordinate in degre
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


def dayly_average(d):
    """
    Takes the dayly average of the light curve to reduce the noise.
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
    :param dec: Dec coordinate in degre
    :type dec: float
    """
    lc = download_light_curve(ra, dec)
    
    pl.clf()
    sp = pl.subplot(211)
#    sp.errorbar(lc['MJD'],
#                lc['Mag'],
#                lc['Magerr'],
#                fmt='.k',
#                capsize=2)
    avg = dayly_average(lc)
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
    
    print('normal', vari_index(lc['Mag'].values, lc['Magerr'].values))
    print('avg', vari_index(avg['Mag'].values, avg['Magerr'].values))
    print('avg smoothed', vari_index(avg_smooth, avg['Magerr'].values))
    pl.show()


def test():
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    s = SkyCoord(155.032761*u.deg, 9.467287*u.deg)
    plot_light_curve(s.ra.degree, s.dec.degree)
