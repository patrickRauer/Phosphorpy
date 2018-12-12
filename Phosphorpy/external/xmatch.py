"""
This script wraps the astroquery.xmatch.XMatch class to create an interface which is more usable in this program.
"""
from astroquery.xmatch import XMatch
from astropy.table import Table
from astropy import units as u

from ..config.survey_data import SURVEY_DATA
import numpy as np
import pandas
import os


def __check_row_id__(data, colnames):
    """
    Checks if the input data has a column 'row_id'. If not it will add such a column with unique integer id's

    :param data: The input data
    :type data: astropy.table.Table, pandas.DataFrame
    :param colnames: List with the names of the columns of the input data
    :type colnames: list
    :return: The input data with an additional column 'row_id' if there wasn't such a column before.
    :rtype: astropy.table.Table, pandas.DataFrame
    """
    if 'row_id' not in colnames:
        data['row_id'] = np.linspace(1, len(data), len(data), dtype=np.int32)
    return data


def __write_temp_file__(data, ra_name, dec_name):
    """
    Writes a temporary csv-file with the input data and adds a column with id's if there is no such column

    :param data: The input data
    :type data: astropy.table.Table, pandas.DataFrame
    :param ra_name: The name of the ra column
    :type ra_name: str
    :param dec_name: The name of the Dec column
    :type dec_name: str
    :return:
    """
    coords = ['{}_input'.format(ra_name), '{}_input'.format(dec_name), 'row_id']

    # if the input data are a pandas.DataFrame
    if type(data) == pandas.DataFrame:
        if 'input' not in ra_name:
            data = data.rename({ra_name: '{}_input'.format(ra_name),
                                dec_name: '{}_input'.format(dec_name)},
                               axis='columns')
        data = __check_row_id__(data, data.columns)
        data[coords].to_csv('temp.csv')
    # if the input data are an astropy.table.Table
    elif type(data) == Table:
        if 'input' not in ra_name:
            data.rename_column(ra_name, '{}_input'.format(ra_name))
            data.rename_column(dec_name, '{}_input'.format(dec_name))
        data = __check_row_id__(data, data.colnames)
        data[coords].write('temp.csv', format='ascii.csv')
    # if the input data aren't an astropy.table.Table or a pandas.DataFrame
    # raise a TypeError
    else:
        raise TypeError('Unsupported data type!')


def __output_columns__(survey):
    """
    Reads the magnitude names and additional columns from the configs and format them in the vizier style.

    :param survey: The name of the survey
    :type survey: str
    :return: The required column names in the vizier style
    :rtype: list
    """
    cols = []
    # add the magnitude column names to the column list
    # use the vizier style (mag_name+mag)
    for c in SURVEY_DATA[survey]['magnitude']:
        cols.append('{}mag'.format(c))
        cols.append('e_{}mag'.format(c))

    # check if additional columns are set. If yes, add them to the column list
    if 'columns' in SURVEY_DATA[survey].keys():
        cols.extend(SURVEY_DATA[survey]['columns'])
    coord_cols = []
    if 'coordinate' in SURVEY_DATA[survey].keys():
        coord_cols.extend(SURVEY_DATA[survey]['coordinate'])
    else:
        coord_cols.extend(['RAJ2000', 'DEJ2000'])
    cols.extend(coord_cols)
    cols.append('row_id')
    return cols, coord_cols


def find_suffix(cols):
    for c in cols:
        if 'mag' in c:
            pre = c.split('mag')[0]
            return c.split(pre)[-1]
        elif 'Mag' in c:
            if 'Aper' in c:
                pre = c.split('AperMag')[0]
            else:
                pre = c.split('Mag')[0]
            return c.split(pre)[-1]


def xmatch(data, ra_name, dec_name, survey, max_distance=1.*u.arcsec, blank=False):
    """
    Interface to the astroquery.xmatch.XMatch module

    :param data: The input data
    :type data: astropy.table.Table, pandas.DataFrame
    :param ra_name: The name of the ra column
    :type ra_name: str
    :param dec_name: The name of the Dec column
    :type dec_name: str
    :param survey: The name of the survey
    :type survey: str
    :param max_distance: Maximal distance to the counterpart in the other catalog
    :type max_distance: ﻿astropy.units.quantity.Quantity
    :param blank: True if all columns should return, else False for survey specific columns.
    :type blank: bool
    :return: The results of the catalog query
    :rtype: pandas.DataFrame
    """
    # write the temp file
    __write_temp_file__(data, ra_name, dec_name)

    # use XMatch to download the data
    rs = XMatch.query(cat1=open('temp.csv'),
                      cat2='vizier:{}'.format(SURVEY_DATA[survey]['vizier']),
                      colRA1='{}_input'.format(ra_name),
                      colDec1='{}_input'.format(dec_name),
                      max_distance=max_distance)

    # remove the temp file
    os.remove('temp.csv')
    if blank:
        return rs

    # reduce the number of columns to the required ones
    output_cols, coord_cols = __output_columns__(survey)
    # rs = rs[output_cols]

    rs = rs.to_pandas()

    # handle the labeling of XMatch Gaia
    if survey == 'GAIA':
        conv = {'phot_g_mean_mag': "Gmag",
                'phot_bp_mean_mag': 'BPmag',
                'phot_rp_mean_mag': 'RPmag',
                'parallax': 'Plx',
                'parallax_error': 'e_Plx'}
        rs = rs.rename(conv, axis='columns')
        rs['e_Gmag'] = 2.5*rs['phot_g_mean_flux_error']/(np.log(10) *
                                                         rs['phot_g_mean_flux'])
        rs['e_BPmag'] = 2.5*rs['phot_bp_mean_flux_error']/(np.log(10) *
                                                           rs['phot_bp_mean_flux'])
        rs['e_RPmag'] = 2.5*rs['phot_rp_mean_flux_error']/(np.log(10) *
                                                           rs['phot_rp_mean_flux'])

    else:
        # check if the config-file contains special hints for the labeling
        if 'xmatch' in SURVEY_DATA[survey].keys():
            replacer = SURVEY_DATA[survey]['xmatch']

            # rename the column names to the standard vizier style
            conv = {}
            for c in rs.columns:

                if 'Err' in c:
                    conv[c] = 'e_{}mag'.format(c.replace('{}Err'.format(replacer), '').upper())
                elif 'err' in c:
                    conv[c] = 'e_{}mag'.format(c.replace('{}err'.format(replacer), '').upper())
                elif replacer in c:
                    conv[c] = '{}mag'.format(c.replace(replacer, '').upper())
            rs = rs.rename(conv, axis='columns')

            # special case of UKIDSS J-band
            conv = {}
            for c in rs.columns:

                if '_1' in c:
                    conv[c] = c.replace('_1', '')
            rs = rs.rename(conv, axis='columns')

    # check if output_cols contains the right coordinate names
    if output_cols[-2] not in rs.columns:
        if 'ra' in rs.columns:
            output_cols[-3] = 'ra'
            output_cols[-2] = 'dec'
        elif 'RAdeg' in rs.columns:
            output_cols[-3] = 'RAdeg'
            output_cols[-2] = 'DEdeg'

    rs = rs[output_cols]
    # group by the row id (unique identifier of the input sources)
    # and take mean value of all values without the NaN values
    rs = rs.groupby('row_id')
    rs = rs.aggregate(np.nanmean)

    return rs
