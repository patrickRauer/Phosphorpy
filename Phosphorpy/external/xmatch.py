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
    coords = [ra_name, dec_name, 'row_id']

    # if the input data are a pandas.DataFrame
    if type(data) == pandas.DataFrame:
        data = __check_row_id__(data, data.columns)
        data[coords].to_csv('temp.csv')
    # if the input data are an astropy.table.Table
    elif type(data) == Table:
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


def xmatch(data, ra_name, dec_name, survey, max_distance=2.*u.arcsec):
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
    :return: The results of the catalog query
    :rtype: pandas.DataFrame
    """
    # write the temp file
    __write_temp_file__(data, ra_name, dec_name)

    # use XMatch to download the data
    rs = XMatch.query(cat1=open('temp.csv'),
                      cat2='vizier:{}'.format(SURVEY_DATA[survey]['vizier']),
                      colRA1=ra_name,
                      colDec1=dec_name,
                      max_distance=max_distance)
    # remove the temp file
    os.remove('temp.csv')

    # reduce the number of columns to the required ones
    output_cols, coord_cols = __output_columns__(survey)
    rs = rs[output_cols]

    rs = rs.to_pandas()
    return rs
