#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:07:21 2019

@author: Jean Patrick Rauer
"""
import numpy as np
import pandas as pd
from collections.abc import Iterable
from sklearn.cluster import DBSCAN

from Phosphorpy.core.structure import Table
from Phosphorpy.data.sub.plots.color import create_color_name


def _add_mag(c):
    mag_count = c.count('mag')
    if mag_count == 0:
        out = f'{c[0]}mag - {c[-1]}mag'
    elif mag_count == 2 or mag_count == 1:
        return _add_mag(c.replace('mag', ''))
    else:
        raise ValueError('Can not understand your column names.')
    return out


class Color(Table):
    
    _survey = None
    
    def __init__(self, data, name, mask=None):
        Table.__init__(self, data, name, mask=mask)

    def __str__(self):
        return f'Colors:\n{len(self)}\n'

    def __repr__(self):
        return str(self)
        
    def __get_mask_data__(self, col, minimum, maximum, previous):
        if col not in self.columns.values:
            col = create_color_name(col)
        d = self[col]
        mask = (d < maximum) & (d > minimum)
        self.mask.add_mask(mask,
                           f'Color cut (minimum={minimum}, maximum={maximum}',
                           combine=previous)

    def get_columns(self, cols):
        """
        Returns all columns, which are or include a part of the given names.
        If no columns are found, then a ValueError raises.

        :param cols: The name (parts) of the required columns
        :type cols: Union, str
        :return: DataFrame with the required columns
        :rtype: pd.DataFrame
        """
        if type(cols) == str:
            cols = [cols]
        elif not isinstance(cols, Iterable):
            raise AttributeError('Input must be a string, list, tuple or a set.')

        use_cols = []
        col_names = self.columns
        for c in cols:
            if type(c) != str:
                raise ValueError('Elements must be strings.')
            c = _add_mag(c)
            if c in col_names:
                use_cols.append(c)
            else:
                for column in col_names:
                    if c in column:
                        use_cols.append(column)

        if len(use_cols) == 0:
            raise ValueError('No columns with the given name (parts) found.')

        use_cols = list(np.unique(use_cols))
        return self[use_cols]
        
    def outlier_detection(self):
        """
        Makes an outlier detection based on a DBSCAN.

        NaN values in the colors are replaced with the mean of the color.
        
        .. warning:
            
            In the case of very small amount of data points it can happen
            that all data points are excluded. 
            Use this outlier detection for larger amount of data (N > 1000).
        """
        db = DBSCAN(eps=2, min_samples=5)
        values = self.values.copy()
        for i in range(len(values[0])):
            v = values[:, i]
            v /= np.nanstd(v)
            k = v != v
            v[k] = np.nanmean(v)
            values[:, i] = v
        db.fit(values)
        mask_values = db.labels_ != -1
        mask_values = pd.Series(mask_values, self.index.values)
        self.mask.add_mask(mask_values, 'DBSCAN outlier detection')

    def set_limit(self, col, minimum=-99, maximum=99, previous=True):
        """
        Sets a constrain to the colors and create a new mask of it

        :param col: The columns of the constrain.
        :type col: str, list, tuple
        :param minimum: The minimal value
        :type minimum: float
        :param maximum: The maximal value
        :type maximum: float
        :param previous: True if the last mask must be True too, else False to create a complete new mask.
        :type previous: bool
        :return:
        """
        if type(col) == str:
            self.__get_mask_data__(col, minimum, maximum, previous)
        elif type(col) == list or type(col) == tuple:
            if len(col) == 2:
                for c in self.columns:
                    if col[0] in c and col[1] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            elif len(col) == 3:
                for c in self.columns:
                    if col[0] in c and col[1] in c and col[2] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            else:
                raise ValueError('The list must have 2 or 3 elements. Not more and not less.')
        else:
            raise ValueError('Sorry, I don\'t how I should handle col in the format {}. Please Use a string or '
                             'tuple/list to specify the right color')
