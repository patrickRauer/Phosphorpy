#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:07:21 2019

@author: Jean Patrick Rauer
"""
from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd

from Phosphorpy.core.structure import Table
from Phosphorpy.data.sub.plots.color import create_color_name


class Color(Table):
    
    _survey = None
    
    def __init__(self, data, name, mask=None):
        Table.__init__(self, data, name, mask=mask)

    def __str__(self):
        return 'Colors:\n{}\n'.format(len(self))

    def __repr__(self):
        return str(self)
        
    def __get_mask_data__(self, col, minimum, maximum, previous):
        col = create_color_name(col)
        d = self[col].values
        mask = (d < maximum) & (d > minimum)
        self._mask.add_mask(mask, 'Color cut (minimum={}, maximum={}'.format(minimum,
                                                                             maximum),
                            combine=previous)

    def get_columns(self, cols):
        use_cols = []
        col_names = self.columns
        for c in cols:
            if c in col_names:
                use_cols.append(c)
        return self[use_cols]
        
    def outlier_detection(self):
        """
        Makes an outlier detection based on a DBSCAN
        
        .. warning:
            
            In the case of very small amount of data points it can happen
            that all data points are excluded. 
            Use this outlier detection for larger amount of data (N > 1000).
        """
        db = DBSCAN(eps=2, min_samples=5)
        values = self.values
        for i in range(len(values[0])):
            values[:, i] /= np.std(values[:, i])
        db.fit(self.values)
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
                for c in self.column:
                    if col[0] in c and col[1] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            elif len(col) == 3:
                for c in self.column:
                    if col[0] in c and col[1] in c and col[2] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            else:
                raise ValueError('The list must have 2 or 3 elements. Not more and not less.')
        else:
            raise ValueError('Sorry, I don\'t how I should handle col in the format {}. Please Use a string or '
                             'tuple/list to specify the right color')
