#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 08:43:54 2019

@author: Jean Patrick Rauer
"""
import numpy as np
import pandas as pd

from Phosphorpy.core.functions import power_2_10, subtract
from Phosphorpy.core.structure import Table
from .color import Color
from .flux import Flux


class Magnitude(Table):
    _mag_cols = []
    _err_cols = []
    _survey = None
    
    def __init__(self, data, names, survey='', mask=None):
        Table.__init__(self, data, survey, mask)
        self._set_cols(names)
        self.__remove_sources_without_data()

    def __remove_sources_without_data(self):
        m = np.array([False]*len(self.data))
        for c in self.data.columns:
            d = self.data[c].values
            m |= d == d
        self.data = self.data[m]

    def __str__(self):
        return f'Magnitude of {self.survey_name} with {len(self)} entries\n'

    def set_survey_data(self, survey_data):
        self._survey = survey_data
        
    def _set_cols(self, cols):
        mag = []
        err = []
        for c in cols:
            if 'e_' in c:
                err.append(c)
            else:
                mag.append(c)
        self._mag_cols = mag
        self._err_cols = err

    def apply_extinction_correction(self, correction):
        """
        Applies the extinction correction to the data

        :param correction: A table with the correction values
        :type correction: astropy.table.Table, pandas.DataFrame
        :return:
        """
        rows = np.arange(len(self))
        for c in self._mag_cols:
            self[c] -= correction[c][rows]

    def get_flux(self):
        """
        Returns the flux of the columns
        :return:
        """
        # transform the magnitude in the flux scale
        flux = Flux(self[self._mag_cols].apply(power_2_10).copy(),
                    mask=self.mask, survey=self.survey_name)

        # multiply the flux zero point of all bands of the different surveys
        for m in self._mag_cols:
            fo = self._survey.flux_zero(self.survey_name.lower(), m)
            flux[m] *= fo

        # compute the errors of the fluxes
        errs = self[self._err_cols]
        errs = pd.DataFrame(data=flux.values*errs.values*np.log(10.),
                            columns=errs.columns)

        # return a merged table out of the fluxes and the flux errors
        flux.merge(errs, left_index=True, right_index=True)
        return flux
        
    def get_colors(self, cols=None):
        """
        Computes all possible colors. The first magnitude will be the bluer and the second one will be the redder one,
        except magnitude columns are given, then this order is not fixed.

        .. math:

            C_{i, j} = mag_i - mag_j \\
            \\lambda(band_i) < \\lambda(band_j)

        :param cols:
        :return:
        """
        if cols is None:
            cols = self._mag_cols
        else:
            try:
                cols[0][0][0]
            except TypeError:
                cols = [cols]

        colors = Color(self.data[[]], self.survey_name, mask=self.mask)
        for i, c1 in enumerate(cols):
            for j, c2 in enumerate(cols):
                
                # compute only such colors where the first argument is in a bluer band
                # (assumption is that the columns are sorted in such a way)
                if i < j:
                    cc_name = '{} - {}'.format(c1, c2)
                    if len(self) < 1e4:
                        colors[cc_name] = self[c1].values - self[c2].values
                    else:
                        colors[cc_name] = subtract(self[c1].values, self[c2].values)

        return colors

    def get_columns(self, cols):
        use_columns = []
        for c in self.columns.values:
            if c in cols or c.split('mag')[0] in cols:
                use_columns.append(c)

        return self[use_columns]
    
    def has_full_photometry(self, previous=True):
        """
        Creates a mask of sources, which have photometry values in every band of the survey.

        :param previous: True if the previous mask should be used True, else False. Default is True.
        :type previous: bool
        :return:
        """
        mask = None

        for c in self._mag_cols:

            k = self[c] > -99
            if mask is None:
                mask = k
            else:
                mask = mask & k
        self.mask.add_mask(mask, 'Mask sources with full photometry in {}'.format(self._survey), combine=previous)

    def set_limit(self, band, minimum=99, maximum=-99, previous=True):
        """
        Sets a magnitude limit to the magnitude columns with the name band.

        .. code-block:: python

            # select all sources with a G magnitude between 16 and 18
            # (notice that minimum and maximum are magnitudes, therefore the minimum value is larger than the maximum)
            ds.magnitudes.set_limit('G', minimum=18, maximum=16)


        :param band: The name of the band/magnitude column
        :type band: str
        :param minimum: The minimal magnitude value. Default is 99, which means that no cut will be done.
        :type minimum: float
        :param maximum: The maximal magnitude value. Default is -99, which means that no cut will be done.
        :type maximum: float
        :param previous: True if the previous mask should be used True, else False. Default is True.
        :type previous: bool
        :return:
        """
        if band not in self._mag_cols:
            band_new = band + 'mag'
            if band_new not in self._mag_cols:
                raise ValueError('Magnitude column {} or {} not found.'.format(band, band_new))
            band = band_new
        cut_string = ''
        if minimum != 99:
            cut_string = 'minimum={}'.format(minimum)
        if maximum != -99:
            max_str = 'maximum={}'.format(maximum)
            if cut_string == '':
                cut_string = max_str
            else:
                cut_string = '{} and {}'.format(cut_string, max_str)

        mags = self[band]
        if maximum != -99 and minimum != 99:
            mask = (mags <= minimum) & (mags >= maximum)
        elif maximum != -99:
            mask = mags >= maximum
        elif minimum != 99:
            mask = mags <= 99
        else:
            raise AttributeError('At least one of maximum or minimum must be different to the default value.')
        self.mask.add_mask(mask, 'Apply a magnitude cut on the band {} with {}'.format(band, cut_string),
                           combine=previous)

    @property
    def mag_names(self):
        return self._mag_cols

    @property
    def err_names(self):
        return self._err_cols

    @property
    def survey(self):
        return self._survey
