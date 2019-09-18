#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:03:54 2019

@author: Jean Patrick Rauer
"""

from Phosphorpy.core.structure import Table


class Flux(Table):
    
    _fits = None
    _survey = None
    _mask = None
    
    def __init__(self, data, survey=None, mask=None):
        Table.__init__(self, data, survey, mask)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self.data)

    def get_fluxes(self):
        """
        Returns all fluxes without their errors
        :return:
        """
        names = [n for n in self.data.columns if 'e_' not in n]
        return self.data[names]

    def get_errors(self):
        """
        Returns all errors without their fluxes
        :return:
        """
        names = [n for n in self.data.columns if 'e_' in n]
        return self.data[names]

    def get_index(self):
        """
        Returns the indices of the fluxes

        :return: The indices
        """
        return self.data[[]]

    def get_flux(self, index):
        """
        Returns the fluxes of one source without their errors

        :param index: The index of the source
        :type index: int
        :return:
        """
        return self.get_fluxes().iloc[index]

    def get_error(self, index):
        """
        Returns the errors of a source without the flux

        :param index: The index of the source
        :type index: int
        :return:
        """
        return self.get_errors().iloc[index]

    def get_wavelengths(self):
        pass

