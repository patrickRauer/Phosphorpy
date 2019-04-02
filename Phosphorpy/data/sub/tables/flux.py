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
        names = [n for n in self.data.columns if 'e_' not in n]
        return self.data[names]

    def get_errors(self):
        names = [n for n in self.data.columns if 'e_' in n]
        return self.data[names]

    def get_flux(self, index):
        return self.get_fluxes().iloc[index].values

    def get_error(self, index):
        return self.get_errors().iloc[index].values

    def get_wavelengths(self):
        print(self.survey)
