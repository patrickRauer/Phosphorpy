from astropy import units as u
import numpy as np

from Phosphorpy.data.sub.plots.spectra import SpectraPlot


class Spectra:
    NORMALIZE_MEAN = np.mean
    NORMALIZE_MEDIAN = np.median
    NORMALIZE_MAX = np.max
    NORMALIZE_SUM = np.sum

    _wavelength = None
    _flux = None

    _wavelength_unit = None
    _flux_unit = None

    _lines = None

    _plot = None

    def __init__(self, wavelength=None, flux=None, wavelength_unit=None, flux_unit=None):
        self._wavelength = wavelength
        self._flux = flux

        self.wavelength_unit = wavelength_unit
        self.flux_unit = flux_unit

        self._plot = SpectraPlot(self)

    def normalize(self, kind, inplace=True):
        """
        Normalize the spectra with a certain function.

        .. math::

            flux^* = \frac{flux}{function(flux)}


        :param kind:
            The kind of function, which is used to normalize the spectra.
            Allowed function are the static function or 'mean', 'median', 'max' or 'sum' as ky words.
        :type kind: function, str
        :param inplace:
            True, if the flux should be overwritten with the normalized flux, else False for a new spectra
            object. Default is True.
        :type inplace: bool
        :return: None or if inplace is False, a new spectra object.
        :rtype: None, Spectra
        """
        if type(kind) == str:
            kind = kind.lower()
        fl = self.flux.copy()
        if kind == self.NORMALIZE_MAX or kind == 'max':
            fl /= np.max(fl)
        elif kind == self.NORMALIZE_MEAN or kind == 'mean':
            fl /= np.mean(fl)
        elif kind == self.NORMALIZE_MEDIAN or kind == 'median':
            fl /= np.median(fl)
        elif kind == self.NORMALIZE_SUM or kind == 'sum':
            fl /= np.sum(fl)
        else:
            raise ValueError('Kind of normalization is not allowed.')

        if type(inplace) == bool:
            if inplace:
                self._flux = fl
                self._flux_unit /= self._flux_unit
            else:
                return Spectra(wavelength=self.wavelength.copy(),
                               flux=fl)
        else:
            raise TypeError('Inplace must be a bool-type.')

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def wavelength_unit(self):
        return self._wavelength_unit

    @wavelength_unit.setter
    def wavelength_unit(self, unit):
        if unit is None:
            self._wavelength_unit = u.angstrom
        elif type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        elif self._wavelength_unit is None:
            self._wavelength_unit = unit
        else:
            try:
                factor = (self._wavelength_unit.to(unit))
                self._wavelength *= factor
                self._wavelength_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a length.')

    @property
    def flux(self):
        return self._flux

    @property
    def flux_unit(self):
        return self._flux_unit

    @flux_unit.setter
    def flux_unit(self, unit):
        if unit is None:
            self._flux_unit = u.erg
        elif type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        elif self._flux_unit is None:
            self._flux_unit = unit
        else:
            try:
                factor = (self._flux_unit.to(unit))
                self._flux *= factor
                self._flux_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a flux.')

    @property
    def plot(self):
        return self._plot
