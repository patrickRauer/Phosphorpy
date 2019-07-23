from astropy.table import Table
from astropy import units as u
from astropy import constants as cons
from scipy.optimize import curve_fit
import numpy as np

from Phosphorpy.core.functions import gaus
from Phosphorpy.data.sub.plots.spectra import SpectraPlot
from Phosphorpy.data.sub.spectra.lines import SpectraLine


class Spectra:

    _wavelength = None
    _flux = None
    _header = None
    _lines = None

    _wavelength_unit = None
    _flux_unit = None

    _estimations = None

    _plot = None

    def __init__(self, wavelength=None, flux=None):
        """
        Spectra is a meta class for spectra handling in general

        :param wavelength:
        :param flux:
        """
        if wavelength is not None and flux is not None:
            if len(wavelength) != len(flux):
                raise ValueError('Wavelength and flux must have the same length.')

            self._wavelength = wavelength
            self._flux = flux
        elif wavelength is not flux:
            raise ValueError('Wavelength and flux must be set or non of them.')

        self._wavelength_unit = u.angstrom
        self._flux_unit = ''

        self._lines = (SpectraLine('Ha', 6562.8, 'H_\\alpha'),
                       SpectraLine('CaII-a', 8498, 'CaII_a'),
                       SpectraLine('CaII-b', 8542, 'CaII_b'),
                       SpectraLine('CaII-c', 8662, 'CaII_c'))
        self._plot = SpectraPlot(self)

    @property
    def plot(self):
        return self._plot

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def wavelength_unit(self):
        return self._wavelength_unit

    @property
    def flux(self):
        return self._flux

    @property
    def flux_unit(self):
        return self._flux_unit

    @property
    def header(self):
        return self._header

    @property
    def lines(self):
        return self._lines

    @lines.setter
    def lines(self, value):
        if type(value) != dict:
            raise AttributeError('Line information must be a dict with the keys as lines names and '
                                 'the values as the vacuum wavelength')
        line_data = []
        for v in value:
            d = value[v]
            type_d = type(d)
            if type_d == list or type_d == tuple:
                line_data.append(SpectraLine(v, *d))
            else:
                line_data.append(SpectraLine(v, d))
        self._lines = line_data

    def sub_range(self, min_wave, max_wave):
        """
        Returns a sub-range of the spectra

        :param min_wave: The minimal wavelength
        :type min_wave: float
        :param max_wave: The maximal wavelength
        :type max_wave: float
        :return: wavelength, flux
        """
        wavelength = self.wavelength
        flux = self.flux

        if min_wave < wavelength.min():
            raise ValueError(f'The lower wavelength limit is smaller than the minimal wavelength.\n'
                             f'{min_wave} < {wavelength.min()}')
        if max_wave > wavelength.max():
            raise ValueError(f'The upper wavelength limit is bigger than the maximal wavelength.\n'
                             f'{max_wave}>{wavelength.max()}')

        m = wavelength >= min_wave
        m &= wavelength <= max_wave

        wavelength = wavelength[m]
        flux = flux[m]
        return Spectra(wavelength=wavelength, flux=flux)

    def bin(self, bins, left=False):
        """
        Returns a binned version of the spectra.

        :param bins: The number of bins of the resulting spectra
        :type bins: int
        :param left: True if the cut should on the left side, else False. Default is False.
        :type left: bool
        :return: The binned spectra
        :rtype: Spectra
        """
        if bins <= 0:
            raise ValueError('At least one bin must be left.')

        bin_size = len(self.wavelength) // bins
        l_max = bins*bin_size

        if left:
            wavelength = self.wavelength[-l_max:]
            flux = self.wavelength[-l_max:]
        else:
            wavelength = self.wavelength[:l_max]
            flux = self.wavelength[:l_max]

        wavelength = wavelength.reshape((bins, bin_size))
        flux = flux.reshape((bins, bin_size))

        wavelength = np.sum(wavelength, axis=1)
        flux = np.sum(flux, axis=1)
        return Spectra(wavelength=wavelength, flux=flux)

    def _format_output(self, out, as_velocity):
        """
        Formats the output-table of the line property estimations

        :param out:
        :param as_velocity:
        :return:
        """
        out = out[['name', 'z', 'zerr', 'lb', 'lberr']]

        if as_velocity:
            for c in out.colnames:
                if c != 'name':
                    out[c] = (out[c]*cons.c.to(u.km/u.s)).round(2)
            out.rename_columns(['z', 'zerr'], ['RV', 'RVerr'])
            out['RV'] += self.header['HELIO_RV']
            out['RV'] = out['RV'].round(2)
        else:
            for c in out.colnames:
                if c != 'name':
                    out[c] = out[c].round(7)
        return out

    def _fit_line(self, line):
        """

        :param line: The spectral line information
        :type line: SpectraLine
        :return:
        """
        delta_lambda = 20
        m = (self.wavelength > line.lambda0-delta_lambda) & (self.wavelength < line.lambda0+delta_lambda)
        w = self.wavelength[m]
        fl = self.flux[m]
        fl -= np.mean(fl)
        popt, pcov = curve_fit(gaus, w, fl, p0=[-1, line.lambda0, 10, 0, 0])
        perr = np.sqrt(np.diag(pcov))

        return {'z': popt[1]/line.lambda0-1, 'zerr': perr[1]/line.lambda0,
                'lb': abs(popt[2]/line.lambda0), 'lberr': perr[2]/line.lambda0}

    def estimate_line_properties(self, as_velocity=False, redo=False):
        """
        Estimates the line shift and the line broadening with

        .. math:

            f(\lambda) = a*e^{-\frac{(\lambda-lambda_0)**2}{2*\sigma**2}}+b*\lambda+c

        If the output should be in km/s, the radial velocity got a heliocentric correction.
        The value for the heliocentric correction is taken from the header of the spectra.

        :param as_velocity:
            True, if the line-shift should be in km/s with a heliocentric correction else False.
            Default is False.
        :type as_velocity: bool
        :param redo:
            True, if the estimations should be done again, else False. Default is False.
            It does not have an effect, if no estimations are done before.
        :type redo: bool
        :return: Table with the results of the lines
        :rtype: astropy.table.Table
        """
        if not redo and self._estimations is not None:
            return self._estimations
        out = []
        for l in self.lines:
            line_data = {'name': l.name}
            line_data.update(self._fit_line(l))
            out.append(line_data)
        out = Table(rows=out)

        out = self._format_output(out, as_velocity)
        self._estimations = out.copy()
        return out
