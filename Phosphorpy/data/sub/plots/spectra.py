from astropy.units.core import Unit
from astropy.units import Quantity
from astropy import units as u
import pylab as pl
import numbers

from Phosphorpy.data.sub.spectra.lines import SpectraLine


def _adjust_units(spectral_line, wavelength_unit):
    if isinstance(spectral_line, numbers.Number):
        return (spectral_line*u.AA).to(wavelength_unit).value
    elif isinstance(spectral_line, Quantity):
        return spectral_line.to(wavelength_unit).value
    elif isinstance(spectral_line, SpectraLine):
        return spectral_line.get_lambda0(wavelength_unit).value
    else:
        raise ValueError('Line information not usable.')


def _label(label_func, label, unit):
    if type(unit) == Quantity:
        unit = unit.unit

    if type(unit) == Unit:
        label_func(f'{label} [{unit.to_string("latex")}]')

    elif unit != '':
        label_func(f'{label} [${unit}$]')
    else:
        label_func(label)


class SpectraPlot:
    _spectra = None

    def __init__(self, spectra):
        """

        :param spectra:
        :type spectra: Phosphorpy.data.sub.spectra.spectra.Spectra
        """
        self._spectra = spectra

    def _labels(self, sp):
        _label(sp.set_xlabel, 'wavelength', self._spectra.wavelength_unit)
        _label(sp.set_ylabel, 'flux', self._spectra.flux_unit)

    def spectra(self):
        """
        Plot the complete spectra

        :return:
        """
        pl.clf()
        sp = pl.subplot()
        sp.step(self._spectra.wavelength, self._spectra.flux, '-k')
        self._labels(sp)

        sp.set_xlim(self._spectra.wavelength.min(), self._spectra.wavelength.max())
        pl.show()

    def line(self, spectral_line, delta_wavelength=30):
        """

        :param spectral_line: The line around which the plot should be created

        :type spectral_line: int, float, Quantity, SpectralLine
        :param delta_wavelength:
        :return:
        """

        line_wavelength = _adjust_units(spectral_line, self._spectra.wavelength_unit)

        wavelength, flux = self._spectra.sub_range(line_wavelength-delta_wavelength,
                                                   line_wavelength+delta_wavelength)

        pl.clf()
        sp = pl.subplot()

        sp.step(wavelength, flux, '-k', where='mid')

        sp.plot([line_wavelength, line_wavelength],
                [flux.min(), flux.max()], '--r')

        self._labels(sp)
        sp.set_xlim(wavelength.min(),
                    wavelength.max())
        pl.show()
