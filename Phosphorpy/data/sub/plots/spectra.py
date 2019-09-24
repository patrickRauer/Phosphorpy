import pylab as pl
import numpy as np


class SpectraPlot:
    _spectra = None

    def __init__(self, spectra):
        self._spectra = spectra

    def spectra(self, spec_id, path='', min_wavelength=None, max_wavelength=None):

        spec = self._spectra.get_spectra_by_ud(spec_id)
        pl.clf()
        sp = pl.subplot()

        wave = spec.wavelength
        flux = spec.flux

        if min_wavelength is not None:
            m = wave >= min_wavelength
            flux = flux[m]
            wave = wave[m]

        if max_wavelength is not None:
            m = wave <= max_wavelength
            flux = flux[m]
            wave = wave[m]

        sp.step(wave, flux, '-k')

        sp.set_xlabel('wavelength')
        sp.set_ylabel('flux')

        if path == '':
            pl.show()
        else:
            pl.savefig(path)
