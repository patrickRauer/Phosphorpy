import pylab as pl


class SpectraPlot:
    _spectra = None

    def __init__(self, spectra):
        """
        SpectraPlot os the plotting environment of the :class:`Phosphorpy.data.sub.spectra.Spectra` class.

        :param spectra: The spectra
        :type spectra: Phosphorpy.data.sub.spectra.Spectra
        """
        self._spectra = spectra

    def spectra(self, path='', min_wavelength=None, max_wavelength=None):
        """
        Basic spectra plot

        :param path:
            Path, where to save the figure or an empty string, if no save is wanted.
            Default is an empty string.
        :type path: str
        :param min_wavelength:
            The minimal wavelength to plot or None if no limit is wanted.
            Default is None.
        :type min_wavelength: None, float
        :param max_wavelength:
            The maximal wavelength to plot or None, if no limit is wanted.
        :type max_wavelength: None, float
        :return:
        """
        pl.clf()
        sp = pl.subplot()

        wave = self._spectra.wavelength
        flux = self._spectra.flux

        if min_wavelength is not None or max_wavelength is not None:
            spec = self._spectra.cut(min_wavelength=min_wavelength,
                                     max_wavelength=max_wavelength,
                                     inplace=False)
            wave = spec.wavelength
            flux = spec.flux

        sp.step(wave, flux, '-k')

        sp.set_xlabel(f'wavelength [${self._spectra.wavelength_unit.to_string("latex")}$]')
        sp.set_ylabel(f'flux [${self._spectra.flux_unit.to_string("latex")}$]')

        if path == '':
            pl.show()
        else:
            pl.savefig(path)
