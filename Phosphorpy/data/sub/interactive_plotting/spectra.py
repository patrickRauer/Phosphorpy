from collections.abc import Iterable
import numpy as np
import warnings

from Phosphorpy.data.sub.interactive_plotting.interactive_plotting import HVPlot

try:
    import holoviews as hv
except ImportError:
    hv = None


class SpectraPlot(HVPlot):
    _opts = None
    _spectra = None

    def __init__(self, spectra):
        """
        SpectraPlot os the plotting environment of the :class:`Phosphorpy.data.sub.spectra.Spectra` class.

        :param spectra: The spectra
        :type spectra: Phosphorpy.data.sub.spectra.Spectra
        """
        self._spectra = spectra
        self._opts = dict(
            xlabel='wavelength',
            ylabel='flux'
        )

    def spectra(self, path='', min_wavelength=None, max_wavelength=None, normalize=False,
                **hv_kwargs):
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
        :param normalize:
            True if the spectra should be normalized to its median value, else False.
            Default is False.
        :type normalize: bool
        :return:
        """

        if path != '':
            warnings.warn('Variable \'path\' will be ignored in the interactive plotting environment. Use the'
                          ' save buttons instead.')

        if min_wavelength is not None or max_wavelength is not None:
            spec = self._spectra.cut(min_wavelength=min_wavelength,
                                     max_wavelength=max_wavelength,
                                     inplace=False)
            wave = spec.wavelength
            flux = spec.flux
        else:
            wave = self._spectra.wavelength
            flux = self._spectra.flux

        if normalize:
            flux /= np.nanmedian(flux)

        m = np.abs(flux) < 1e8

        wave_c = wave[m].copy()
        flux_c = flux[m].copy()

        graph = hv.Curve((
            wave_c, flux_c
        ))

        graph = self._hover(graph)

        if len(self._spectra.fit) > 0:
            for fit in self._spectra.fit:
                graph *= self._hover(hv.Curve((wave, fit(wave))))

        opts = self._opts.copy()
        if hv_kwargs is not None:
            opts.update(hv_kwargs)
        graph = graph.opts(
            **opts
        )
        return graph


class SpectraListPlot(HVPlot):
    _opts = None
    _spectra_list = None

    def __init__(self, spectra_list):
        """
        SpectraPlot os the plotting environment of the :class:`Phosphorpy.data.sub.spectra.Spectra` class.

        :param spectra_list: The spectra
        :type spectra_list: Phosphorpy.data.sub.spectra.SpectraList
        """
        self._spectra_list = spectra_list
        self._opts = dict(
            xlabel='wavelength',
            ylabel='flux'
        )

    def spectra(self, index, path='', min_wavelength=None, max_wavelength=None, normalize=False,
                **hv_kwargs):
        if type(index) == int:
            self._spectra_list.get_by_id(index)[0][0].plot.spectra(path, min_wavelength, max_wavelength,
                                                                   **hv_kwargs)
        elif isinstance(index, Iterable):

            opts = self._opts.copy()
            if hv_kwargs is not None:
                opts.update(hv_kwargs)

            graph = None
            for i in index:
                specs = self._spectra_list.get_by_id(i)
                for spec_id in range(len(specs)):
                    spec = specs[spec_id]
                    if type(spec) == tuple:
                        spec = spec[0]
                    spec = spec.cut(min_wavelength=min_wavelength,
                                    max_wavelength=max_wavelength,
                                    inplace=False)
                    wave = spec.wavelength
                    flux = spec.flux

                    m = np.abs(flux) < 1e8
                    if normalize:
                        flux /= np.nanmedian(flux[m])

                    if graph is None:
                        graph = hv.Curve(
                            (
                                wave[m],
                                flux[m]
                            ),
                            label=f'{i}')
                        graph = self._hover(graph)
                    else:
                        g = hv.Curve(
                            (
                                wave[m],
                                flux[m]
                            ),
                            label=f'{i}')
                        graph *= self._hover(g)

            graph = graph.opts(
                **opts
            )
            return graph
        else:
            raise ValueError('Only integer or iterables like tuple or list with are allowed.')
