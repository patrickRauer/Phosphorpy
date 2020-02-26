import numpy as np

from Phosphorpy.data.sub.interactive_plotting.interactive_plotting import HVPlot

try:
    import holoviews as hv
except ImportError:
    hv = None


class AstrometryPlot(HVPlot):
    _astrometry = None

    def __init__(self, astrometry):
        self._astrometry = astrometry

    def proper_motion(self, path='', cos_correction=False, **hv_kwargs):
        """
        Makes a plot of the proper motion

        :param path: Path to the location
        :type path: str
        :param cos_correction: True if the pmra should be correct with cos(deg), else False. Default is False.
        :type cos_correction: bool
        :return:
        """
        pms = self._astrometry.proper_motion(cos_correction)

        graph = hv.ErrorBars(pms, 'pmra', ['pmdec', 'pmra_err', 'pmdec_err'])

        for i in range(1, self._astrometry.mask.get_mask_count()):
            m = self._astrometry.mask.get_mask(i).copy()
            m.index.name = pms.index.name
            dpms = pms[m]
            graph *= self._hover(hv.Scatter(dpms, 'pmra', 'pmdec'))

        if cos_correction:
            xlabel = '$\\mu_\\alpha^*$ [mas/yr]'
        else:
            xlabel = '$\\mu_\\alpha$ [mas/yr]'
        ylabel = '$\\mu_\\delta$ [mas/yr]'

        graph = graph.opts(
            xlabel=xlabel,
            ylabel=ylabel,
            **hv_kwargs
        )
        return graph

    def _hist(self, x, xlabel, path='', **hv_kwargs):
        """
        Makes a histogram plot of the input data

        :param x: The input data
        :type x: Union, numpy.ndarray
        :param xlabel: The x-label
        :type xlabel: str
        :param path: Path to the location
        :type path: str
        :return:
        """
        hist_range = [np.nanmin(x), np.nanmax(x)]
        hist, edge = np.histogram(x, bins='auto', range=hist_range)

        graph = hv.Histogram((edge, hist)).opts(color='k')
        graph = self._hover(graph)
        for i in range(1, self._astrometry.mask.get_mask_count()):
            hist, edge = np.histogram(x[self._astrometry.mask.get_mask(i).values],
                                      bins='auto', range=hist_range)
            graph *= self._hover(hv.Histogram((edge, hist)))

        graph = graph.opts(xlabel=xlabel, ylabel='count', **hv_kwargs)

        return graph

    def proper_motion_hist(self, path='', **hv_kwargs):
        """
        Plots a histogram of the total proper motion

        :param path: The path to the location
        :type path: str
        :return:
        """
        pms = self._astrometry.proper_motion(True)
        pm = np.hypot(pms['pmra'], pms['pmdec'])
        self._hist(pm, '$\\mu$ [mas/yr]', path, **hv_kwargs)

    def parallax_hist(self, path='', **hv_kwargs):
        """
        Plots a histogram of the parallax

        :param path: The path to the location
        :type path: str
        :return:
        """
        parallax = self._astrometry.data['parallax']
        self._hist(parallax, 'p [mas]', path, **hv_kwargs)

    def parallax(self, path='', **hv_kwargs):
        """
        Plot the parallax vs its error

        :param path: The path to the location
        :type path: str
        :return:
        """

        graph = hv.Scatter(self._astrometry.data, 'parallax', 'parallax_error')

        for i in range(1, self._astrometry.mask.get_mask_count()):
            m = self._astrometry.mask.get_mask(i)
            graph *= hv.Scatter(self._astrometry.data[m], 'parallax', 'parallax_error')

        graph = graph.opts(
            xlabel='p [mas]',
            ylabel='$\\sigma_p$ [mas]',
            **hv_kwargs
        )

        return graph

    @staticmethod
    def holoviews():
        if hv is not None:
            return True
        else:
            return False
