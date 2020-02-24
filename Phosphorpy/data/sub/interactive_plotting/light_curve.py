from matplotlib.cm import get_cmap
from collections.abc import Iterable
import numpy as np
import holoviews as hv

SURVEY = ('CSS', 'PTF', 'ZTF')


def _plot_light_curve(lc, label=None, color=None, hover=True):
    """
    Plots the light curve of a specific target with different markers for every survey

    :param sp: The plotting environment
    :param lc: The light curve data
    :type lc: DataFrame
    :param label: The name of the target
    :type label: str
    :param color: The color of the light curve in the figure
    :return:
    """
    if color is None:
        color = 'k'

    if label is None:
        label = ''

    unique_surveys = np.unique(lc['survey'])
    markers = ['.', 'd', '+', 'x']
    graph = None
    for s in unique_surveys:
        m = markers[s-1]
        l = lc[lc['survey'] == s]
        error = hv.ErrorBars(l, 'mjd', ['mag', 'magerr'])
        g = hv.Scatter(l, 'mjd', 'mag').opts(size=5)

        if hover:
            g = g.opts(tools=['hover'])
        g = error*g

        if graph is None:
            graph = g
        else:
            graph *= g
    return graph


class LightCurvePlot:
    _opts = None
    _light_curve = None

    def __init__(self, light_curve):
        self._light_curve = light_curve
        self._opts = dict(
            xlabel='MJD [days]',
            ylabel='mag',
            tools=['hover']
        )

    def plot_light_curve(self, light_curve_id, min_mjd=None, max_mjd=None, path=None,
                         **hv_kwargs):
        """
        Plots the light curve of the CSS

        :param light_curve_id:
            The ID of the light curve source or a tuple of ID's, if multiple light curves should be plotted
        :type light_curve_id: int, tuple, list
        :param path: Path to the storage place. Default is None, which means that the plot will be shown only.
        :type path: str
        :return:
        """

        if min_mjd is None:
            min_mjd = self._light_curve.light_curves['mjd'].min()

        if max_mjd is None:
            max_mjd = self._light_curve.light_curves['mjd'].max()

        opts = self._opts.copy()
        if hv_kwargs is not None:
            opts.update(hv_kwargs)

        if type(light_curve_id) is int:
            lc = self._light_curve.light_curves[self._light_curve.light_curves['row_id'] == light_curve_id]
            lc = lc[(lc['mjd'] >= min_mjd) & (lc['mjd'] <= max_mjd)]
            graph = _plot_light_curve(lc)

        elif isinstance(light_curve_id, Iterable):
            colors = get_cmap('Set1').colors
            graph = None
            for i, lci in enumerate(light_curve_id):
                lc = self._light_curve.light_curves[self._light_curve.light_curves['row_id'] == lci]
                lc = lc[(lc['mjd'] >= min_mjd) & (lc['mjd'] <= max_mjd)]
                g = _plot_light_curve(lc, str(lci),
                                      colors[i % len(colors)])
                if graph is None:
                    graph = g
                else:
                    graph *= g
        else:
            raise ValueError('\'licht_curve_id\' must be an integer or an iterable object of integers.')

        graph = graph.opts(
            **opts
        )
        return graph

    def light_curve(self, light_curve_id, min_mjd=None, max_mjd=None, path=None, **hv_kwargs):
        self.plot_light_curve(light_curve_id, min_mjd=min_mjd, max_mjd=max_mjd, path=path, **hv_kwargs)
