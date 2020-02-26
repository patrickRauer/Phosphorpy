from matplotlib.cm import get_cmap
from collections.abc import Iterable
import numpy as np
import pylab as pl

SURVEY = ('CSS', 'PTF', 'ZTF')


def _plot_light_curve(sp, lc, label=None, color=None):
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
    for s in unique_surveys:
        m = markers[s-1]
        l = lc[lc['survey'] == s]
        sp.errorbar(l['mjd'].values,
                    l['mag'].values,
                    l['magerr'].values,
                    fmt=m, capsize=2, alpha=0.2, color=color)
        sp.scatter(l['mjd'].values, l['mag'].values, marker=m,
                   label=f'{label} {SURVEY[s-1]}', color=color)


class LightCurvePlot:

    _light_curve = None

    def __init__(self, light_curve):
        self._light_curve = light_curve

    def plot_light_curve(self, light_curve_id, min_mjd=None, max_mjd=None, path=None):
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

        pl.clf()
        sp = pl.subplot()
        if type(light_curve_id) is int:
            lc = self._light_curve.get_light_curve(light_curve_id).light_curves
            lc = lc[(lc['mjd'] >= min_mjd) & (lc['mjd'] <= max_mjd)]
            _plot_light_curve(sp, lc)

        elif isinstance(light_curve_id, Iterable):
            colors = get_cmap('Set1').colors
            for i, lci in enumerate(light_curve_id):
                lc = self._light_curve.get_light_curve(lci).light_curves
                lc = lc[(lc['mjd'] >= min_mjd) & (lc['mjd'] <= max_mjd)]
                _plot_light_curve(sp, lc, str(lci),
                                  colors[i % len(colors)])
            pl.legend(loc='best')

        sp.invert_yaxis()

        sp.set_xlabel('MJD [days]')
        sp.set_ylabel('mag')

        if path is not None:
            pl.savefig(path)
        else:
            pl.show()

    def light_curve(self, light_curve_id, min_mjd=None, max_mjd=None, path=None):
        self.plot_light_curve(light_curve_id, min_mjd=min_mjd, max_mjd=max_mjd, path=path)
