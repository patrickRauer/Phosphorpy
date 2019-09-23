import numpy as np
import pylab as pl


def _plot_light_curve(sp, lc, label=None):
    if label is None:
        fmt = 'k'
    else:
        fmt = ''

    unique_surveys = np.unique(lc['survey'])
    markers = ['.', '+', 'd', 'x']
    for s, m in zip(unique_surveys, markers[:len(unique_surveys)]):
        l = lc[lc['survey'] == s]
        print(s, l)
        sp.errorbar(l['mjd'].values,
                    l['mag'].values,
                    l['magerr'].values,
                    fmt=m+fmt, capsize=2, alpha=0.2)
        sp.scatter(l['mjd'].values, l['mag'].values, marker=m,
                   label=label)


class LightCurvePlot:

    _light_curve = None

    def __init__(self, light_curve):
        self._light_curve = light_curve

    def plot_light_curve(self, light_curve_id, path=None):
        """
        Plots the light curve of the CSS

        :param light_curve_id:
            The ID of the light curve source or a tuple of ID's, if multiple light curves should be plotted
        :type light_curve_id: int, tuple, list
        :param path: Path to the storage place. Default is None, which means that the plot will be shown only.
        :type path: str
        :return:
        """
        pl.clf()
        sp = pl.subplot()
        print()
        if type(light_curve_id) is int:
            lc = self._light_curve.light_curves[self._light_curve.light_curves['row_id'] == light_curve_id]
            _plot_light_curve(sp, lc)

        elif type(light_curve_id) is tuple or type(light_curve_id) is list:
            for lci in light_curve_id:
                lc = self._light_curve.light_curves[self._light_curve.light_curves['row_id'] == lci]
                _plot_light_curve(sp, lc, str(lci))
            pl.legend(loc='best')

        sp.invert_yaxis()

        sp.set_xlabel('MJD [days]')
        sp.set_ylabel('mag')

        if path is not None:
            pl.savefig(path)
        else:
            pl.show()

    def light_curve(self, light_curve_id, path=None):
        self.plot_light_curve(light_curve_id, path=path)
