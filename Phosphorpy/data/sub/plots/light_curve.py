import pylab as pl


class LightCurvePlot:

    _light_curve = None

    def __init__(self, light_curve):
        self._light_curve = light_curve

    def plot_light_curve(self, light_curve_id, path=None):
        pl.clf()
        sp = pl.subplot()
        if type(light_curve_id) is int:
            lc = self._light_curve.light_curves[self._light_curve.light_curves['InputID'] == light_curve_id]
            sp.errorbar(lc['MJD'].values,
                        lc['Mag'].values,
                        lc['Magerr'].values,
                        fmt='.k', capsize=2)
        elif type(light_curve_id) is tuple or type(light_curve_id) is list:
            for lci in light_curve_id:
                lc = self._light_curve.light_curves[self._light_curve.light_curves['InputID'] == lci]
                sp.errorbar(lc['MJD'].values,
                            lc['Mag'].values,
                            lc['Magerr'].values,
                            fmt='.', capsize=2, label=str(lci))
            pl.legend(loc='best')

        sp.invert_yaxis()

        sp.set_xlabel('MJD [days]')
        sp.set_ylabel('mag')
        if path is not None:
            pl.savefig(path)
        else:
            pl.show()
