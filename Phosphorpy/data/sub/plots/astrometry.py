import numpy as np
import pylab as pl


class AstrometryPlot:
    _astrometry = None

    def __init__(self, astrometry):
        self._astrometry = astrometry

    def proper_motion(self, path='', cos_correction=False):
        """
        Makes a plot of the proper motion

        :param path: Path to the location
        :type path: str
        :param cos_correction: True if the pmra should be correct with cos(deg), else False. Default is False.
        :type cos_correction: bool
        :return:
        """
        pms = self._astrometry.proper_motion(cos_correction)

        pl.clf()
        sp = pl.subplot()
        sp.errorbar(pms['pmra'], pms['pmdec'], xerr=pms['pmra_err'], yerr=pms['pmdec_err'], fmt='.', capsize=2)
        for i in range(1, self._astrometry.mask.get_mask_count()):
            m = self._astrometry.mask.get_mask(i).copy()
            m.index.name = pms.index.name
            dpms = pms[m]
            sp.errorbar(dpms['pmra'], dpms['pmdec'], xerr=dpms['pmra_err'], yerr=dpms['pmdec_err'], fmt='.', capsize=2)

        if cos_correction:
            sp.set_xlabel('$\\mu_\\alpha^*$ [mas/yr]')
        else:
            sp.set_xlabel('$\\mu_\\alpha$ [mas/yr]')
        sp.set_ylabel('$\\mu_\\delta$ [mas/yr]')

        if path != '':
            pl.savefig(path)
        else:
            pl.show()

    def _hist(self, x, xlabel, path=''):
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
        pl.clf()
        sp = pl.subplot()
        hist_range = [np.nanmin(x), np.nanmax(x)]
        sp.hist(x, bins='auto', range=hist_range, histtype='step', color='k')
        for i in range(1, self._astrometry.mask.get_mask_count()):
            sp.hist(x[self._astrometry.mask.get_mask(i)].values,
                    bins='auto', range=hist_range, histtype='step')
        sp.set_xlabel(xlabel)
        sp.set_ylabel('count')

        if path != '':
            pl.savefig(path)
        else:
            pl.show()

    def proper_motion_hist(self, path=''):
        """
        Plots a histogram of the total proper motion

        :param path: The path to the location
        :type path: str
        :return:
        """
        pms = self._astrometry.proper_motion(True)
        pm = np.hypot(pms['pmra'], pms['pmdec'])
        self._hist(pm, '$\\mu$ [mas/yr]', path)

    def parallax_hist(self, path=''):
        """
        Plots a histogram of the parallax

        :param path: The path to the location
        :type path: str
        :return:
        """
        parallax = self._astrometry.data['parallax']
        self._hist(parallax, 'p [mas]', path)

    def parallax(self, path=''):
        """
        Plot the parallax vs its error

        :param path: The path to the location
        :type path: str
        :return:
        """
        pl.clf()
        sp = pl.subplot()
        parallax = self._astrometry.data['parallax']
        parallax_error = self._astrometry.data['parallax_error']
        sp.scatter(parallax, parallax_error, marker='.', c='k')

        for i in range(1, self._astrometry.mask.get_mask_count()):
            m = self._astrometry.mask.get_mask(i)
            sp.scatter(parallax[m].values, parallax_error[m].values, marker='.')

        sp.set_xlabel('p [mas]')
        sp.set_ylabel('$\\sigma_p$ [mas]')

        if path != '':
            pl.savefig(path)
        else:
            pl.show()
