import numpy as np
import pandas as pd
from astropy.table import Table

from Phosphorpy.data.sub.plots.light_curve import LightCurvePlot
from Phosphorpy.external.css import download_light_curves


class LightCurves:
    _stat_columns = None
    _stat_operations = None

    _light_curves = None
    _average = None

    _plot = None

    def __init__(self, coordinates=None, light_curves=None):
        if coordinates is not None:
            self._light_curves = download_light_curves(coordinates['ra'],
                                                       coordinates['dec'])
        elif light_curves is not None:
            self._light_curves = light_curves
        else:
            raise ValueError('coordinates or light curves must be given.')

        self._stat_columns = ['InputID', 'Mag', 'Magerr', 'RA', 'Decl', 'MJD']
        self._stat_operations = [np.mean, np.median, np.std, np.min, np.max, 'count']

        self._plot = LightCurvePlot(self)

    def __str__(self):
        out = ''.join(['Number of light curves: {}\n',
                       'with {} entries.'])
        out = out.format(len(np.unique(self._light_curves['InputID'])), len(self._light_curves))
        return out

    def stats(self):
        """
        Computes the statistics of the light curves

        :return: The statistics of the light curves
        :rtype: pandas.DataFrame
        """
        print(self.light_curves.columns)
        return self.light_curves[self._stat_columns].groupby('InputID').aggregate(self._stat_operations)

    def average(self, dt_max=1, overwrite=False):
        """
        Averages the light curve

        :param dt_max: The maximal time difference between two data points
        :type dt_max: float
        :param overwrite: True, if the previous results should be overwritten, else False. Default is False.
        :type overwrite: bool
        :return: The averaged light curves
        :rtype: LightCurves
        """
        if dt_max < 0:
            raise ValueError('\'dt_max must be larger than 0')

        if self._average is not None and overwrite:
            return self._average

        out = []
        for lc_id in np.unique(self._light_curves['InputID']):
            lc = self._light_curves[self._light_curves['InputID'] == lc_id]
            dt = lc['MJD'][1:].values - lc['MJD'][:-1].values
            p = np.where(dt > dt_max)[0]+1
            start = 0
            mags = []
            errs = []
            mjds = []
            ra = []
            dec = []
            for k in p:
                l = lc[start: k]
                err_sq = 1/l['Magerr'].values**2
                err_sq_sum = 1./np.sum(err_sq)
                m = np.sum(l['Mag'].values*err_sq)*err_sq_sum
                e = np.sum(err_sq*l['Magerr'].values)*err_sq_sum
                mags.append(m)
                errs.append(e)
                mjds.append(np.sum(l['MJD'].values*err_sq)*err_sq_sum)
                ra.append(np.sum(l['RA'].values*err_sq)*err_sq_sum)
                dec.append(np.sum(l['Decl'].values*err_sq)*err_sq_sum)

            out.append(pd.DataFrame({'Mag': mags, 'Magerr': errs, 'MJD': mjds,
                                     'InputID': len(errs)*[lc['InputID'].values[0]],
                                     'RA': ra, 'Decl': dec}))
        self._average = LightCurves(light_curves=pd.concat(out))
        return self._average

    def get_light_curve(self, index):
        """
        Returns the data of the light curve with the given index

        :param index: The index of the light curve
        :rtype index: int
        :return:
        """
        return LightCurves(light_curves=self._light_curves[self._light_curves['InputID'] == index])

    @property
    def light_curves(self):
        return self._light_curves

    @property
    def stat_columns(self):
        return self._stat_columns

    @stat_columns.setter
    def stat_columns(self, value):
        if type(value) == str:
            if value in self.light_curves.columns:
                self._stat_columns = value
            else:
                raise ValueError(f'{value} is not a column.')
        elif type(value) == (list or tuple):
            for v in value:
                if v not in self.light_curves.columns:
                    raise ValueError(f'{v} is not a columns.')
        else:
            raise ValueError('Str or tuple/list are allowed.')

    @property
    def plot(self):
        return self._plot

    def write(self, path, format='fits'):
        """
        Writes the light curve data to a file

        :param path: The path to the file
        :type path: str
        :param format: The format of the file ('fits' or 'csv')
        :type format: str
        :return:
        """
        if format == 'fits':
            Table.from_pandas(self.light_curves).write(path, overwrite=True)
        elif format == 'csv':
            self.light_curves.to_csv(path)
        else:
            raise ValueError(f'{format} is not supported.')

    @staticmethod
    def read(path, format='fits'):
        """
        Reads the light curve data from a fits or csv file

        :param path: The path to the file
        :type path: str
        :param format: The format of the file (fits or csv)
        :type format: str
        :return: The light curve data from the file
        :rtype: LightCurves
        """
        if format == 'fits':
            return LightCurves(light_curves=Table.read(path).to_pandas())
        elif format == 'csv':
            return LightCurves(light_curves=pd.read_csv(path, index_col=0))
        else:
            raise ValueError(f'{format} is not supported.')
