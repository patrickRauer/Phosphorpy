import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import warnings

from Phosphorpy.external import zwicky
from Phosphorpy.data.sub.plots.light_curve import LightCurvePlot
from Phosphorpy.external.css import download_light_curves

try:
    from Phosphorpy.data.sub.interactive_plotting.light_curve import LightCurvePlot as LightCurvePlotHV
except ImportError:
    LightCurvePlotHV = None


def _average_light_curve(lc, dt_max):
    """
    Averages the given light curve

    :param lc: The input light curve
    :type lc: DataFrame
    :param dt_max: The box size of the average
    :type dt_max: float
    :return: The averaged light curve
    :rtype: DataFrame
    """
    dt = lc['mjd'][1:].values - lc['mjd'][:-1].values
    p = np.where(dt > dt_max)[0]+1
    start = 0
    mags = []
    errs = []
    mjds = []
    ra = []
    dec = []
    for k in p:
        l = lc[start: k]
        err_sq = 1/l['magerr'].values**2
        err_sq_sum = 1./np.sum(err_sq)
        m = np.sum(l['mag'].values*err_sq)*err_sq_sum
        e = np.sum(err_sq*l['magerr'].values)*err_sq_sum
        mags.append(m)
        errs.append(e)
        mjds.append(np.sum(l['mjd'].values*err_sq)*err_sq_sum)
        ra.append(np.sum(l['ra'].values*err_sq)*err_sq_sum)
        dec.append(np.sum(l['dec'].values*err_sq)*err_sq_sum)
    s = lc['survey'].values[0]
    avg = pd.DataFrame({'mag': mags, 'magerr': errs, 'mjd': mjds,
                        'row_id': len(errs)*[lc.index.values[0]],
                        'ra': ra, 'dec': dec, 'survey': np.linspace(s, s, len(errs))})
    return avg.set_index('row_id')


class LightCurves:
    _stat_columns = None
    _stat_operations = None
    _survey_names = None

    _light_curves = None
    _average = None

    _plot = None

    def __init__(self, coordinates=None, light_curves=None, surveys=None):
        if coordinates is not None:
            if surveys is None:
                surveys = ['css', 'ptf', 'ztf']
            else:
                if type(surveys) == str:
                    surveys = [surveys]
                surveys = [s.lower() for s in surveys]
            self._survey_names = np.array(surveys)
            out = []
            if 'css' in surveys:
                try:
                    css_lc = download_light_curves(coordinates['ra'],
                                                   coordinates['dec'])[['InputID', 'ID', 'Mag', 'Magerr',
                                                                        'RA', 'Decl', 'MJD']]
                    css_lc['survey'] = 1
                    css_lc = css_lc.rename({'InputID': 'row_id', 'ID': 'oid',
                                            'Mag': 'mag', 'Magerr': 'magerr', 'RA': 'ra',
                                            'Decl': 'dec', 'MJD': 'mjd'}, axis='columns')
                    out.append(css_lc)
                except:

                    warnings.warn('CSS light curve(s) are not available.')

            if type(coordinates) == pd.DataFrame:
                index = coordinates.index.values
            else:
                index = coordinates.data.index.values

            if 'ptf' in surveys:
                try:
                    ptf_lc = zwicky.download_ptf(coordinates['ra'], coordinates['dec'],
                                                 index=index)
                    ptf_lc['survey'] = 2
                    out.append(ptf_lc)
                except:
                    warnings.warn('PTF light curve(s) are not available.')

            if 'ztf' in surveys:
                try:
                    zwicky_lc = zwicky.download_light_curve(coordinates['ra'], coordinates['dec'],
                                                            index=index)
                    zwicky_lc['survey'] = 3
                    out.append(zwicky_lc)
                except:
                    warnings.warn('ZTF light curve(s) are not available.')

            if len(out) == 0:
                raise ValueError('No light curves found for the coordinate and survey combination.')

            self._light_curves = pd.concat(out)
            self._light_curves.set_index('row_id', inplace=True)
        elif light_curves is not None:
            self._light_curves = light_curves
        else:
            raise ValueError('coordinates or light curves must be given.')

        self._stat_columns = ['survey', 'mag', 'magerr', 'ra', 'dec', 'mjd']
        self._stat_operations = [np.mean, np.median, np.std, np.min, np.max, 'count']

        self._plot = LightCurvePlot(self)

        if LightCurvePlotHV.holoviews():
            self._hv_plot = LightCurvePlotHV(self)

    def __str__(self):
        out = ''.join(['Number of light curves: {}\n',
                       'with {} entries.'])
        out = out.format(len(np.unique(self._light_curves.index.values)), len(self._light_curves))
        return out

    def __repr__(self):
        return self.__str__()

    def survey_id2name(self, sid):
        return self._survey_names[sid-1]

    def stats(self):
        """
        Computes the statistics of the light curves

        :return: The statistics of the light curves
        :rtype: pandas.DataFrame
        """
        return self.light_curves[self._stat_columns].groupby(['row_id', 'survey']).aggregate(self._stat_operations)

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
        dt_min = np.min(self._light_curves['mjd'].values[1:]-self._light_curves['mjd'].values[:-1])
        if dt_max < 0:
            raise ValueError('dt_max must be larger than 0')
        elif dt_max < dt_min:
            raise ValueError(f'dt_max must be larger than the minimal difference between two data-points'
                             f'({dt_min})')

        if self._average is not None and overwrite:
            return self._average

        out = []
        for lc_id in np.unique(self._light_curves.index.values):
            lc = self._light_curves.loc[lc_id]
            for s in np.unique(lc['survey']):

                out.append(_average_light_curve(lc[lc['survey'] == s], dt_max))
        self._average = LightCurves(light_curves=pd.concat(out))
        return self._average

    def get_light_curve(self, index):
        """
        Returns the data of the light curve with the given index

        :param index: The index of the light curve
        :rtype index: int
        :return:
        """
        return LightCurves(light_curves=self._light_curves.loc[index])

    def align_light_curves(self, inplace=True):
        """
        Aligns the median of the light curves from the different surveys to median of the CSS light curve

        .. warning::

            Because of the simple median alignment, the aligned light curves can have systematic effects.
            For example, if a light curve as a general trend over the time of one survey, the correction will
            fail. This can happen, if light curves of QSO's are used or if a transient is present in one of the
            light curves without enough coverage before and after it.

        :param inplace: True, if the stored light curve should be overwritten, else False. Default is True.
        :type inplace: bool
        :return: The aligned light curves, if inplace is True. Otherwise None.
        :rtype: LightCurves, None
        """
        lc = self._light_curves.copy()
        stats = self.stats()
        mag0 = 0
        out = []
        for index in stats.index.values:
            mag_med = stats.loc[index][('mag', 'median')]
            lc_mask = (lc.index.values == index[0]) & (lc['survey'] == index[1])
            if index[1] == 1:
                mag0 = mag_med
                out.append(lc[lc_mask])
            else:
                lc_s = lc[lc_mask]
                lc_s.loc[:, 'mag'] -= mag_med-mag0
                out.append(lc_s)
        out = pd.concat(out)

        if inplace:
            self._light_curves = out
        else:
            return LightCurves(light_curves=out)

    def to_time_series(self, index):
        raise NotImplementedError()

    def to_astropy_table(self):
        """
        Returns the light curve data as an astropy Table

        :return: The light curve data points
        :rtype: Table
        """
        return Table.from_pandas(self.light_curves)

    def to_bin_table_hdu(self):
        return fits.BinTableHDU(self.to_astropy_table())

    def get_data(self, with_names=False):
        """
        Returns the raw data of the light curves.

        :param with_names:
            True if the survey ID's should be replaced by the names of the surveys, else False.
            Default is False.
        :type with_names: bool
        :return: The raw data of the light curves
        :rtype: DataFrame
        """
        data = self._light_curves.copy()
        if with_names:
            sid = data['survey'].values
            del data['survey']
            data['survey'] = self.survey_id2name(sid)
        return data

    @property
    def light_curves(self):
        return self._light_curves

    @property
    def data(self):
        return self.light_curves

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

    @property
    def hvplot(self):
        return self._hv_plot

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
