from .table import DataTable
from Phosphorpy.data.sub.plots.flux import FluxPlot
from pandas import DataFrame
import numpy as np


def get_column_names(length):
    # create columns names 'a', 'b', 'c', ...
    cols = []
    for i in range(97, 98 + length):
        cols.append(chr(i))
    return cols


class FluxTable(DataTable):
    """
    FluxTable is the interface to interact with the fluxes of the sources.
    """

    _fits = None
    _data = None
    _survey = None

    def __init__(self, data, survey_head=None, mask=None):
        """

        :param data:
        :param survey_head: Survey information
        :type survey_head: SearchEngine.data.sub.magnitudes.SurveyData
        """
        DataTable.__init__(self, mask=mask)
        self._data = data
        self._survey = survey_head
        self._plot = FluxPlot(self)

    def fit_polynomial(self, degree, error_weighted=True,
                       x_log=False, y_log=False,
                       norm=False):
        d = self._data[self._survey.all_magnitudes()].values
        err = self._data[self._survey.all_error_names()].values

        wavelengths = self._survey.get_all_wavelengths()

        if x_log:
            wavelengths = np.log10(wavelengths)

        if y_log:
            err = np.abs(1/d)*err
            d = np.log10(d)

        if norm:
            maxi = np.max(d, axis=1)
            maxi = np.transpose(np.tile(maxi, (d.shape[1], 1)))
            d = d/maxi
            err = err/maxi

        fits = []
        if error_weighted:
            for row, err_row in zip(d, err):
                fit = np.polyfit(wavelengths, row, degree, w=1/err_row)
                fits.append(fit)
        else:
            for row in d:
                fit = np.polyfit(wavelengths, row, degree)
                fits.append(fit)

        self._fits = DataFrame(data=np.array(fits), columns=get_column_names(degree))

    def fit_blackbody(self):
        pass

    def __str__(self):
        return str(self._data)

    @property
    def fit(self):
        return self._fits
