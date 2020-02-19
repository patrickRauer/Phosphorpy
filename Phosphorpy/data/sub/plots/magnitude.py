import numpy as np
import pylab as pl
import warnings

def _hist(sp, x, bins, histtype, label=''):

    sp.hist(x, bins=bins, histtype=histtype, label=label, range=[np.nanmin(x),
                                                                 np.nanmax(x)])


class MagnitudePlot:
    _data = None

    def __init__(self, data):
        self._data = data

    def hist(self, cols=None, survey=None, path=''):
        """
        Plots the histogram(s) of the different magnitude(s).

        :param cols: A list of magnitude names. Default is None which means that all magnitudes are taken.
        :type cols: list, str
        :param survey: Name of the survey or None. If None all surveys with the cols names are used.
        :type survey: str
        :param path: Path to the save place. Default is an empty string, which means that the plot will be shown, only.
        :type path: str
        :return:
        """
        if survey is not None:
            d = self._data.get_survey_data(survey)
        elif len(self._data.data) == 1:
            d = self._data.data[0]
        else:
            d = self._data.get_columns(cols)
        all_columns = d.columns.values
        # if no columns are given, take all columns from the dataset
        if cols is None:
            cols = all_columns
        elif type(cols) == str:
            cols = [cols]

        # check every input column if it contains 'mag' in its name, if not add 'mag'
        # to the end of the column name to convert it to the default style
        for i, c in enumerate(cols):
            if 'mag' not in c:
                cols[i] = f'{c}mag'

        pl.clf()
        sp = pl.subplot()

        column_not_found = 0
        # go through all given magnitudes
        for c in cols:
            # check if the magnitude is in the data. If not, raise an error
            if c not in all_columns:
                warnings.warn(f'Magnitude ({c}) is not in the Dataset!')
                column_not_found += 1
                continue

            if self._data.mask.get_mask_count() == 0:
                _hist(sp, d[c],
                      bins='auto', histtype='step',
                      label=c.split('mag')[0])
            else:
                for i in range(self._data.mask.get_mask_count()):
                    _hist(sp, d[c][self._data.mask.get_mask(i).values],
                          bins='auto', histtype='step', label=c.split('mag')[0])

        if column_not_found == len(cols):
            raise ValueError('At least oe column must be in one survey available.')
        
        # set the axis-labels
        sp.set_xlabel('magnitude')
        sp.set_ylabel('counts')
        pl.legend(loc='best')

        # if a path is given, save the figure, else show the figure
        if path != '':
            pl.savefig(path)
        else:
            pl.show()
