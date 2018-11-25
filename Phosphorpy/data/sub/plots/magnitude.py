import pylab as pl


def _hist(sp, x, bins, histtype, label=''):
    sp.hist(x, bins=bins, histtype=histtype, label=label)


class MagnitudePlot:
    _data = None

    def __init__(self, data):
        self._data = data

    def hist(self, cols=None, path=''):
        """
        Plots the histogram(s) of the different magnitude(s).

        :param cols: A list of magnitude names. Default is None which means that all magnitudes are taken.
        :type cols: list
        :param path: Path to the save place. Default is an empty string, which means that the plot will be shown, only.
        :type path: str
        :return:
        """
        all_columns = self._data.get_names()
        # if no columns are given, take all columns from the dataset
        if cols is None:
            cols = all_columns

        pl.clf()
        sp = pl.subplot()

        # go through all given magnitudes
        for c in cols:
            # check if the magnitude is in the data. If not, raise an error
            if c not in all_columns:
                raise KeyError('Magnitude is not in the dataset!')
            _hist(sp, self._data.data[c], bins='auto', histtype='step', label=c)

        # set the axis-labels
        sp.set_xlabel('magnitude')
        sp.set_ylabel('frequency')
        pl.legend(loc='best')

        # if a path is given, save the figure, else show the figure
        if path != '':
            pl.savefig(path)
        else:
            pl.show()
