import sqlite3

import numpy as np
from astropy.io import fits
from astropy.table import Table

from Phosphorpy.core.structure import Mask


class DataTable:

    _mask = None
    _head = None
    _data = None

    _plot = None
    _hv_plot = None
    _q = [0.15, 0.25, 0.75, 0.85]

    _category = None

    def __init__(self, data=None, mask=None, category='table'):
        """
        Basic data table class
        """
        if data is not None:
            self._data = data
            if mask is None:
                mask = Mask(len(self._data))
            self._mask = mask
        elif mask is not None:
            self._mask = mask
        self._category = category

    def set_mask(self, mask):
        if type(mask) != Mask:
            raise ValueError('Mask must be a mask object.')
        self._mask = mask
        # if type(self.data) == list:
        #     for d in self.data:
        #         d.set_mask(mask)

    def stats(self):
        """
        Returns basic statistics (mean, median, std (unbiased), min, max) of the magnitudes

        :return: A DataFrame with the resulting statistics
        :rtype: pandas.core.frame.DataFrame
        """
        st = self.apply([np.mean, np.median, np.std, np.min, np.max])
        st = st.append(self.data.quantile(self.q))
        return st

    def apply(self, func):
        """
        Applies a function to the DataFrame. See `pandas apply method
        <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.apply.html>` for details.

        :param func: The to applied function
        :type: function, list
        :return:
            pandas Series or a DataFrame, in the case of multiple function, with the results of the applied function
        :rtype: pandas.core.series.Series, pandas.core.frame.DataFrame
        """
        return self.data.apply(func)

    def apply_on_ndarray(self, func):
        return func(self._data.values)

    def apply_on_dataframe(self, func):
        """
        Same as :meth:`apply`
        :param func:
        :return:
        """
        return self.apply(func)

    def remove_unmasked_data(self):
        """
        Removes all unmasked (mask[i] == False) from the data
        :return:
        """
        if self._data is not None:
            latest_mask = self._mask.get_latest_mask()
            latest_mask = latest_mask.align(self._data, fill_value=False)[0]

            self._data = self._data[latest_mask]

    def select_nan(self, column):
        """
        Select all rows with a NaN value
        :return:
        """
        if type(column) != str:
            for c in column:
                self.select_nan(c)
        d = self._data[column]
        self.mask.add_mask(d == d, f'Mask all NaN values in {column}')

    def min(self):
        """
        Returns the minimal values of all magnitudes (the magnitudes can come from different sources)

        :return: A pandas series with the column names as indices and the minimal values
        :rtype: pandas.core.series.Series
        """
        return self.apply(np.min)

    def max(self):
        """
        Returns the maximal values of all magnitudes (the magnitudes can come from different sources)

        :return: A pandas series with the column names as indices and the minimal values
        :rtype: pandas.core.series.Series
        """
        return self.apply(np.max)

    @property
    def shape(self):
        return len(self._data), len(self._data.columns)

    @property
    def head(self):
        return self._head

    @property
    def data(self):
        return self._data

    @property
    def plot(self):
        return self._plot

    @plot.setter
    def plot(self, value):
        raise AttributeError('Replacing plot is not allowed!')

    def hvplot(self):
        return self._hv_plot

    @property
    def hvplot(self):
        return self._hv_plot

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, value):
        if np.min(value) < 0 or np.max(value) > 1:
            raise ValueError('Minimal value must be larger or equal 0 and the maximal value must be smaller or equal 1')
        self._q = value

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def columns(self):
        if type(self._data) == list:
            out = []
            for d in self._data:
                out.extend(d.columns)
        else:
            return self._data.columns

    def to_astropy_table(self, category='table', one=True):
        """
        Returns the data as an astropy.table.Table

        :return: the data
        :rtype: astropy.table.Table
        """
        if one or type(self.data) != list:
            d = Table.from_pandas(self.combine())
            d.meta['category'] = category
        else:
            d = [Table.from_pandas(data) for data in self.data]

        return d

    def to_bin_table_hdu(self, category='table'):
        t = self.to_astropy_table(category, one=False)
        if type(t) == list:
            return [fits.BinTableHDU(k) for k in t]

        return fits.BinTableHDU(t)

    def _write(self, path, function):
        paths = []
        if type(self.data) == list:
            ending = path.split('.')[-1]
            path = path.split(f'.{ending}')[0]
            for d in self.data:
                name = d.survey_name
                p = f'{path}_{name}.{ending}'

                d.__getattr__(function)(p)
                paths.append(p)
        else:
            self.data.__getattr__(function)(path)
            paths.append(path)
        return paths

    def _to_sql(self, path, **kwargs):
        """
        Writes the data to sql database
        :param path: Path to the sql-database
        :type path: str
        :param kwargs: Additional arguments for pandas to_sql
        :return:
        """
        con = sqlite3.Connection(path)
        self._data.to_sql(self._category, con, **kwargs)
        return path

    def write(self, path, data_format='parquet', **kwargs):
        """
        Writes the data to a file

        :param path: Path to the save place
        :type path: str
        :param data_format:
            The format of the data-file. Current supported types are 'parquet', 'csv', 'sql', 'latex' and 'fits'.
        :return:
        """
        data_format = data_format.lower()
        if data_format == 'parquet':
            return self._write(path, 'to_parquet')
        elif data_format == 'csv':
            return self._write(path, 'to_csv')
        elif data_format == 'sql':
            return self._to_sql(path, **kwargs)
        elif data_format == 'latex':
            return self._write(path, 'to_latex')
        elif data_format == 'fits':
            return self.to_astropy_table().write(path, overwrite=True)
        else:
            raise ValueError(f'Format {data_format} is not supported.')

    def combine(self):
        """
        Combines all sub tables into one big table.
        If no sub tables exists, the data are returned
        :return:
        """
        if type(self.data) != list:
            return self.data
        out = None
        for d in self.data:
            if out is None:
                out = d.copy()
            else:
                out = out.join(d)
        return out

    def has_data(self):
        if self._data is not None:
            if type(self._data) == list and len(self._data) == 0:
                return False
            return True
        else:
            return False

    def __str__(self):
        return str(self._data)

    def __getitem__(self, item):
        if type(item) == str:
            return self._data[item]
        return self._data.loc[item]

    def __len__(self):
        return len(self._data)
