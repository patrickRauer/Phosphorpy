from astropy.table import Table
import numpy as np
import pandas as pd


class Mask:
    _mask = None
    _desc = None

    def __init__(self, length):
        self._mask = [pd.Series(length*[True], np.arange(length))]
        self._desc = ['initialization']

    def __str__(self):
        string = ''
        for i in range(self.get_mask_count()):
            string += 'Mask No {}: {}\n'.format(i, self.get_description(i))
        return string

    def add_mask(self, mask, description='', combine=True):
        """
        Adds a new mask to the storage.

        :param mask:
            The new mask with the size of the complete dataset or with the size of passed rows in the previous mask.
        :type mask: pandas.Series
        :param description: The description of the mask. Default is an empty string.
        :type description: str
        :param combine: True if the previous mask should be used (have True values), too, else False. Default is True.
        :type combine: bool
        :return:
        """
        if combine:
            if len(self._mask) > 0:
                # Aligns the new incoming mask with the previous mask
                # because if a survey doesn't contain all sources the incoming mask
                # will have different size
                # join left because the original mask is the largest by definition
                # fill_value is True because missing values can not be selected
                mask = self._mask[-1].align(mask, fill_value=True, join='left')[1]
                mask &= self._mask[-1]

        self._mask.append(mask)
        self._desc.append(description)

    def get_latest_mask(self):
        """
        Returns the latest mask.

        :return: The latest mask
        :rtype: numpy.ndarray
        """
        return self.get_mask(-1)

    def get_latest_description(self):
        """
        Returns the latest description.

        :return: The description text
        :rtype: int
        """
        return self.get_description(-1)

    def get_mask(self, level):
        """
        Returns the description of the mask at the corresponding level.

        :param level: The level of the mask.
        :type level: int
        :return: The mask of the data
        :rtype: numpy.ndarray
        """
        return self._mask[level]

    def get_description(self, level):
        """
        Returns the description of the mask at the corresponding level.

        :param level: The level of the mask.
        :type level: int
        :return: The description of the mask
        :rtype: str
        """
        return self._desc[level]

    def get_mask_count(self):
        """
        Returns the number of masks.

        :return: The number of masks
        :rtype: int
        """
        return len(self._mask)

    def remove_mask(self, mask_id):
        """
        Removes a mask from the list.

        :param mask_id:
        :return:
        """
        pass

    def reset_mask(self):
        """
        Deletes all masks
        :return:
        """
        self._mask = []
        self._desc = []

    @property
    def mask(self):
        return self.get_latest_mask()

    @property
    def description(self):
        return self.get_latest_description()

    def __call__(self, *args, **kwargs):
        # if no mask was set
        if len(self._mask) == 0:
            # return the input data
            return args[0]
        # if at least one mask was set
        else:
            # return the mask data
            return args[0][self.get_latest_mask()]


class DataTable:

    _mask = None
    _head = None
    _data = None

    _plot = None
    _q = [0.15, 0.25, 0.75, 0.85]

    def __init__(self, mask=None):
        """
        Basic data table class
        """
        self._mask = mask

    def set_mask(self, mask):
        self._mask = mask
        if type(self.data) == list:
            print(self.data)
            for d in self.data:
                print(type(d))
                d.set_mask(mask)

    def stats(self):
        """
        Returns basic statistics (mean, median, std, min, max) of the magnitudes

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
        self.apply(func)

    def remove_unmasked_data(self):
        """
        Removes all unmasked (mask[i] == False) from the data
        :return:
        """
        if self._data is not None:
            self._data = self._data[self._mask.get_latest_mask()]

    def select_nan(self, column):
        """
        Select all rows with a NaN value
        :return:
        """
        d = self._data[column].values
        self.mask.add_mask(d == d, f'Mask all NaN values in {column}')

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

    def to_astropy_table(self, category='table'):
        """
        Returns the data as an astropy.table.Table

        :return: the data
        :rtype: astropy.table.Table
        """
        d = Table.from_pandas(self._data)
        d.meta['category'] = category
        return d

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

    def write(self, path, data_format='parquet'):
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
            return self._write(path, 'to_sql')
        elif data_format == 'latex':
            return self._write(path, 'to_latex')
        elif data_format == 'fits':
            return self._write(path, 'fits')
        else:
            raise ValueError(f'Format {data_format} is not supported.')

    def __str__(self):
        return str(self._data)

    def __getitem__(self, item):
        return self._data[item]
