from astropy.table import Table
import numpy as np


class Mask:
    _mask = []
    _desc = []

    def __init__(self):
        pass

    def add_mask(self, mask, description='', combine=True):
        """
        Adds a new mask to the storage.

        :param mask:
            The new mask with the size of the complete dataset or with the size of passed rows in the previous mask.
        :type mask: numpy.ndarray
        :param description: The description of the mask. Default is an empty string.
        :type description: str
        :param combine: True if the previous mask should be used (have True values), too, else False. Default is True.
        :type combine: bool
        :return:
        """
        if combine:
            if len(self._mask) > 0:
                mask = mask & self._mask[-1]
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

    @property
    def mask(self):
        return self.get_latest_mask()

    @property
    def description(self):
        return self.get_latest_description()


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
            self.data.to_parquet(path)
        elif data_format == 'csv':
            self.data.to_csv(path)
        elif data_format == 'sql':
            self.data.to_sql(path)
        elif data_format == 'latex':
            self.data.to_latex(path)
        elif data_format == 'fits':
            Table.from_pandas(self.data).write(path)
        else:
            raise ValueError('Format {} is not supported.'.format(data_format))

    def __str__(self):
        return str(self._data)
