from astropy.table import Table as AstroTable
import numpy as np
import pandas as pd
import os


class Mask:
    _mask = None
    _desc = None

    def __init__(self, length, ids=None):
        if ids is None:
            ids = np.arange(length)+1
        elif len(ids) != length:
            raise AttributeError('"length" and the length of "ids" must be the same value.')

        self._mask = [pd.Series(length*[True], ids)]
        self._desc = ['initialization']

    def __str__(self):
        string = ''
        for i in range(self.get_mask_count()):
            string += f'Mask No {i}: {self.get_description(i)}\n'
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
        if type(mask) != pd.Series:
            raise ValueError(f'Mask must be a pandas Series but mask has the type {type(mask)}')

        if len(mask) == 0:
            raise ValueError('Must must have at least one element.')
        if combine:
            if len(self._mask) > 0:
                # Aligns the new incoming mask with the previous mask
                # because if a survey doesn't contain all sources the incoming mask
                # will have different size
                # join left because the original mask is the largest by definition
                # fill_value is True because missing values can not be selected
                mask = self._mask[-1].align(mask, fill_value=True, join='left')[1]
                mask &= self._mask[-1]
        else:
            if len(self._mask[0]) != len(mask):
                raise ValueError('New mask must have the same length as the initial one.')

        self._mask.append(mask)
        self._desc.append(description)

    def get_latest_mask(self):
        """
        Returns the latest mask.

        :return: The latest mask
        :rtype: pd.Series
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
        :rtype: pd.Series
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
        self._mask.remove(self._mask[mask_id])
        self._desc.remove(self._desc[mask_id])

    def reset_mask(self):
        """
        Deletes all masks and sets the mask back to initialization status
        :return:
        """
        self._mask = self._mask[:1]
        self._desc = self._desc[:1]

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

    def __len__(self):
        return len(self._mask)


class Head:
    _name = None

    def __init__(self, name=None):
        if name is None:
            name = ''
        self._name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value is not None and type(value) == str:
            self._name = value
        else:
            raise ValueError('New name must be a string.')


class Table:
    __data = None
    __head = None
    __mask = None

    def __init__(self, d, name, mask):
        if d is None:
            d = pd.DataFrame()
        self.__data = d
        self.__head = Head(name=name)

        if mask is None:
            mask = Mask(len(d))
        self.__mask = mask

    def __getitem__(self, item):
        return self.__data[item]

    def __setitem__(self, key, value):
        self.__data.loc[:, key] = value

    def __len__(self):
        return len(self.__data)

    # def __getattr__(self, item):
    #     if item in dir(self):
    #         return self.__getattr__(item)
    #     return self.__data.__getattribute__(item)

    # def __setattr__(self, key, value):
    #     print('key', key)
    #     if self.__data is None or '_' in key[0] or 'data' in key:
    #         print('to class')
    #         super().__setattr__(key, value)
    #     else:
    #         print('to data')
    #         self.__data.__setattr__(key, value)

    def set_mask(self, mask):
        """
        Sets a new mask to the table

        :param mask: The new mask

        :return:
        """
        self.__mask = mask

    def select_columns(self, columns):
        """
        Select a subset of columns

        :param columns: The name of the columns
        :type columns: Union
        :return:
        """
        self.__data = self.__data[columns]

    def rename(self, name_map, axis=None):
        """
        Rename of columns

        :param name_map: The rename map
        :type name_map: dict
        :param axis:
        :return:
        """
        self.__data.rename(name_map, axis=axis)

    def merge(self, right, left_index=False, right_index=False):
        self.__data = self.__data.merge(right, left_index=left_index, right_index=right_index)

    def has_name(self, name):
        return name == self.survey_name or name.lower() == self.survey_name

    @property
    def survey_name(self):
        return self.__head.name

    @property
    def mask(self):
        return self.__mask

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, d):
        self.__data = d

    @property
    def survey(self):
        return self.__head

    @property
    def columns(self):
        return self.__data.columns

    @property
    def values(self):
        return self.__data.values

    @property
    def index(self):
        return self.__data.index

    def copy(self):
        return self.__data.copy()

    def write(self, path, data_format='parquet', **kwargs):
        path, extension = os.path.splitext(path)
        path = f'{path}_{self.survey_name}.{extension}'
        if data_format == 'parquet':
            self.data.to_parquet(path)
        elif data_format == 'sql':
            self.data.to_sql(self.survey_name, kwargs['connection'])
        elif data_format == 'csv':
            self.data.to_csv(path)
        elif data_format == 'fits':
            t = AstroTable.from_pandas(self.data)
            try:
                self.survey.add_to_meta(t.meta, self.survey_name)
            except AttributeError:
                pass
            t.write(path, overwrite=True)
        else:
            raise ValueError('Unknown format')
        return path
