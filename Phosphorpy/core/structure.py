import pandas as pd


class Head:
    _name = None

    def __init__(self, name=None):
        self._name = name

    @property
    def name(self):
        return self._name


class Table:
    __data = None
    __head = None

    def __init__(self, d, name, mask):
        if d is None:
            d = pd.DataFrame()
        self.__data = d
        self.__head = Head(name=name)
        self.__mask = mask

    def __getitem__(self, item):
        return self.__data[item]

    def __setitem__(self, key, value):
        self.__data[key] = value

    def __len__(self):
        return len(self.__data)

    def __getattr__(self, item):
        return self.__data.__getattribute__(item)

    def __setattr__(self, key, value):
        if self.__data is None or '_' in key[0]:
            super().__setattr__(key, value)
        else:
            self.__data.__setattr__(key, value)

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
        self.__data.rename(name_map, axis)

    def merge(self, right, left_index=False, right_index=False):
        self.__data = self.__data.merge(right, left_index=left_index, right_index=right_index)

    @property
    def survey_name(self):
        return self.__head.name

    @property
    def mask(self):
        return self.__mask

    @property
    def data(self):
        return self.__data

    @property
    def survey(self):
        return self.__head
