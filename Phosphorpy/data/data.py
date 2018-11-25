from Phosphorpy.data.sub.magnitudes import Magnitude
from Phosphorpy.data.sub.colors import Color
from Phosphorpy.data.sub.coordinates import CoordinateTable
from Phosphorpy.external.vizier import query_by_name
from Phosphorpy.data.sub.plots.plot import Plot
from astropy.table import Table
from pandas import DataFrame
import numpy as np


class DataSet:
    _index = None
    _head = None
    _mask = None
    _coordinates = None
    _magnitudes = None
    _colors = None
    _flux = None
    _plot = None

    def __init__(self, data=None, index=None, coordinates=None, magnitudes=None, colors=None):
        """
        Standard data class for a survey like dataset. It requires a data file or at least coordinates and magnitudes.
        If no other data are given, it will try to compute the colors and index the sources.

        :param data: An input dataset in a numpy.ndarray style
        :type data: numpy.ndarray, astropy.table.Table, pandas.DataFrame
        :param index: A list of index, they must have the same length as coordinates, magnitudes and colors
        :type index: numpy.ndarray
        :param coordinates: A list with coordinates, they must have the same length as index, magnitudes and colors
        :type coordinates: numpy.ndarray
        :param magnitudes: A list with magnitudes, the must have the same length as index, coordinates and colors
        :type magnitudes: numpy.ndarray
        :param colors: A list with colors, they must have the same length as index, coordinates and magnitudes
        :type colors: numpy.ndarray
        """
        if data is not None:
            if type(data) == np.ndarray:
                cols = data.dtype.names
            elif type(data) == Table:
                cols = data.colnames
            elif type(data) == DataFrame:
                cols = data.columns
            else:
                raise TypeError('Unsupported data type: {}'.format(type(data)))

            if 'index' in cols:
                self._index = np.array(data['index'], dtype=np.int32)
            else:
                self._index = np.linspace(0, len(data), len(data), dtype=np.int32)

            self._coordinates = CoordinateTable(data)
            self._magnitudes = Magnitude(data)

        # if no data are given but coordinates, magnitudes and maybe colors
        elif coordinates is not None and magnitudes is not None and colors is not None:
            if coordinates.shape[0] != magnitudes.shape[0]:
                raise AttributeError('Coordinates and magnitudes do not have the same length!')
            if coordinates != type(CoordinateTable):
                coordinates = CoordinateTable(coordinates)

            self.coordinates = coordinates

            if magnitudes != type(Magnitude):
                magnitudes = Magnitude(magnitudes)

            self.magnitudes = magnitudes

            if type(colors) is bool and colors:
                self.colors = self.magnitudes.get_colors()
            else:
                if type(colors) != Color:
                    colors = Color(colors)

                self.colors = Color(colors)
            if index is not None:
                self.index = index
            else:
                self.index = np.arange(1, coordinates.shape[0], 1)
        else:
            raise AttributeError('data or at least coordinates and magnitudes are required!')

        self._plot = Plot(self)

    def __return_masked__(self, attribute):
        """
        Returns the input attribute with a applied mask, if a mask is set.

        :param attribute: One attribute of this class object
        :type attribute: CoordinateTable, Magnitude
        :returns: The attribute with an applied mask, if a mask is set. If no mask is set, the whole attribute
        :rtype: numpy.ndarray
        """
        if self._mask is not None:
            return attribute[self._mask.mask]
        else:
            return attribute

    @property
    def coordinates(self):
        return self.__return_masked__(self._coordinates)

    @coordinates.setter
    def coordinates(self, value):
        # todo: implement mask system
        self._coordinates = value

    @property
    def magnitudes(self):
        return self.__return_masked__(self._magnitudes)

    @magnitudes.setter
    def magnitudes(self, value):
        # todo: implement mask system
        self._magnitudes = value

    @property
    def colors(self):
        if self._colors is None:
            self._colors = self._magnitudes.get_colors()
        return self.__return_masked__(self._colors)

    @colors.setter
    def colors(self, value):
        # todo: implement mask system
        self._colors = value

    @property
    def flux(self):
        if self._flux is None:
            self._flux = self._magnitudes.get_flux()
        return self._flux

    @property
    def index(self):
        return self.__return_masked__(self._index)

    @index.setter
    def index(self, value):
        # todo: implement mask system
        self._index = value

    @property
    def plot(self):
        return self._plot

    def add_row(self, coordinate, magnitude, index=None, color=None):
        # todo: implement add row
        pass

    def remove_unmasked_data(self):
        """
        Removes all unmasked items from the dataset and sets the mask back to None.

        :return:
        """
        self._index = self.__return_masked__(self._index)
        self._magnitudes = self.__return_masked__(self._magnitudes)
        self._coordinates = self.__return_masked__(self._coordinates)
        self._colors = self.__return_masked__(self._colors)
        self._mask = None

    def __get_attribute__(self, item):
        """
        Returns the corresponding attribute to item, if the item doesn't match with
        one of the names, a KeyError will raise.

        :param item:
            The name of the attribute. Allowed strings are 'index', 'coordinate', 'coordinates',
            'magnitude', 'magnitudes', 'color', 'colors' and 'mask'
        :type item: str
        :returns: The corresponding attribute
        :rtype: numpy.ndarray
        """
        if item == 'index':
            return self.index
        elif item == 'coordinates' or item == 'coordinate':
            return self.coordinates
        elif item == 'magnitudes' or item == 'magnitude':
            return self.magnitudes
        elif item == 'colors' or item == 'color':
            return self.colors
        elif item == 'mask':
            return self._mask
        else:
            error = 'Key {} not known! Possible option are index, coordinates, magnitudes and colors.'.format(item)
            raise KeyError(error)

    def __get_row__(self, item):
        """
        Returns the values of the data element or row wise
        :param item: A slice or an int to indicate the required data
        :type item: slice, int
        :return:
        """
        if type(item) == slice:
            # todo: switch to DataTable representation
            out = self.coordinates.data[item].merge(self.magnitudes.data[item],
                                                    left_index=True,
                                                    right_index=True)
            if self.colors is not None:
                out = out.merge(self.colors[item],
                                left_index=True,
                                right_index=True)

            return out

        elif type(item) == int:
            # TODO: implement item wise returns
            pass

    def __getitem__(self, item):

        # if the item is a string, return one of the attributes of the object
        if type(item) == str:
            return self.__get_attribute__(item)
        # if the item is a slice or an int, return the corresponding data
        else:
            return self.__get_row__(item)

    def load_from_vizier(self, name):
        d = query_by_name(name, self.coordinates.to_table())
        self._magnitudes.add_survey_mags(d, name)
        return DataSet(d)

    def write(self, path):
        # todo: implement writing
        pass

    @staticmethod
    def read_from_file(name):
        # todo: implement from file
        pass

    @staticmethod
    def load_coordinates(path, format='fits', ra_name='ra', dec_name='dec'):
        """
        Creates a DataSet object from a file with coordinates.

        :param path: The path to the file
        :type path: str
        :param format: The format of the file
        :type format: str
        :param ra_name: The name of the RA column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :return: The DataSet object with the coordinates
        :rtype: DataSet
        """
        if format == 'fits':
            coords = Table.read(path)
            coords = coords[[ra_name, dec_name]]
            coords = coords.to_pandas()
        elif format == 'csv':
            coords = DataFrame.from_csv(path)
        else:
            raise ValueError('Format \'{}\' is not supported.'.format(format))

        return DataSet(coords)
