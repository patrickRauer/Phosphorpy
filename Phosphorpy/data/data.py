from Phosphorpy.data.sub.magnitudes import Magnitude, SurveyData
from Phosphorpy.data.sub.colors import Colors
from Phosphorpy.data.sub.coordinates import CoordinateTable
from Phosphorpy.data.sub.flux import FluxTable
from Phosphorpy.external.vizier import query_by_name, constrain_query
from Phosphorpy.data.sub.plots.plot import Plot
from Phosphorpy.data.sub.table import Mask
from Phosphorpy.external.image import PanstarrsImage, SDSSImage
from astropy.table import Table
from astropy.io import fits
import pandas as pd
from pandas import DataFrame
import numpy as np
import zipfile
import os


def add_to_zip(zi, data, name, format='csv'):
    """
    Writes a file to zip file

    :param zi: The zip-file object

    :param data: The data file to write
    :type data: Phosphorpy.data.sub.table.DataTable, Phosphorpy.data.sub.magnitudes.Survey
    :param name: Name of the file
    :type name: str
    :param format: The format of the file
    :type format: str
    :return:
    """
    data.write('./{}'.format(name), data_format=format)
    zi.write('./{}'.format(name), name)
    os.remove('./{}'.format(name))


def read_from_zip(zi, name):
    """
    Reads a file from a zip archive and converts to a pandas DataFrame

    :param zi: The zip-file object
    :type zi: zipfile.ZipFile
    :param name: Name of the member file
    :type name: str
    :return: The red item
    :rtype: pandas.DataFrame
    """
    zi.extract(name, '.')
    ending = name.split('.')[-1]
    if ending == 'fits':
        d = Table.read(name).to_pandas()
    elif ending == 'csv':
        d = pd.read_csv(name)
    elif ending == 'ini':
        return SurveyData.read(name)
    else:
        raise ValueError('Format {} is not supported.'.format(format))
    os.remove(name)
    return d


class DataSet:
    _index = None
    _head = None
    _mask = None
    _coordinates = None
    _magnitudes = None
    _colors = None
    _flux = None
    _plot = None

    def __init__(self, data=None, index=None, coordinates=None, magnitudes=None, colors=None, flux=None,
                 survey=None):
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
        self._mask = Mask()
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
            self._coordinates = CoordinateTable(data, mask=self._mask)
            self._magnitudes = Magnitude(data, mask=self._mask)

        # if no data are given but coordinates, magnitudes and maybe colors
        # todo: coordinates and magnitudes are the minimal requirement, colors and flux are optional
        elif coordinates is not None and magnitudes is not None and colors is not None:
            if coordinates.shape[0] != magnitudes.shape[0]:
                raise AttributeError('Coordinates and magnitudes do not have the same length!')
            if type(coordinates) != CoordinateTable:
                coordinates = CoordinateTable(coordinates, mask=self._mask)

            self.coordinates = coordinates

            if type(magnitudes) != Magnitude:
                magnitudes = Magnitude(magnitudes, mask=self._mask)

            self.magnitudes = magnitudes

            if type(colors) is bool and colors:
                self.colors = self.magnitudes.get_colors()
            else:
                if type(colors) != Colors:
                    self.colors = Colors(colors, mask=self._mask)

            if type(flux) != FluxTable:
                flux = FluxTable(flux)

            self._flux = flux

            if index is not None:
                self.index = index
            else:
                self.index = np.arange(1, len(coordinates), 1)
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
        if self._mask.get_mask_count() > 0:
            return attribute[self._mask.mask]
        else:
            return attribute

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        # todo: implement mask system
        self._coordinates = value

    @property
    def magnitudes(self):
        """

        :return: MagnitudeTable of the DataSet
        :rtype: Magnitude
        """
        return self._magnitudes

    @magnitudes.setter
    def magnitudes(self, value):
        # todo: implement mask system
        self._magnitudes = value

    @property
    def colors(self):
        if self._colors is None:
            self._colors = self._magnitudes.get_colors()
        return self._colors

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
        return self._index

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
        """
        Load new photometric data from Vizier

        :param name: Name of the survey
        :type name: str
        :return:
        """
        d = query_by_name(name, self.coordinates.to_table())
        self._magnitudes.add_survey_mags(d, name)
        # return DataSet(d)

    def images(self, survey, source_id, directory=''):
        """
        Download images from SDSS or Pan-STARRS and create a colored image out of them

        :param survey: SDSS or Pan-STARRS
        :type survey: str
        :param source_id: The id of the source
        :type source_id: int
        :return:
        """
        survey = survey.lower()
        if survey == 'sdss':
            s = SDSSImage()
        elif survey == 'panstarrs' or survey == 'ps' or survey == 'ps1' or survey == 'pan-starrs':
            s = PanstarrsImage()
        else:
            raise ValueError('{} is not allowed. Only SDSS or PANSTARRS as surveys are allowed.')
        coord = self.coordinates.as_sky_coord(source_id)
        if directory != '':
            directory = directory + coord.to_string('hmsdms')+'.png'
        s.get_color_image(coord, directory)

    def all_images(self, survey, directory=''):
        """
        Downloads images of all sources and create a colored image for every source.

        :param survey: SDSS or Pan-STARRS
        :type survey: str
        :param directory: The path to the directory, where the images should be stored
        :type directory: str
        :return:
        """
        m = self._mask.get_latest_mask()
        for i in range(len(self.coordinates)):
            if not m[i]:
                continue
            self.images(survey, i, directory)

    def __combine_all__(self):
        """
        Combines all tables to one big table

        :return: The combined table
        :rtype: pandas.DataFrame
        """
        comb = self.coordinates.to_table().join(self.magnitudes, how='outer')
        if self._flux is not None:
            comb = comb.join(self._flux, how='outer')

        if self._colors is not None:
            comb = comb.join(self._colors, how='outer')
        return comb

    def __write_as_zip__(self, path, format):
        """
        Creates a zip file with the data of the DataSet

        :param path: Path to the save place
        :type path: str
        :param format: Format of the internal data-files
        :type format: str
        :return:
        """
        with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as zi:
            add_to_zip(zi, self.coordinates, 'coordinates.{}'.format(format), format=format)
            add_to_zip(zi, self.magnitudes, 'magnitudes.{}'.format(format), format=format)
            add_to_zip(zi, self.magnitudes.survey, 'survey.ini', format='init')
            if self._colors is not None:
                add_to_zip(zi, self.colors, 'colors.{}'.format(format), format=format)
            if self._flux is not None:
                add_to_zip(zi, self.flux, 'flux.{}'.format(format), format=format)
                # todo: add plots to the zip file if some of them were made.

    def __write_as_fits__(self, path):
        """
        Writes the current table into one fits file.

        :param path: The path to the save place
        :type path: str
        """
        # todo: change to multi layer
        hdu_list = [self.coordinates.to_astropy_table('coordinates'),
                    self.magnitudes.to_astropy_table('magnitudes')]

        if self._colors is not None:
            hdu_list.append(self.colors.to_astropy_table('colors'))
        if self._flux is not None:
            hdu_list.append(self.flux.to_astropy_table('flux'))
        hdu_list = fits.HDUList(hdu_list)
        hdu_list.writeto(path, overwrite=True)
        # combined = self.__combine_all__()
        # combined = Table.from_pandas(combined)
        # combined.write(path, format='fits')

    def __write_as_csv__(self, path):
        """
        Writes the current table into one csv file.

        :param path: The path to the save place
        :type path: str
        """
        combined = self.__combine_all__()
        combined.to_csv(path)

    def __write_as_report__(self, path):
        # todo: implement report writing
        a = self.plot
        raise AttributeError('Report is not implemented yet')
        pass

    def write(self, path, format='zip', in_zip_format='fits'):
        """
        Writes the current data to a file (-system)

        :param path: The path to the save place
        :type path: str
        :param format: The format of the saving (zip, fits, csv or report)
        :type format: str
        :param in_zip_format: Format of the data-files in the zip-file. Only required, if format is 'zip'
        :type in_zip_format: str
        :return:
        """
        format = format.lower()
        if format == 'zip':
            self.__write_as_zip__(path, in_zip_format)
        elif format == 'fits':
            self.__write_as_fits__(path)
        elif format == 'csv':
            self.__write_as_csv__(path)
        elif format == 'report':
            self.__write_as_report__(path)
        else:
            raise ValueError('Format {} is not supported.'.format(format))

    def __str__(self):
        s = 'Number of entries: {}\n'.format(len(self.coordinates))
        s += 'Available surveys:\n'
        if self._magnitudes is not None:
            if self._magnitudes._survey is not None:
                for su in self._magnitudes._survey.get_surveys():
                    s += '\t{}\n'.format(su)
            else:
                s += '\tno survey data set\n'
        else:
            s += '\tNone\n'
        if self._colors is not None:
            s += 'Colors computed: True\n'
        else:
            s += 'Colors computed: False\n'
        if self._flux is not None:
            s += 'Flux computed: True\n'
        else:
            s += 'Flux computed: False\n'
        return s

    @staticmethod
    def read_from_file(name, format='zip'):
        """
        Reads a DateSet from a file

        :param name: The name/path of the file
        :type name: str
        :param format: The format of the file
        :type format: str
        :return: The DataSet from the data of the file
        :rtype: DataSet
        """
        if format == 'zip':
            d = {}
            head = None
            with zipfile.ZipFile(name) as zip_file:
                for file_name in zip_file.namelist():
                    if '.ini' not in file_name:
                        d[file_name.split('.')[0]] = read_from_zip(zip_file, file_name)
                    else:
                        head = read_from_zip(zip_file, file_name)
            ds = DataSet(**d)
            ds._magnitudes._survey = head
            return ds
        elif format == 'fits':
            with fits.open(name) as fi:
                d = {}
                for i, f in enumerate(fi):
                    if i > 0:
                        data = f.data
                        head = f.header
                        d[head['category']] = data
                return DataSet(**d)

        elif format == 'csv':
            raise ValueError('Automatic csv reading not implemented')
            # todo: implement csv reading
            pass
        else:
            raise ValueError('Format {} is not supported.'.format(format))

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

    @staticmethod
    def from_vizier(name, **constrains):
        """
        Stars a constrain query on a vizier catalog and convert the output to a DataSet object.

        :param name: The name of the survey
        :type name: str
        :param constrains:
            A dict of constrains with the key as the column name in the vizier catalog and the corresponding
            value as the actual constrain
        :return: A DataSet object with the results of the query
        :rtype: DataSet
        """
        d = constrain_query(name, **constrains)
        d = d.to_pandas()
        ds = DataSet(d)
        ds._magnitudes.add_survey_mags(d, name)
        return ds
