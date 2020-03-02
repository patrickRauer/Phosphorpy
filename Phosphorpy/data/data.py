from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from http.client import RemoteDisconnected

from Phosphorpy.data.sub.colors import Colors
from Phosphorpy.data.sub.coordinates import CoordinateTable
from Phosphorpy.data.sub.flux import FluxTable
from Phosphorpy.data.sub.light_curve import LightCurves
from Phosphorpy.data.sub.magnitudes import MagnitudeTable as Magnitude, SurveyData
from Phosphorpy.data.sub.plots.plot import Plot
from Phosphorpy.data.sub.table import Mask
from Phosphorpy.external.image import PanstarrsImage, SDSSImage
from Phosphorpy.external.spectra import get_spectra, get_all_spectra
from Phosphorpy.external.vizier import query_by_name, query_simbad, constrain_query
from Phosphorpy.data.sub.astrometry import AstrometryTable

try:
    from extinction.extinction import get_extinctions
except ImportError:
    from Phosphorpy.external.extinction import get_extinctions

import pandas as pd
from pandas import DataFrame
import numpy as np
import zipfile
import os
import warnings
import sys


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
    p = f'./{name}'
    paths = data.write(p, data_format=format)
    if paths is None:
        zi.write(p, name)
        os.remove(p)
    else:
        for p in paths:
            zi.write(p, p.split('./')[-1])
            os.remove(p)


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
        d = pd.read_csv(name, index_col=0)
    elif ending == 'ini':
        return SurveyData.read(name)
    else:
        raise ValueError(f'Format {format} is not supported.')
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
    _astrometry = None

    _light_curves = None
    _spectra = None
    _plot = None

    def __init__(self, data=None, index=None, coordinates=None, magnitudes=None, colors=None, flux=None):
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
        try:
            self._mask = Mask(len(data))
        except TypeError:
            self._mask = Mask(len(coordinates))
        if data is not None:
            if type(data) == np.ndarray:
                cols = data.dtype.names
            elif type(data) == Table:
                cols = data.colnames
            elif type(data) == DataFrame:
                cols = data.columns
            else:
                raise TypeError(f'Unsupported data type: {type(data)}')

            if 'index' in cols:
                self._index = np.array(data['index'], dtype=np.int32)
            else:
                self._index = np.linspace(0, len(data), len(data), dtype=np.int32)
            self._coordinates = CoordinateTable(data, mask=self._mask)
            self._magnitudes = Magnitude(data, mask=self._mask)

        # if no data are given but coordinates, magnitudes and maybe colors
        elif coordinates is not None and magnitudes is not None:
            # if coordinates.shape[0] != magnitudes.shape[0]:
            #     raise AttributeError('Coordinates and magnitudes do not have the same length!')
            if type(coordinates) != CoordinateTable:
                coordinates = CoordinateTable(coordinates, mask=self._mask)

            self._coordinates = coordinates

            if type(magnitudes) != Magnitude:
                magnitudes = Magnitude(magnitudes, mask=self._mask)

            self._magnitudes = magnitudes

            if colors is not None:
                if type(colors) is bool and colors:
                    self._colors = self.magnitudes.get_colors()
                else:
                    if type(colors) != Colors:
                        self._colors = Colors(colors, mask=self._mask)
                    else:
                        self._colors = colors

            if flux is not None:
                if type(flux) != FluxTable:
                    flux = FluxTable(flux)

                self._flux = flux

            if index is not None:
                self._index = index
            else:
                self._index = np.arange(1, len(coordinates), 1)
        elif coordinates is not None:
            self._coordinates = CoordinateTable(coordinates, mask=self._mask)
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

    def set_astrometry(self, astrometry):
        self._astrometry = astrometry

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def magnitudes(self):
        """

        :return: MagnitudeTable of the DataSet
        :rtype: Magnitude
        """
        return self._magnitudes

    @property
    def colors(self):
        if self._colors is None:
            self._colors = self._magnitudes.get_colors()
        return self._colors

    @property
    def flux(self):
        if self._flux is None:
            self._flux = self._magnitudes.get_flux()
        return self._flux

    @property
    def index(self):
        return self._index

    @property
    def plot(self):
        return self._plot

    @property
    def astrometry(self):
        if self._astrometry is None:
            self._astrometry = AstrometryTable.load_astrometry(self._coordinates)
        return self._astrometry

    @property
    def spectra(self):
        if self._spectra is None:
            self._spectra = get_all_spectra(
                self.coordinates.to_sky_coord(),
                self.coordinates.data.index.values
            )
        return self._spectra

    @property
    def mask(self):
        return self._mask

    def get_spectra(self, coordinate=None, ra=None, dec=None):
        """
        Returns the spectra close to the given coordinates.

        :param coordinate: Coordinates (or close to) of the required ID.
        :type coordinate: SkyCoord
        :param ra: The RA component of the coordinates
        :type ra: float
        :param dec: The Dec component of the coordinates
        :type dec: float
        :return: The closest spectra
        :rtype: SpectraList
        """
        min_id, distance = self.coordinates.get_closest_source_id(coordinate, ra, dec)
        return self.spectra.get_by_id(min_id)

    def get_light_curve(self, coordinate=None, ra=None, dec=None):
        """
        Returns the light curve close to the given coordinates.

        :param coordinate: Coordinates (or close to) of the required ID.
        :type coordinate: SkyCoord
        :param ra: The RA component of the coordinates
        :type ra: float
        :param dec: The Dec component of the coordinates
        :type dec: float
        :return: The closest light curve
        :rtype: LightCurves
        """
        min_id, distance = self.coordinates.get_closest_source_id(coordinate, ra, dec)
        return self.light_curves.get_light_curve(min_id)

    def remove_unmasked_data(self):
        """
        Removes all unmasked items from the dataset and sets the mask back to None.

        :return:
        """
        # self._index = self.__return_masked__(self._index)

        self._coordinates.remove_unmasked_data()
        if self._magnitudes is not None:
            self._magnitudes.remove_unmasked_data()
        if self._colors is not None:
            self._colors.remove_unmasked_data()
        if self._flux is not None:
            self._colors.remove_unmasked_data()
        if self._astrometry is not None:
            self._astrometry.remove_unmasked_data()
        self._coordinates.mask.reset_mask()

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
            error = f'Key {item} not known! Possible option are index, coordinates, magnitudes and colors.'
            raise KeyError(error)

    def __get_row__(self, item):
        """
        Returns the values of the data element or row wise
        :param item: A slice or an int to indicate the required data
        :type item: slice, int
        :return:
        """
        if type(item) == slice:
            out = self.coordinates.data[item].merge(self.magnitudes.data[item],
                                                    left_index=True,
                                                    right_index=True)
            if self.colors is not None:
                out = out.merge(self.colors[item],
                                left_index=True,
                                right_index=True)

            return out

        elif type(item) == int:
            out = self.coordinates[item]
            if self._magnitudes.has_data():
                out = out.concat(self.magnitudes[item], axis=1)
            if self._colors is not None and self._colors.has_data():
                out = out.concate(self.colors[item], axis=1)

            if self._flux is not None and self._flux.has_data():
                out = out.concate(self.flux[item], axis=1)
            return out

    def __getitem__(self, item):

        # if the item is a string, return one of the attributes of the object
        if type(item) == str:
            return self.__get_attribute__(item)
        # if the item is a slice or an int, return the corresponding data
        else:
            return self.__get_row__(item)

    def __load_from_vizier__(self, name):
        d = query_by_name(name, self.coordinates.to_table())
        self._magnitudes.add_survey_mags(d, name.lower())
        self._colors = None
        self._flux = None

    def load_from_vizier(self, name):
        """
        Load new photometric data from Vizier.
        At the moment the following surveys are available

        * optical (also keyword for all optical surveys)
            * SDSS
            * Pan-STARRS
            * KiDS
            * GAIA
            * APASS (not included in 'optical')
        * NIR (also keyword for all NIR surveys)
            * 2MASS
            * VIKING
            * UKIDSS
        * UV (no keyword, use instead GALEX directly)
            * GALEX
        * IR (no keyword, use instead WISE directly)
            * WISE

        .. code:: python
            # ds is a DataSet object with a set of coordinates
            # download the available photometric data from the SDSS survey
            ds.load_from_vizier('sdss')

            # load also all NIR photometric data
            ds.load_from_vizier('NIR')

        :param name:
            A single name or a list of names of surveys or one of the keywords 'optical' or 'nir' to
            download all optical or NIR surveys
        :type name: str, list, tuple
        :return:
        """
        if type(name) == str:
            if name.lower() == 'optical':
                self.load_from_vizier(['sdss', 'ps', 'kids', 'gaia'])
            elif name.lower() == 'nir':
                self.load_from_vizier(['2mass', 'ukidss', 'viking'])
            elif name.lower() == 'all':
                self.load_from_vizier(['galex', 'optical', 'nir', 'wise'])
            else:
                self.__load_from_vizier__(name)
        elif type(name) == list or type(name) == tuple:
            for n in name:
                self.load_from_vizier(n)

    def load_spectra(self, survey):
        """
        Loads all available spectra for the set coordinates

        :param survey:
            The name of the survey from which the spectra should be taken.
            Currently, available surveys are 'SDSS', 'LAMOST' and 'GAMA'
        :return:
        """
        specs = get_spectra(self.coordinates.to_sky_coord(), survey)

        if self._spectra is None:
            self._spectra = specs
        else:
            self._spectra.merge(specs)

        return self._spectra

    def add_magnitudes(self, data, survey_info, ra_name='ra', dec_name='dec'):
        """
        Adds magnitudes from a different source to the DataSet

        :param data: The new data with coordinates and the magnitudes
        :type data: pd.DataFrame, Table
        :param survey_info:
        :param ra_name: The name of the R.A. column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :return:
        """
        # check the type of data
        if type(data) != (pd.DataFrame or Table):
            raise ValueError('pandas DataFrame or astropy Table is required')
        elif type(data) == Table:
            data = data.to_pandas()

        # rename the columns if they aren't the default names
        if ra_name != 'ra' or dec_name != 'dec':
            data = data.rename({ra_name: 'ra', dec_name: 'dec'})
        # TODO: implement
        data = self.coordinates.match(data)
        pass

    def get_simbad_data(self):
        return query_simbad(self.coordinates)

    def images(self, survey, source_id, directory='', bands=None, size=None, smooth=0):
        """
        Download images from SDSS or Pan-STARRS and create a colored image out of them

        :param survey: SDSS or Pan-STARRS
        :type survey: str
        :param source_id: The id of the source
        :type source_id: int
        :param directory: The path to the directory, where the images should be saved
        :type directory: str
        :param bands: Individual filter combination. Default is None, which means the default combination is used.
        :type bands: None, tuple, list
        :param size: The wanted size of the image
        :type size: float, astropy.units.Quantity
        :param smooth: Number of smooths. Default is 0.
        :type smooth: int
        :return:
        """
        survey = survey.lower()

        # if sdss is the selected survey
        if survey == 'sdss':
            s = SDSSImage()
        # if pan-starrs is the selected survey
        elif survey == 'panstarrs' or survey == 'ps' or survey == 'ps1' or survey == 'pan-starrs':
            s = PanstarrsImage()
        else:
            raise ValueError('{} is not allowed. Only SDSS or PANSTARRS as surveys are allowed.')
        coord = self.coordinates.to_sky_coord(source_id)

        # if the directory string is not empty
        if directory != '':
            # if the last character is not '/'
            if directory[-1] != '/':
                # add '/' to the end of the directory string
                directory = directory + '/'
            directory = directory + coord.to_string('hmsdms')+'.png'
        try:
            s.get_color_image(coord, directory, bands=bands, size=size, smooth=smooth)
        except ValueError:
            warnings.warn("Image is not available", UserWarning)
        except RemoteDisconnected:
            warnings.warn('Connection interrupted.')

    def all_images(self, survey, directory='', bands=None, size=None, smooth=0):
        """
        Downloads images of all sources and create a colored image for every source.

        :param survey: SDSS or Pan-STARRS
        :type survey: str
        :param directory: The path to the directory, where the images should be stored
        :type directory: str
        :param bands: Individual filter combination. Default is None, which means the default combination is used.
        :type bands: None, tuple, list
        :param size: The wanted size of the image
        :type size: float, astropy.units.Quantity
        :param smooth: Number of smoothings. Default is 0.
        :type smooth: int
        :return:
        """
        try:
            m = self._mask.get_latest_mask().values
        except IndexError:
            m = np.ones(len(self.coordinates))
        for i in range(len(self.coordinates)):
            if not m[i]:
                continue
            self.images(survey, i, directory, bands=bands, size=size, smooth=smooth)

    @property
    def light_curves(self):
        """
        Downloads all light curves from CSS
        :return:
        """
        if self._light_curves is None:
            self._light_curves = LightCurves(self._coordinates)
        return self._light_curves

    def correct_extinction(self):
        """
        Estimates and applies the extinction correction based on
        Schlafly & Finkbeiner 2011 (ApJ 737, 103) to the magnitudes.
        If fluxes are colors are already computed, then they will be computed
        again after the correction.
        :return:
        """
        if self.magnitudes.data is None:
            raise AttributeError('No magnitudes are set.')
        for s in self.magnitudes.survey.get_surveys():
            extinc = get_extinctions(self.coordinates['ra'],
                                     self.coordinates['dec'],
                                     wavelengths=self.magnitudes.survey.get_survey_wavelengths(s),
                                     filter_names=self.magnitudes.survey.get_survey_magnitude(s))
            self.magnitudes.apply_extinction_correction_to_survey(extinc, s)
        # extinc = get_extinctions(self.coordinates['ra'],
        #                          self.coordinates['dec'],
        #                          wavelengths=self.magnitudes.survey.get_all_wavelengths(),
        #                          filter_names=self.magnitudes.survey.all_magnitudes())
        # self.magnitudes.apply_extinction_correction(extinc)

        # recompute flux and colors if they are already computed
        if self._flux is not None:
            self._flux = None
            self._flux = self._magnitudes.get_flux()

        if self._colors is not None:
            self._colors = None
            self._colors = self._magnitudes.get_colors()

    def reset_masks(self):
        """
        Resets the mask

        :return:
        """
        if self._magnitudes is not None:
            self.magnitudes.set_mask(self._mask)
        if self._flux is not None:
            self.flux.set_mask(self._mask)
        if self._colors is not None:
            self.colors.set_mask(self._mask)
        self.coordinates.set_mask(self._mask)

    def __combine_all__(self):
        """
        Combines all tables to one big table

        :return: The combined table
        :rtype: pandas.DataFrame
        """
        comb = self.coordinates.to_table()
        if self.magnitudes is not None and self.magnitudes.data is not None:
            comb = comb.join(self.magnitudes.combine(), how='outer')
        if self._flux is not None:
            fluxes = self.flux.combine()
            fluxes = fluxes.rename(
                dict(zip(
                    fluxes.columns.values,
                    [f.replace('mag', 'flux') for f in fluxes.columns.values]
                )),
                axis='columns'
            )

            comb = comb.join(fluxes, how='outer')

        if self._colors is not None:
            comb = comb.join(self.colors.combine(), how='outer')
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
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as zi:
            add_to_zip(zi, self.coordinates, f'coordinates.{format}', format=format)
            if self._magnitudes.data is not None:
                add_to_zip(zi, self.magnitudes, f'magnitudes.{format}', format=format)
                add_to_zip(zi, self.magnitudes.survey, 'survey.ini', format='init')

            if self._colors is not None:
                add_to_zip(zi, self.colors, f'colors.{format}', format=format)

            if self._flux is not None:
                add_to_zip(zi, self.flux, f'flux.{format}', format=format)

            if self._light_curves is not None:
                pass

            if self._spectra is not None:
                add_to_zip(zi, self.spectra, f'spectra/', format=format)

    def __write_as_fits__(self, path):
        """
        Writes the current table into one fits file.

        :param path: The path to the save place
        :type path: str
        """
        hdu_list = [fits.PrimaryHDU(),
                    self.coordinates.to_bin_table_hdu('coordinates')]

        if self.magnitudes is not None and self.magnitudes.data is not None:
            hdu_list.extend(self.magnitudes.to_bin_table_hdu('magnitudes'))

        if self._colors is not None:
            hdu_list.extend(self.colors.to_bin_table_hdu('colors'))

        if self._flux is not None:
            hdu_list.extend(self.flux.to_bin_table_hdu('flux'))

        if self._light_curves is not None:
            hdu_list.extend(self.light_curves.to_bin_table_hdu())

        if self._spectra is not None:
            self._spectra.write(f'{os.path.dirname(path)}/spectra/')

        hdu_list = fits.HDUList(hdu_list)
        hdu_list.writeto(path, overwrite=True)

    def __write_as_csv__(self, path):
        """
        Writes the current table into one csv file.

        :param path: The path to the save place
        :type path: str
        """
        combined = self.__combine_all__()
        combined.to_csv(path)

        if self._spectra is not None:
            self._spectra.write(f'{os.path.dirname(path)}/spectra/')

    def write(self, path, format='zip', in_zip_format='fits'):
        """
        Writes the current data to a file (-system).

        For example:

        .. code:: python

            ds.write('my_data.zip')

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
        # elif format == 'report':
        #     self.__write_as_report__(path)
        else:
            raise ValueError(f'Format {format} is not supported.')

    def __str__(self):
        s = f'Number of entries: {len(self.coordinates)}\n'
        s += 'Available surveys:\n'
        if self._magnitudes is not None:
            if self._magnitudes.survey is not None:
                for su in self._magnitudes.survey.get_surveys():
                    s += f'\t{su}\n'
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

    def __sizeof__(self):
        total_size = 56
        if self._coordinates is not None:
            total_size += sys.getsizeof(self._coordinates)
        return total_size

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
        if not os.path.exists(name):
            raise FileNotFoundError(f'FileNotFoundError: [Errno 2] No such file or directory: {name}')
        if format == 'zip':
            d = {}
            head = None
            with zipfile.ZipFile(name) as zip_file:
                mags = None
                flux = FluxTable()
                colors = Colors()
                for file_name in zip_file.namelist():

                    if '.ini' not in file_name:
                        if 'magnitude' in file_name:
                            if mags is None:
                                mags = Magnitude(read_from_zip(zip_file, file_name),
                                                 survey=file_name.split('.')[0].split('magnitudes_')[-1])
                            else:
                                mags.add_survey_mags(read_from_zip(zip_file, file_name),
                                                     file_name.split('.')[0].split('magnitudes_')[-1])
                        elif 'flux' in file_name:
                            flux.add_fluxes(read_from_zip(zip_file, file_name),
                                            file_name.split('.')[0].split('magnitudes_')[-1])
                        elif 'colors' in file_name:
                            colors.add_colors(read_from_zip(zip_file, file_name),
                                              file_name.split('.')[0].split('magnitudes_')[-1])
                        elif 'coordinates' in file_name:
                            d['coordinates'] = read_from_zip(zip_file, file_name)
                        elif 'light_curves' in file_name:
                            # todo: implement light curve load
                            pass
                    else:
                        head = read_from_zip(zip_file, file_name)
                d['magnitudes'] = mags
                d['flux'] = flux
                d['colors'] = colors
            ds = DataSet(**d)
            ds._magnitudes._survey = head
            ds.reset_masks()
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
            raise NotImplementedError('Automatic csv reading not implemented')
            # todo: implement csv reading
            pass
        else:
            raise ValueError(f'Format {format} is not supported.')

    @staticmethod
    def load_coordinates(path, format='fits', ra_name='ra', dec_name='dec'):
        """
        Creates a DataSet object from a file with coordinates.

        .. code-block::

            from Phosphorpy import DataSet

            ds = DataSet.load_coordinates('./my_coordinates.fits', ra_name='RAJ2000', dec_name='DEJ2000')

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
            coords = coords.rename({ra_name: 'ra', dec_name: 'dec'}, axis='columns')
        elif format == 'csv':
            coords = pd.read_csv(path)[[ra_name, dec_name]]
            coords = coords.rename({ra_name: 'ra', dec_name: 'dec'}, axis='columns')
        else:
            raise ValueError(f'Format \'{format}\' is not supported.')

        return DataSet(coords)

    @staticmethod
    def from_coordinates(data=None, ra=None, dec=None):
        """
        Creates a DataSet object from a set of coordinates.
        The coordinates can be an astropy Table, astropy SkyCoord or a pandas DataFrame.
        Also ra and dec can be a list/array of coordinates.

        At least data or ra and dec must be given.

        :param data: The coordinates for the DataSet
        :type data: Table, SkyCoord, DataFrame
        :param ra: The RA-components of the coordinates or a string with the name of the RA column
        :type ra: None, ndarray, str
        :param dec: The Dec-components of the coordinates or a string with the name of the Dec column
        :type dec: None, ndarray, str
        :return: The DataSet object with the given coordinates
        :rtype: DataSet
        """
        if data is None:
            if ra is None or dec is None:
                raise ValueError('Data or ra and dec must be not None!')

            if len(ra) != len(dec):
                raise ValueError('ra and dec must have the same length.')

            data = DataFrame({'ra': ra, 'dec': dec})
        else:
            if type(data) == Table:
                data = data.to_pandas()
            elif type(data) == SkyCoord:
                data = DataFrame({'ra': data.ra.degree, 'dec': data.dec.degree})
            elif type(data) != DataFrame:
                raise NotImplementedError('Only astropy Table, astropy SkyCoord and pandas DataFrame are implement.')

            if 'ra' in data.columns and 'dec' in data.columns:
                data = data[['ra', 'dec']]
            elif ra is None or dec is None:
                raise ValueError('If the input data does not contain \'ra\' and \'dec\' ra and dec must'
                                 ' be the names of the corresponding columns.')
            else:
                data = data[[ra, dec]]
                data = data.rename({ra: 'ra', dec: 'dec'}, axis='columns')

        return DataSet(data)

    @staticmethod
    def from_vizier(name, **constrains):
        """
        Stars a constrain query on a vizier catalog and convert the output to a DataSet object.

        .. code-block::

            from Phosphorpy import DataSet

            ds = DataSet.from_vizier('SDSS', {'rmag': '<17', 'class': 'STAR'})

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
