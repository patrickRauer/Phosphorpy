from astropy.table import Table
from astropy.modeling import models, fitting
from astropy import units as u
import numpy as np
import glob
import os
import numbers

from Phosphorpy.data.sub.plots.spectra import SpectraPlot, SpectraListPlot

try:
    from Phosphorpy.data.sub.interactive_plotting.spectra import SpectraPlot as SpectraPlotHV
    from Phosphorpy.data.sub.interactive_plotting.spectra import SpectraListPlot as SpectraListPlotHV
except ImportError:
    SpectraPlotHV = None
    SpectraListPlotHV = None


class SpectraList:

    _spectra = None
    _ids = None
    _plot = None

    def __init__(self, spectra=None, ids=None):
        """
        SpectraList provides the handling of multiple spectra

        :param spectra: Optional: A spectra or a list of spectra
        :type spectra: None, Spectra, list
        :param ids:
            Optional: ID's of the spectra. If they are given, they must have the same length as the given spectra.
        :type ids: int, list
        """
        if spectra is None:
            self._spectra = []
            self._ids = []
        elif type(spectra) == list or type(spectra) == tuple:
            self._spectra = list(spectra)
            if ids is None:
                self._ids = list(np.arange(len(spectra)))
            else:
                if len(spectra) == len(ids):
                    self._ids = list(ids)
                else:
                    raise ValueError('If spectra and ids are given, they must have the same length.')
        else:
            self._spectra = [spectra]
            if ids is not None:
                self._ids = [ids]
            else:
                self._ids = [0]

        self._plot = SpectraListPlot(self)

        if SpectraListPlotHV.holoviews():
            self._hv_plot = SpectraListPlotHV(self)

    def __len__(self):
        return len(self._spectra)

    def __str__(self):
        return f'SpectraList with {len(self)} spectra.'

    def __repr__(self):
        return str(self)

    def __getitem__(self, item):
        """

        :param item:
        :return:
        :rtype: Spectra, int
        """
        if isinstance(item, numbers.Integral):
            return self._spectra[item], self._ids[item]
        elif type(item) == list or type(item) == tuple or type(item) == np.ndarray:
            spec_out = [
                self._spectra[i] for i in item
            ]
            ids_out = [
                self._ids[i] for i in item
            ]
            return SpectraList(spec_out, ids_out)

    def write(self, path, data_format='fits', overwrite=True):
        """
        Writes the spectra of this SpectraList to files

        :param path: The path to the directory, where the files should be placed.
        :type path: str
        :param data_format: The format of the spectra files
        :type data_format: str
        :param overwrite: True, if an existing file should be overwritten, else False. Default is True
        :type overwrite: bool
        :return:
        """
        if not os.path.exists(path):
            os.makedirs(path)
        elif not os.path.isdir(path):
            raise ValueError('path must be a directory and not a file.')

        paths = []
        for spec, index in zip(self._spectra, self._ids):
            paths.append(spec.write(
                f'{path}{spec.survey}_{index}.{data_format}',
                data_format=data_format, overwrite=overwrite
            ))
        return paths

    @staticmethod
    def read(path, format='fits', wavelength_name='wavelength', flux_name='flux', survey_key='survey'):
        """
        Reads all spectra from the directory

        :param path: Path to the directory with the spectra
        :type path: str
        :param format: The format of the files of the spectra.
        :type format: str
        :param wavelength_name: The name of the wavelength column
        :type wavelength_name: str
        :param flux_name: The name of the flux column
        :type flux_name: str
        :param survey_key: The key of the survey information
        :type survey_key: str
        :return: The SpectraList with all the red spectra
        :rtype: SpectraList
        """
        spec_list = SpectraList()
        for f in glob.glob(f'{path}*.{format}'):
            spec = Spectra.read(f, format=format,
                                wavelength_name=wavelength_name, flux_name=flux_name,
                                survey_key=survey_key)
            spec_list.append(spec)
        return spec_list

    def get_by_id(self, index):
        """
        Returns the spectra with the required index

        :param index: The index of the required spectra
        :type index: int
        :return:
            SpectraList with the required spectra. The list can contain more than one spectra, if multiple spectra
            with the required ID are found.
            If no spectra is found with the required ID, the method will return None.
        :rtype: SpectraList, None
        """
        if isinstance(index, numbers.Integral):
            con = index in self._ids
        else:
            con = np.isin(index, np.array(self._ids)).any()

        if con:
            p = np.where(np.array(self._ids) == index)[0]
            spec_list = SpectraList()
            for i in p:
                spec_list.append(*(self[i]))
            return spec_list

    def get_ids(self):
        """
        Returns all IDs of the downloaded spectra. Every ID is returned only once but multiple spectra could
        be available.
        :return: The IDs of the spectra
        :rtype: ndarray
        """
        return np.unique(self._ids)

    def append(self, spectra, spec_id=-1):
        """
        Appends a new spectra to the spectra list.

        :param spectra: The new spectra
        :type spectra: LamostSpectra
        :param spec_id: The ID of the spectra
        :type spec_id: int
        :return:
        """
        self._spectra.append(spectra)
        if spec_id > -1:
            self._ids.append(spec_id)
        else:
            self._ids.append(len(self._ids))

        if spectra.index == -1:
            spectra.index = self._ids[-1]

    def merge(self, second):
        """
        Merges a second SpectraList into this one

        :param second: The other SpectraList
        :type second: SpectraList
        :return:
        """
        if second is None:
            return

        if type(second) != SpectraList:
            raise ValueError('second must be of the type \'SpectraList\'')

        for i in range(len(second)):
            temp = second[i]
            self.append(temp[0].copy(), temp[1])

    # def estimate_line_properties(self, as_velocity=False, redo=False):
    #     """
    #     Estimates the line properties of all spectra in the list
    #
    #     :param as_velocity:
    #         True, if the line shift should be returned as a radial velocity in km/s, else False to
    #         get lambda/lambda0-1.
    #         Default is False.
    #     :type as_velocity: bool
    #     :param redo:
    #         True, if old results should be ignored, else False.
    #         Default is False.
    #     :type redo: bool
    #     :return: The line properties of all spectra
    #     :rtype: dict
    #     """
    #     out = {}
    #     for s in self._spectra:
    #         out[s.obs_id] = s.estimate_line_properties(as_velocity=as_velocity,
    #                                                    redo=redo)
    #     return out
    #
    # def as_dataframe(self, as_velocity=False, redo=False):
    #     """
    #     Returns the main information of the stored spectra as a pandas DataFrame
    #
    #     :param as_velocity:
    #     :param redo:
    #     :return: The main information as a dataframe
    #     :rtype: DataFrame
    #     """
    #     out = []
    #     for s, d_id in zip(self._spectra, self._ids):
    #         properties = s.estimate_line_properties(as_velocity=as_velocity,
    #                                                 redo=redo)
    #         properties['obsID'] = s.obs_id
    #         properties['ID'] = d_id
    #         out.append(properties)
    #     out = vstack(out).to_pandas()
    #     return out.set_index('ID')

    @property
    def plot(self):
        return self._plot

    @property
    def hvplot(self):
        return self._hv_plot


class Spectra:
    NORMALIZE_MEAN = np.mean
    """
    Static variable for the normalization with the mean value
    """
    NORMALIZE_MEDIAN = np.median
    """
    Static variable for the normalization with the median value
    """
    NORMALIZE_MAX = np.max
    """
    Static variable for the normalization with the max value
    """
    NORMALIZE_SUM = np.sum
    """
    Static variable for the normalization with the sum
    """

    _wavelength = None
    _flux = None

    _wavelength_unit = None
    _flux_unit = None

    _survey = None
    _index = -1
    _meta = None

    _fits = None

    _lines = None

    _plot = None
    _hv_plot = None

    def __init__(self, wavelength=None, flux=None, wavelength_unit=None, flux_unit=None, survey=None,
                 index=-1, meta=None):
        """
        Spectra is the basic class to handle different kinds of spectra on a basic level
        without any specific functionality related to any certain spectra

        :param wavelength:
            The wavelength values of the spectra or None.
            As wavelength values a Union (tuple, list, array) is allowed or a Quantity.
        :type wavelength: Union, Quantity
        :param flux:
            The flux values of the spectra or None.
            As flux values a Union (tuple, list, array) is allowed or a Quantity.
        :type flux: Union, Quantity
        :param wavelength_unit:
            The units of the wavelengths, if the wavelengths are given. If no unit is given, angstrom are
            assumed to be the wavelength unit.
        :type wavelength_unit: Unit
        :param flux_unit:
            The units of the flux, if the flux is given. If no unit is given, ergs are assumed to be
            the flux unit.
        :type flux_unit: Unit
        :param survey: Name of the survey. Default is None.
        :type survey: str, None
        :param index: Index of the spectra. Default is -1.
        :type index: int
        """
        if type(wavelength) == u.Quantity:
            wavelength_unit = wavelength.unit
            wavelength = wavelength.value

        if type(flux) == u.Quantity:
            flux_unit = flux.unit
            flux = flux.value

        self._wavelength = wavelength
        self._flux = flux

        self.wavelength_unit = wavelength_unit
        self.flux_unit = flux_unit

        self._survey = survey
        self._index = index
        self._meta = meta

        self._plot = SpectraPlot(self)
        self._fits = []

        if SpectraPlotHV.holoviews():
            self._hv_plot = SpectraPlotHV(self)

    def __str__(self):
        return f'Spectra with wavelength between {self.min_wavelength} and {self.max_wavelength} with a' \
               f' wavelength resolution of {np.round(self.resolution_wavelength, 2)}.'

    def __repr__(self):
        return str(self)

    def write(self, path, data_format='fits', overwrite=True):
        """
        Writes the spectral data to a file

        :param path: The path to the file.
        :type path: str
        :param data_format: The format of the output file. Current supported formats are fits and csv.
        :type data_format: str
        :param overwrite:
            True if an existing file should be overwritten, else False.
            Default is False.
        :return:
        """
        tab = Table()
        tab['wavelength'] = self.wavelength
        tab['wavelength'].unit = self.wavelength_unit
        tab['flux'] = self.flux
        tab['flux'].unit = self.flux_unit
        tab.meta['survey'] = self.survey
        tab.meta['index'] = self._index

        if data_format == 'fits':
            data_format = 'fits'
        elif data_format == 'csv':
            data_format = 'ascii.csv'
        else:
            raise ValueError(f'Format: {format} is not supported.\nChoose \'fits\' or \'csv\'.')

        tab.write(path, format=data_format, overwrite=overwrite)
        return path

    @staticmethod
    def read(path, format='fits', wavelength_name='wavelength', flux_name='flux', survey_key='survey'):
        """
        Reads a spectra from a file

        :param path: Path to the file.
        :type path: str
        :param format: Format of the file. Current supported formats are fits and csv
        :type format: str
        :param wavelength_name: The name of the wavelength column
        :type wavelength_name: str
        :param flux_name: The name of the flux column
        :type flux_name: str
        :param survey_key: The name of the header entry with the name of the survey/origin
        :type survey_key: str
        :return: The Spectra-object with the data from the file
        :rtype: Spectra
        """
        if format == 'fits':
            tab = Table.read(path, format='fits')
        elif format == 'csv':
            tab = Table.read(path, format='ascii.csv')
        else:
            raise ValueError(f'Format: {format} is not supported.\nChoose \'fits\' or \'csv\'')

        try:
            survey = tab.meta[survey_key]
        except KeyError:
            survey = None

        try:
            index = tab.meta['index']
        except KeyError:
            index = -1
        return Spectra(
            wavelength=tab[wavelength_name],
            flux=tab[flux_name],
            wavelength_unit=tab[wavelength_name].unit,
            flux_unit=tab[flux_name].unit,
            survey=survey, index=index
        )

    def copy(self):
        """
        Make a copy of the Spectra class

        :return: A copy of the current Spectra-object
        :rtype: Spectra
        """
        return Spectra(
            wavelength=self.wavelength.copy(),
            flux=self.flux.copy(),
            wavelength_unit=self.wavelength_unit,
            flux_unit=self.flux_unit,
            survey=self.survey, index=self.index
        )

    def normalize(self, kind, inplace=True):
        """
        Normalize the spectra with a certain function.

        .. math::

            flux^* = \frac{flux}{function(flux)}


        :param kind:
            The kind of function, which is used to normalize the spectra.
            Allowed function are the static function or 'mean', 'median', 'max' or 'sum' as ky words.
        :type kind: function, str
        :param inplace:
            True, if the flux should be overwritten with the normalized flux, else False for a new spectra
            object. Default is True.
        :type inplace: bool
        :return: None or if inplace is False, a new spectra object.
        :rtype: None, Spectra
        """
        if type(kind) == str:
            kind = kind.lower()
        fl = self.flux.copy()
        if kind == self.NORMALIZE_MAX or kind == 'max':
            norm = np.max(fl)
        elif kind == self.NORMALIZE_MEAN or kind == 'mean':
            norm = np.mean(fl)
        elif kind == self.NORMALIZE_MEDIAN or kind == 'median':
            norm = np.median(fl)
        elif kind == self.NORMALIZE_SUM or kind == 'sum':
            norm = np.sum(fl)
        else:
            raise ValueError('Kind of normalization is not allowed.')

        if norm == 0:
            raise ValueError('Normalization not possible. Norm is 0.')

        fl /= norm

        if type(inplace) == bool:
            if inplace:
                self._flux = fl
                self._flux_unit /= self._flux_unit
            else:
                return Spectra(wavelength=self.wavelength.copy(),
                               flux=fl)
        else:
            raise TypeError('Inplace must be a bool-type.')

    def cut(self, min_wavelength=None, max_wavelength=None, inplace=True):
        """
        Make a cutout of the spectra between the minimal and maximal wavelengths.
        If nether minimal wavelength or maximal wavelength is given,
        it will make a copy of the original data.

        :param min_wavelength:
            The minimal wavelength as a float, then the we assume the same units as the wavelength,
            a Quantity with the dimension length or None.
        :type min_wavelength: float, Quantity, None
        :param max_wavelength:
            The maximal wavelength as a float, then the we assume the same units as the wavelength,
            a Quantity with the dimension length or None.
        :type max_wavelength: float, Quantity, None
        :param inplace:
            True, if the current spectra should be overwritten with the cutout, else False to get
            a new Spectra object. Default is True.
        :type inplace: bool
        :return:
        """

        wave = self.wavelength.copy()
        flux = self.flux.copy()

        if min_wavelength is not None:
            # if the minimal wavelength is a astropy Quantity, align the units and apply the limit then
            if type(min_wavelength) == u.Quantity:
                min_wavelength = min_wavelength.to(self.wavelength_unit).value
            m = wave >= min_wavelength
            flux = flux[m]
            wave = wave[m]

        if max_wavelength is not None:
            # if the maximal wavelength is a astropy Quantity, align the units and apply the limit then
            if type(max_wavelength) == u.Quantity:
                max_wavelength = max_wavelength.to(self.wavelength_unit).value
            m = wave <= max_wavelength
            flux = flux[m]
            wave = wave[m]

        if type(inplace) == bool:
            if inplace:
                self._wavelength = wave
                self._flux = flux
            else:
                return Spectra(wavelength=wave, flux=flux,
                               wavelength_unit=self.wavelength_unit,
                               flux_unit=self.flux_unit)
        else:
            raise ValueError('inplace must be a bool.')

    def fit_line(self, model=None):
        """
        Fit's the spectra/line with the model. The model must be a fittable object from
        the astropy.modeling.model classes.
        If no model is given, a 1D gaussian is used.

        :param model: The fittable model or None. Default is None, which means that a 1D gaussian is used.
        :return: The resulting fit
        """
        if model is None:
            model = models.Gaussian1D(stddev=10)
        else:
            try:
                model.fittable
            except AttributeError:
                raise ValueError(f'Model must be an instance of astropy\'s Fittable1DModel and not {type(model)}')

        fitter = fitting.SLSQPLSQFitter()
        fit_rs = fitter(model, self.wavelength, self.flux)
        self._fits.append(fit_rs)
        return fit_rs

    def fit_gauss(self, guesses):
        """
        Fit's the current spectra/line with a 1D gaussian

        :param guesses: Initial guesses for the gaussian
        :type guesses: dict
        :return: The resulting fit
        """
        gauss = models.Gaussian1D(**guesses)
        return self.fit_line(model=gauss)

    def fit_double_gauss(self, guesses, guesses2):
        """
        Fit's the current spectra/line(s) with two 1D gaussian's

        :param guesses: Initial guesses for the first gaussian
        :type guesses: dict
        :param guesses2: Initial guesses for the second gaussian
        :type guesses2: dict
        :return: The resulting fit
        """
        dgauss = models.Gaussian1D(**guesses)+models.Gaussian1D(**guesses2)
        return self.fit_line(model=dgauss)

    def estimate_line_properties(self):
        pass

    @property
    def wavelength(self):
        """
        The wavelengths of spectra
        :return:
        """
        return self._wavelength.copy()

    @property
    def wavelength_unit(self):
        """
        The wavelength unit
        :return:
        """
        return self._wavelength_unit

    @wavelength_unit.setter
    def wavelength_unit(self, unit):
        """
        Sets a new wavelength unit
        :param unit: The unit of the dimension length
        :type unit: Unit
        :return:
        """
        # if no unit is given, assume angstrom
        if unit is None:
            # if no previous unit was set
            if self._wavelength_unit is None:
                self._wavelength_unit = u.angstrom
            # if wavelength has a unit, recursive call to change the unit to angstrom
            else:
                self.wavelength_unit = u.angstrom
        elif not isinstance(unit, u.Unit):
            raise ValueError('The new unit must be a astropy unit.')
        elif self._wavelength_unit is None:
            self._wavelength_unit = unit
        else:
            try:
                # estimate the conversion factor and apply it to the wavelengths
                factor = (self._wavelength_unit.to(unit))
                self._wavelength *= factor
                self._wavelength_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a length.')

    @property
    def flux(self):
        """
        The flux values
        :return:
        """
        return self._flux

    @property
    def flux_unit(self):
        """
        The flux unit
        :return:
        """
        return self._flux_unit

    @flux_unit.setter
    def flux_unit(self, unit):
        """
        Sets a new flux unit
        :param unit: The unit of the dimension of a flux
        :type unit: Unit
        :return:
        """
        if unit is None:
            self._flux_unit = u.erg
        elif not isinstance(unit, u.Unit):
            raise ValueError('The new unit must be a astropy unit.')
        elif self._flux_unit is None:
            self._flux_unit = unit
        else:
            try:
                # estimate the conversion factor and apply it to the fluxes
                factor = (self._flux_unit.to(unit))
                self._flux *= factor
                self._flux_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a flux.')

    @property
    def plot(self):
        """
        The plotting environment
        :return:
        """
        return self._plot

    @property
    def hvplot(self):
        return self._hv_plot

    @property
    def min_wavelength(self):
        return self.wavelength.min()

    @property
    def max_wavelength(self):
        return self.wavelength.max()

    @property
    def resolution_wavelength(self):
        return np.mean(self._wavelength[1:]-self._wavelength[:-1])

    @property
    def min_flux(self):
        return self.flux.min()

    @property
    def max_flux(self):
        return self.flux.max()

    @property
    def fit(self):
        return self._fits

    @property
    def survey(self):
        return self._survey

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, value):
        self._index = value
