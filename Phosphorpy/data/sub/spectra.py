from astropy.table import vstack
from astropy.modeling import models, fitting, Fittable1DModel
from astropy import units as u
import numpy as np

from Phosphorpy.data.sub.plots.spectra import SpectraPlot


class SpectraList:

    _spectra = None
    _ids = None

    def __init__(self, spectra=None):
        if spectra is None:
            self._spectra = []
            self._ids = []
        elif type(spectra) == list or type(spectra) == tuple:
            self._spectra = list(spectra)
        else:
            self._spectra = [spectra]
            self._ids = [0]

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

    def estimate_line_properties(self, as_velocity=False, redo=False):
        """
        Estimates the line properties of all spectra in the list

        :param as_velocity:
            True, if the line shift should be returned as a radial velocity in km/s, else False to
            get 1-lambda/lambda0.
            Default is False.
        :type as_velocity: bool
        :param redo:
            True, if old results should be ignored, else False.
            Default is False.
        :type redo: bool
        :return: The line properties of all spectra
        :rtype: dict
        """
        out = {}
        for s in self._spectra:
            out[s.obs_id] = s.estimate_line_properties(as_velocity=as_velocity,
                                                       redo=redo)
        return out

    def as_dataframe(self, as_velocity=False, redo=False):
        """
        Returns the main information of the stored spectra as a pandas DataFrame

        :param as_velocity:
        :param redo:
        :return: The main information as a dataframe
        :rtype: DataFrame
        """
        out = []
        for s, d_id in zip(self._spectra, self._ids):
            properties = s.estimate_line_properties(as_velocity=as_velocity,
                                                    redo=redo)
            properties['obsID'] = s.obs_id
            properties['ID'] = d_id
            out.append(properties)
        out = vstack(out).to_pandas()
        return out.set_index('ID')

    def __len__(self):
        return len(self._spectra)


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

    _fits = None

    _lines = None

    _plot = None

    def __init__(self, wavelength=None, flux=None, wavelength_unit=None, flux_unit=None):
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

        self._plot = SpectraPlot(self)
        self._fits = []

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
        # if both are None, raise an error because it does not make any sense
        if min_wavelength is max_wavelength:
            raise ValueError('At least one of minimal wavelength or maximal wavelength must be given.')

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
        gauss = models.Gaussian1D(**guesses)
        return self.fit_line(model=gauss)

    def fit_double_gauss(self, guesses, guesses2):
        dgauss = models.Gaussian1D(**guesses)+models.Gaussian1D(**guesses2)
        return self.fit_line(model=dgauss)

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
    def min_wavelength(self):
        return self.wavelength.min()

    @property
    def max_wavelength(self):
        return self.wavelength.max()

    @property
    def min_flux(self):
        return self.flux.min()

    @property
    def max_flux(self):
        return self.flux.max()

    @property
    def fit(self):
        return self._fits
