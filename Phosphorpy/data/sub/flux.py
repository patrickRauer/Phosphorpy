from .table import DataTable
from Phosphorpy.data.sub.plots.flux import FluxPlot
from Phosphorpy.data.sub.tables.flux import Flux
from astropy import units as u
from astropy.units.quantity import Quantity
from pandas import DataFrame
from scipy.optimize import curve_fit
import numpy as np
import inspect


def get_column_names(length):
    # create columns names 'a', 'b', 'c', ...
    cols = []
    for i in range(97, 98 + length):
        cols.append(chr(i))
    return cols


def _apply_wavelength_limits(wavelengths, d, err, lower_limit, upper_limit):
    """
    Applies wavelength limits to the data

    :param wavelengths: The wavelengths of the d columns
    :type wavelengths: Union
    :param d:
        The data, rows are the different samples and columns are the data-points
        for every sample
    :type d: Union
    :param err:
        The errors of the data. The shape must be the same as d's shape.
    :type err: Union
    :param lower_limit: The lower wavelength limit
    :type lower_limit: float, Quantity
    :param upper_limit: The upper wavelength limit
    :type upper_limit: float, Quantity
    :return: The wavelengths, data and errors in the wavelength limit
    :rtype: Union, Union, Union
    """
    if lower_limit is not None:
        if type(lower_limit) != Quantity:
            lower_limit = lower_limit*u.angstrom

        wavelength_mask = (wavelengths*u.angstrom > lower_limit)
        d = d[:, wavelength_mask]
        err = err[:, wavelength_mask]
        wavelengths = wavelengths[wavelength_mask]

    if upper_limit is not None:
        if type(upper_limit) != Quantity:
            upper_limit = upper_limit*u.angstrom
        wavelength_mask = (wavelengths*u.angstrom < upper_limit)
        d = d[:, wavelength_mask]
        err = err[:, wavelength_mask]
        wavelengths = wavelengths[wavelength_mask]
    return wavelengths, d, err


def _error_weighted_fitting(wavelengths, d, err, degree):
    """
    Fit all data with a polynomial. if the data or the errors contain NaN values,
    then they will be excluded from the fitting process.

    If the fit fails, all values are filled with -99999.

    :param wavelengths: The wavelengths of the d columns
    :type wavelengths: Union
    :param d:
        The data, rows are the different samples and columns are the data-points
        for every sample
    :type d: Union
    :param err:
        The errors of the data. The shape must be the same as d's shape.
    :type err: Union
    :param degree: The degree of the polynomial
    :type degree: int
    :return: The coefficients from the fit's.
    :rtype: list
    """
    fits = []
    for row, err_row in zip(d, err):
        try:
            m = (row > -9999) & (err_row > -999)
            fit = np.polyfit(wavelengths[m], row[m], degree, w=1/err_row[m])
            fits.append(fit)
        except np.linalg.LinAlgError:
            fits.append((degree+1)*[-99999])
        except TypeError:
            fits.append((degree+1)*[-99999])
    return fits


def _fitting(wavelengths, d, degree):
    """
    Fit all data with a polynomial. if the data contain NaN values, then they will
    be excluded from the fitting process.

    If the fit fails, all values are filled with -99999.

    :param wavelengths: The wavelengths of the d columns
    :type wavelengths: Union
    :param d:
        The data, rows are the different samples and columns are the data-points
        for every sample
    :type d: Union
    :param degree: The degree of the polynomial
    :type degree: int
    :return: The coefficients from the fit's.
    :rtype: list
    """
    fits = []
    for row in d:
        try:
            m = (row > -9999)
            fit = np.polyfit(wavelengths, row[m], degree)
            fits.append(fit)
        except np.linalg.LinAlgError:
            fits.append((degree+1)*[-99999])
        except TypeError:
            fits.append((degree+1)*[-99999])
    return fits


def _normalize(d, err, norm):
    """
    If norm is true, normalize the data by their maximum value otherwise do
    nothing.

    :param d:
        The data, rows are the different samples and columns are the data-points
        for every sample
    :type d: Union
    :param err:
        The errors of the data. The shape must be the same as d's shape.
    :type err: Union
    :param norm:
        True if the data should be normalized by their maximum value else False.
        Default is False
    :return:
        If norm is True, the normalized data, if norm is False, the input data
    :rtype: Union, Union
    """
    if norm:
        maxi = np.nanmax(d, axis=1)
        maxi = np.transpose(np.tile(maxi, (d.shape[1], 1)))
        d = d/maxi
        err = err/maxi
    return d, err


def _rescale(wavelengths, d, err, x_log, y_log):
    """
    Rescales the wavelengths and/or the flux to a logarithmic scale.

    :param wavelengths: The wavelengths of the d columns
    :type wavelengths: Union
    :param d:
        The data, rows are the different samples and columns are the data-points
        for every sample
    :type d: Union
    :param err:
        The errors of the data. The shape must be the same as d's shape.
    :param x_log: If True the wavelengths will become logarithmic.
    :type x_log: bool
    :param y_log: If True, the flux and its errors will become logarithmic.
    :type y_log: bool
    :return:
        If x_log and/or y_log are true, the logarithmic values, else
        the input values
    :rtype: Union, Union, Union
    """
    if x_log:
        wavelengths = np.log10(wavelengths)

    if y_log:
        err = np.abs(1 / d) * err
        d = np.log10(d)

    return wavelengths, d, err


class FluxTable(DataTable):
    """
    FluxTable is the interface to interact with the fluxes of the sources.
    """

    _fits = None
    _data = None
    _survey = None

    def __init__(self, data=None, survey_head=None, mask=None):
        """

        :param data:
        :param survey_head: Survey information
        :type survey_head: SearchEngine.data.sub.magnitudes.SurveyData
        """
        DataTable.__init__(self, mask=mask)
        if data is None:
            data = []
        self._data = data
        self._survey = survey_head
        self._plot = FluxPlot(self)

    def add_fluxes(self, data, survey_name):
        fl = Flux(data, survey=survey_name)
        self.data.append(fl)

    def get_all_fluxes(self):
        flux = None
        for d in self.data:
            if flux is None:
                flux = d.get_fluxes()
            else:
                flux = flux.merge(d.get_fluxes(), how='outer', right_index=True, left_index=True)
        return flux.values

    def get_all_errors(self):
        errs = None
        for d in self.data:
            if errs is None:
                errs = d.get_errors()
            else:
                errs = errs.merge(d.get_errors(), how='outer', right_index=True, left_index=True)
        return errs.values

    def _get_sed_data(self):
        """
        Returns the data of the SED's in structured way.
        First wavelength, then flux and as the last component the errors

        :return: wavelength, flux, flux error
        :rtype: Union, Union, Union
        """
        d = self.get_all_fluxes()
        err = self.get_all_errors()

        wavelengths = self._survey.get_all_wavelengths()

        return wavelengths, d, err

    def get_sed(self, index):
        fluxes = []
        errs = []
        for d in self.data:
            survey_flux, survey_err = d.get_sed(index)
            fluxes.extend(survey_flux)
            errs.extend(survey_err)
        return np.array(fluxes), np.array(errs)

    def fit_polynomial(self, degree, error_weighted=True,
                       x_log=False, y_log=False,
                       norm=False, lower_limit=None, upper_limit=None):
        """
        Fit's all flux SED's with polynomial

        If the fit fails, the values of the coefficients are set to -99999.

        :param degree: The degree of the polynomial
        :type degree: int
        :param error_weighted:
            True if the error should be used as weights, else False.
            Default is True.
        :type error_weighted: bool
        :param x_log:
            True if the wavelengths should be logarithmic, else False.
            Default is False.
        :type x_log: bool
        :param y_log:
            The if the flux should be logarithmic, else False.
            Default is False.
        :type y_log: bool
        :param norm:
            The if the SED should be normalized by its maximum value, else False.
            Default is False.
        :type norm: bool
        :param lower_limit:
            The lower wavelengths limit or None if no lower limit should be applied.
            Default is None.
        :type lower_limit: None, float, Quantity
        :param upper_limit:
            The upper wavelengths limit or None if no upper limit should be applied.
            Default is None.
        :type upper_limit: None, float, Quantity
        :return: The results of the fit's
        :rtype: DataFrame
        """
        wavelengths, d, err = self._get_sed_data()

        # apply wavelength cuts
        wavelengths, d, err = _apply_wavelength_limits(wavelengths, d, err,
                                                       lower_limit, upper_limit)

        wavelengths, d, err = _rescale(wavelengths, d, err, x_log, y_log)

        d, err = _normalize(d, err, norm)

        if error_weighted:
            fits = _error_weighted_fitting(wavelengths, d, err, degree)
        else:
            fits = _fitting(wavelengths, d, degree)

        self._fits = DataFrame(data=np.array(fits),
                               columns=get_column_names(degree))

        return self._fits

    def func_fit(self, func, error_weighted=True,
                 x_log=False, y_log=False,
                 norm=False, lower_limit=None, upper_limit=None,
                 **curve_fit_parameters):
        """
        Fit's the wavelength and fluxes with the given function and stores the results in fits

        :param func: The fitting function in the style of scipy.optimize.curve_fit
        :param error_weighted:
            True if the error should be used as weights, else False.
            Default is True.
        :type error_weighted: bool
        :param x_log:
            True if the wavelengths should be logarithmic, else False.
            Default is False.
        :type x_log: bool
        :param y_log:
            The if the flux should be logarithmic, else False.
            Default is False.
        :type y_log: bool
        :param norm:
            The if the SED should be normalized by its maximum value, else False.
            Default is False.
        :type norm: bool
        :param lower_limit:
            The lower wavelengths limit or None if no lower limit should be applied.
            Default is None.
        :type lower_limit: None, float, Quantity
        :param upper_limit:
            The upper wavelengths limit or None if no upper limit should be applied.
            Default is None.
        :type upper_limit: None, float, Quantity
        :param curve_fit_parameters: Additional parameters for the curve_fit function
        :return: The results of the fit's
        :rtype: DataFrame
        """
        wavelengths, d, err = self._get_sed_data()

        # apply wavelength cuts
        wavelengths, d, err = _apply_wavelength_limits(wavelengths, d, err,
                                                       lower_limit, upper_limit)

        wavelengths, d, err = _rescale(wavelengths, d, err, x_log, y_log)

        d, err = _normalize(d, err, norm)

        colnames = inspect.getfullargspec(func)[0][1:]

        if error_weighted:
            fits = []
            for row in d:
                try:
                    m = (row > -9999)
                    fit, cov = curve_fit(func, wavelengths, row[m], sigma=err,
                                         **curve_fit_parameters)
                    fits.append(fit)
                except RuntimeError:
                    fits.append(len(colnames)*[-99999])
                except TypeError:
                    fits.append(len(colnames)*[-99999])
        else:
            fits = []
            for row in d:
                try:
                    m = (row > -9999)
                    fit, cov = curve_fit(func, wavelengths, row[m],
                                         **curve_fit_parameters)
                    fits.append(fit)
                except RuntimeError:
                    fits.append(len(colnames)*[-99999])
                except TypeError:
                    fits.append(len(colnames)*[-99999])

        self._fits = DataFrame(data=np.array(fits),
                               columns=colnames)

        return self._fits

    def create_decomposition(self, n_components):
        if self._fits is None:
            raise AttributeError('No fitting date found. Call first a fitting method.')
        
    def decomposition(self, n_components, model):
        if self._fits is None:
            raise AttributeError('No fitting date found. Call first a fitting method.')

    def fit_blackbody(self):
        pass

    def __str__(self):
        return str(self._data)

    @property
    def fit(self):
        return self._fits
