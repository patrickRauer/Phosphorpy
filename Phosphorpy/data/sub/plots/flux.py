from scipy.optimize import curve_fit
import numpy as np
import pylab as pl


class FluxPlot:
    _data = None

    def __init__(self, data):
        """

        :param data:
        :type data: Phosphorpy.data.sub.flux.FluxTable
        """
        self._data = data

    def sed(self, index, fit=None, path='', x_log=False, y_log=False, legend=False):
        """
        Plot the SED of the source with the given index.

        :param index: The index of the source
        :type index: int
        :param fit:
            A function to fit, a string with the name of a predefined function or an integer for a polynomial fit.
        :type fit: function, str, int
        :param path:
            Path to the place where the figure should be saved. Default is an empty string, which means that
            the figure will be shown only.
        :type path: str
        :param x_log: True if the x-axis should be in log scale, else False. Default is False.
        :type x_log: bool
        :param y_log: True if the y-axis should be in log-scale, else False. Default is False
        :type y_log: bool
        :param legend: True if the legend should be shown, else False. Default is False.
        :type legend: bool
        :return:
        """
        print(index, fit, x_log, y_log)
        pl.clf()
        sp = pl.subplot()

        # take the data of the given index
        data_loc = self._data.data.loc[index]

        # iterate over all surveys
        for s in self._data._survey.get_surveys():
            sp.errorbar(self._data._survey.get_survey_wavelengths(s),
                        data_loc[self._data._survey.get_survey_magnitude(s)].values,
                        data_loc[self._data._survey.get_survey_error_names(s)].values,
                        fmt='.', capsize=2,
                        label=s)

        data_loc = data_loc[self._data._survey.all_magnitudes()].values
        mask = data_loc > 0
        if fit is not None:
            if type(fit) == str:
                if fit == 'blackbody':
                    pass
                else:
                    raise ValueError('Unknown function\'{}\''.format(fit))

            # if the fit-parameter is a function
            # use the function as fitting function
            elif callable(fit):
                wavelengths = self._data._survey.get_all_wavelengths()
                popt, pcov = curve_fit(fit,
                                       wavelengths[mask],
                                       data_loc[mask])
                waves = np.linspace(np.min(wavelengths),
                                    np.max(wavelengths),
                                    1000)
                # plot the fitted curve as a black dashed line
                sp.plot(waves, fit(waves, *popt), '--k', legend='fit')
            elif type(fit) == int:
                wavelengths = self._data._survey.get_all_wavelengths()
                print(data_loc.shape)
                print(mask.shape)
                print(wavelengths.shape)
                print(wavelengths[mask])
                print(data_loc[mask])
                fit = np.polyfit(wavelengths[mask], data_loc[mask], fit)
                poly = np.poly1d(fit)
                waves = np.linspace(np.min(wavelengths),
                                    np.max(wavelengths),
                                    1000)
                sp.plot(waves, poly(waves), '--k')

            else:
                raise ValueError('Only functions and strings are allowed as fit-parameter values.')

        sp.set_xlabel('wavelength [$\\AA$]')
        sp.set_ylabel('flux')

        if legend:
            pl.legend(loc='best')

        if x_log:
            sp.set_xscale('log')

        if y_log:
            sp.set_yscale('log')

        # if no path is given, show the figure
        if path == '':
            pl.show()
        # if a path is given, save the figure
        else:
            pl.savefig(path)
