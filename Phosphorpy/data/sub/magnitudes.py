from pandas import DataFrame
from astropy.table import Table
from .table import DataTable
from .plots.magnitude import MagnitudePlot
from .colors import Colors
from .flux import FluxTable
import numpy as np
import configparser
import numba
import warnings

from armapy.svo import get_survey_filter_information


def power_2_10(x):
    return np.power(10., -x/2.5)


@numba.vectorize
def subtract(a, b):
    """
    Numba implementation of a subtraction

    :param a:
    :param b:
    :return:
    """
    return a-b


def get_prefix(cols, unique_col):
    """
    Extracts a prefix from the column names

    :param cols: The name of the columns
    :type cols: list
    :param unique_col: A unique column
    :type unique_col: str
    :return: The prefix of all columns
    :rtype: str
    """
    for c in cols:
        if unique_col in c:
            return c.split(unique_col)[0]
    return ''


def get_survey_cols(cols, s_cols, prefix):
    """
    Returns the names of the magnitudes in a specific survey

    :param cols:
    :param s_cols:
    :param prefix:
    :return:
    """
    out = []
    prefix_cond = len(prefix) == 0
    for i, c in enumerate(cols):
        for j, s in enumerate(s_cols):
            # if the survey magnitude name is in and the prefix is in
            # of if there is no prefix and the magnitude name is equal to the survey name
            if ((s in c and prefix in c and not prefix_cond) or
                    (prefix_cond and (c == s or 'mag' in c))):
                out.append(c)
                break

    cols_new = []
    for i, c in enumerate(cols):
        if c not in out:
            cols_new.append(c)
    return cols_new, out


def guess_surveys(mag_cols):
        """
        Try to find matches of the magnitude column names with a few
        implement large sky surveys.

        example:
            SDSS has the column u, g, r, i and z and PAN-STARRS g, r, i, z, y
            so both have 5 bands but PAN-STARRS doesn't have u and SDSS
            doesn't have y.
        :return:
        """
        pan_starrs = 0
        two_mass = 0
        sdss = 0
        kids = 0
        viking = 0
        wise = 0
        apass = 0
        gaia = 0
        galex = 0
        ukidss = 0

        # count how often a band could be part of a survey
        for mag in mag_cols:
            # ignore error columns or J2000 coordinate columns
            if 'e_' in mag or 'J2000' in mag or 'deg' in mag:
                continue
            # strip mag postfix
            mag = mag.split('mag')[0]
            if 'u' in mag:
                kids += 1
                sdss += 1
            elif 'g' in mag or 'r' in mag or 'i' in mag:
                if 'p' in mag:
                    apass += 1
                else:
                    pan_starrs += 1
                    sdss += 1
                    kids += 1

            elif 'z' in mag:
                pan_starrs += 1
                sdss += 1
            elif 'y' in mag:
                pan_starrs += 1
            elif 'Y' in mag or 'Z' in mag:
                viking += 1
                ukidss += 1
            elif 'J' in mag or 'H' in mag or 'K' in mag:
                viking += 1
                two_mass += 1
                ukidss += 1
            elif 'W' in mag:
                wise += 1
            elif 'G' in mag or 'RP' in mag or 'BP' in mag:
                gaia += 1
            elif 'FUV' in mag or 'NUV' in mag:
                galex += 1
            elif 'B' in mag or 'V' in mag:
                apass += 1
        # TODO: use data from local file
        # decide which survey it is
        if pan_starrs == 5 or (pan_starrs > 5 and 'y' in mag_cols):
            kids -= 3
            sdss -= 4

            prefix = get_prefix(mag_cols, 'y')

            s_cols = ['g', 'r', 'i', 'z', 'y']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            return cols

        if sdss == 5 or sdss == 9:
            kids -= 4

            prefix = get_prefix(mag_cols, 'z')

            s_cols = ['u', 'g', 'r', 'i', 'z']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            return cols

        if kids == 4:
            prefix = get_prefix(mag_cols, 'u')

            s_cols = ['u', 'g', 'r', 'i']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            return cols

        if viking == 5 or (viking > 5 and 'Z' in mag_cols):
            prefix = get_prefix(mag_cols, 'Z')

            s_cols = ['Z', 'Y', 'J', 'H', 'K']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            two_mass -= 3
            ukidss -= 4
            return cols

        if ukidss == 3:
            # todo: update after J band is included
            prefix = get_prefix(mag_cols, 'Y')

            s_cols = ['Y', 'H', 'K']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            two_mass -= 2
            return cols

        if two_mass == 3:
            prefix = get_prefix(mag_cols, 'H')
            s_cols = ['J', 'H', 'K']
            mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
            return cols

        if wise == 4:
            cols = []
            for i, c in enumerate(mag_cols):
                if 'W' in c:
                    cols.append(c)
            return cols

        if gaia == 3:
            cols = []
            n_cols = []
            for i, c in enumerate(mag_cols):
                if 'G' in c or 'BP' in c or 'RP' in c:
                    cols.append(c)
                else:
                    n_cols.append(c)
            return cols

        if apass == 5:
            cols = []
            n_cols = []
            for i, c in enumerate(mag_cols):
                if 'gp' in c or 'rp' in c or 'ip' in c or 'B' in c or 'V' in c:
                    cols.append(c)
                else:
                    n_cols.append(c)
            return cols

        if galex == 2:
            cols = []
            n_cols = []
            for i, c in enumerate(mag_cols):
                if 'FUV' in c or 'NUV' in c:
                    cols.append(c)
                else:
                    n_cols.append(c)
            return cols


def mag_cols_only(cols):
    """
    Returns only the columns names of the magnitudes by excluding columns with 'e_' in it.
    :param cols: A list with all column names
    :type cols: list
    :return: A list with the magnitude names only
    :rtype: list
    """
    out = []
    for c in cols:
        if 'e_' not in c:
            out.append(c)
    return out


class SurveyData:
    """
    Class to handle the properties of different surveys, like providing the information of the different
    bands (flux zero point, effective wavelength, etc).
    """
    _survey_cols = None
    _properties = {}

    def __init__(self, name=None, mag_cols=None):
        if name is not None:
            mag_cols = mag_cols_only(mag_cols)
            self._survey_cols = {name: mag_cols}
            self.__add_survey_properties__(name, mag_cols)
        else:
            self._survey_cols = {}

    def __add_survey_properties__(self, survey, cols):
        survey_properties = {}
        for c in cols:
            c = c.split('mag')[0]
            survey_properties[c] = get_survey_filter_information(survey, c)
        self._properties[survey] = survey_properties

    def __getitem__(self, item):
        return self._survey_cols[item]

    def all_magnitudes(self):
        """
        Returns all magnitude column names
        :return:
        """
        out = []
        for s in self._survey_cols:
            out.extend(self._survey_cols[s])
        return out

    def get_survey_magnitude(self, survey_name):
        """
        Returns the names of the magnitudes of the different surveys

        :param survey_name: The name of the survey
        :type survey_name: str
        :return: A list with the names of the magnitudes
        :rtype: list
        """
        return self._survey_cols[survey_name]

    def all_error_names(self):
        """
        Returns all error column names
        :return:
        """
        out = []
        for s in self._survey_cols:
            for c in self._survey_cols[s]:
                out.append('e_{}'.format(c))
        return out

    def get_survey_error_names(self, survey_name):
        out = self.get_survey_magnitude(survey_name)
        o = []
        for k in out:
            o.append('e_{}'.format(k))
        return o

    def get_surveys(self):
        """
        Returns the names of the surveys
        :return:
        """
        return self._survey_cols.keys()

    def add_survey(self, name, mag_cols):
        """
        Adds a new survey with the magnitude names to the current list.

        :param name: The name of the new survey
        :param name: str
        :param mag_cols: A list with the magnitude names
        :type mag_cols: list
        :return:
        """
        mag_cols = mag_cols_only(mag_cols)
        self._survey_cols[name] = mag_cols
        self.__add_survey_properties__(name, mag_cols)

    def effective_wavelength(self, survey, band):
        """
        Returns the effective wavelength of the filter

        :param survey: The name of the survey
        :type survey: str
        :param band: The filter band of the survey
        :type band: str
        :return: The effective wavelength
        :rtype: float
        """
        return self._properties[survey][band]['lambda_eff']

    def get_all_wavelengths(self):
        """
        Returns all effective wavelengths at once.
        The order is the same order as the magnitude columns.

        :return: All effective wavelengths
        :rtype: numpy.ndarray
        """
        out = []
        for s in self._properties:
            survey = self._properties[s]
            for band in survey:
                out.append(self.effective_wavelength(s, band))
        return np.array(out)

    def get_survey_wavelengths(self, survey_name):
        out = []
        survey = self._properties[survey_name]
        for band in survey:
            out.append(self.effective_wavelength(survey_name, band))
        return np.array(out)

    def flux_zero(self, survey, band):
        """
        Returns the flux zero point of the filter.
        For 2MASS, VIKING, GALEX and Gaia the Vega zero point is used and for all the other
        ones, the AB zero point

        :param survey: The name of the survey
        :type survey: str
        :param band: The filter band of the survey
        :type band: str
        :return: The flux zero point in ergs
        :rtype: float
        """
        band = band.split('mag')[0]
        props = self._properties[survey][band]
        # check if in the SVO data the vega system is set as the used system
        if 'vega' in props.keys() and props['vega']:
            return props['Vega_ergs']
        return props['AB_ergs']

    def write(self, path, data_format=''):
        """
        Writes the survey data to a config-file

        :param path: Path to the config file
        :type path: str
        :param data_format: Not used
        :return:
        """
        conf = configparser.ConfigParser()
        for survey in self.get_surveys():
            conf.add_section(survey)
            conf.set(survey, 'magnitudes', str(self.get_survey_magnitude(survey)))
            for props in self._properties[survey]:
                proper = self._properties[survey][props]
                for p in proper:
                    conf.set(survey, '{}_{}'.format(props, p), str(proper[p]))
        with open(path, 'w') as f:
            conf.write(f)
        _survey_cols = None
        _properties = {}

    @staticmethod
    def read(path):
        """
        Reads survey data from a config file

        :param path: Path to the config file
        :type path: str
        :return: The survey data from the config file
        :rtype: SurveyData
        """
        conf = configparser.ConfigParser()
        conf.read(path)

        sd = SurveyData()
        for s in conf.sections():
            mag_cols = conf[s]['magnitudes']
            mag_cols = mag_cols.split('[\'')[-1].split('\']')[0]
            mag_cols = mag_cols.split('\', \'')
            props = {}
            for c in conf[s]:
                if 'magnitudes' not in c:
                    band = c.split('_')[0]
                    key = c.split(band+'_')[-1]
                    if band not in props:
                        props[band] = {}
                    props[band][key] = conf[s][c]

            sd._survey_cols[s] = mag_cols
            sd._properties[s] = props

        return sd

    def __str__(self):
        s = ''
        if self._survey_cols is None:
            return 'No information'
        for k in self._survey_cols:
            cols = self._survey_cols[k]
            s += '{}\n'.format(k)
            for c in cols:
                s += '\t{}\n'.format(c)
        return s


class Magnitude(DataTable):
    """
    Class to handle any interactions with magnitudes of multiple surveys
    """
    _mag_cols = []
    _err_cols = []
    _survey = None

    def __init__(self, data, names=None, survey='', mask=None):
        DataTable.__init__(self, mask=mask)
        self._plot = MagnitudePlot(self)
        if type(data) != DataFrame:
            if type(data) == Table:
                data = data.to_pandas()
            elif type(data) == np.ndarray:
                if len(data.shape) == 1:
                    data = DataFrame.from_records(data)
                else:
                    if names is None:
                        raise ValueError('If the input data are in numpy.ndarray, names are needed!')
                    elif len(names) != data.shape[1]:
                        raise AttributeError('Number of columns does not the number of names!')
                    else:
                        data = DataFrame(data=data, columns=names)

        try:
            cols = guess_surveys(data.columns)
            if cols is not None and survey != '':
                self._survey = SurveyData(survey, cols)
                self._set_cols(cols)
                self._data = data[cols]
        except ValueError:
            warnings.warn('No magnitude columns!')

    def _set_cols(self, cols):
        mag = []
        err = []
        for c in cols:
            if 'e_' in c:
                err.append(c)
            else:
                mag.append(c)
        self._mag_cols = mag
        self._err_cols = err

    def min(self):
        """
        Returns the minimal values of all magnitudes (the magnitudes can come from different sources)

        :return: A pandas series with the column names as indices and the minimal values
        :rtype: pandas.core.series.Series
        """
        return self.apply(np.min)

    def max(self):
        """
        Returns the maximal values of all magnitudes (the magnitudes can come from different sources)

        :return: A pandas series with the column names as indices and the minimal values
        :rtype: pandas.core.series.Series
        """
        return self.apply(np.max)

    def get_names(self):
        """
        Returns the names of the magnitude columns

        :return: The name of magnitude columns
        :rtype: list
        """
        return list(self.data.columns)

    def get_flux(self):
        """
        Returns the flux of the columns
        :return:
        """
        # transform the magnitude in the flux scale
        flux = self.data[self._survey.all_magnitudes()].apply(power_2_10)

        # multiply the flux zero point of all bands of the different surveys
        for s in self._survey.get_surveys():
            for m in self._survey[s]:
                fo = self._survey.flux_zero(s, m)
                flux[m] *= fo

        # compute the errors of the fluxes
        errs = self.data[self._survey.all_error_names()]
        errs = DataFrame(data=flux.values*errs.values*np.log(10.),
                         columns=errs.columns)

        # return a merged table out of the fluxes and the flux errors
        return FluxTable(flux.merge(errs, left_index=True, right_index=True),
                         survey_head=self._survey, mask=self._mask)

    def get_colors(self, cols=None):
        """
        Computes all possible colors. The first magnitude will be the bluer and the second one will be the redder one,
        except magnitude columns are given, then this order is not fixed.

        .. math:

            C_{i, j} = mag_i - mag_j \\
            \\lambda(band_i) < \\lambda(band_j)

        :param cols:
        :return:
        """
        if cols is None:
            cols = []
            for s in self._survey.get_surveys():
                cols.append((s, self._survey.get_survey_magnitude(s)))
        else:
            try:
                cols[0][0][0]
            except TypeError:
                cols = [cols]
        all_cols = cols
        colors = DataFrame()
        survey_c = {}
        if len(self.data) < 1e4:
            for cols in all_cols:
                survey = cols[0]
                cols = cols[-1]
                color_color_cols = []
                for i, c1 in enumerate(cols):
                    for j, c2 in enumerate(cols):
                        # compute only such colors where the first argument is in a bluer band
                        # (assumption is that the columns are sorted in such a way)
                        if i < j:
                            cc_name = '{} - {}'.format(c1, c2)
                            color_color_cols.append(cc_name)
                            colors[cc_name] = self.data[c1] - self.data[c2]
                survey_c[survey] = color_color_cols
        else:
            for cols in all_cols:
                survey = cols[0]
                cols = cols[-1]
                color_color_cols = []
                for i, c1 in enumerate(cols):
                    for j, c2 in enumerate(cols):
                        # compute only such colors where the first argument is in a bluer band
                        # (assumption is that the columns are sorted in such a way)
                        if i < j:
                            cc_name = '{} - {}'.format(c1, c2)
                            color_color_cols.append(cc_name)
                            colors[cc_name] = subtract(self.data[c1].values, self.data[c2].values)
                survey_c[survey] = color_color_cols
        return Colors(colors, mask=self._mask, survey_colors = survey_c)

    def has_full_photometry(self, survey, previous=True):
        """
        Creates a mask of sources, which have photometry values in every band of the survey.

        :param survey: The name of the survey
        :type survey: str
        :param previous: True if the previous mask should be used True, else False. Default is True.
        :type previous: bool
        :return:
        """
        mask = None
        for c in self._survey.get_survey_magnitude(survey):
            k = self._data[c] > -99
            if mask is None:
                mask = k
            else:
                mask = mask & k
        self.mask.add_mask(mask, 'Mask sources with full photometry in {}'.format(survey), combine=previous)

    def set_limit(self, band, minimum=99, maximum=-99, previous=True):
        """
        Sets a magnitude limit to the magnitude columns with the name band.

        :param band: The name of the band/magnitude column
        :type band: str
        :param minimum: The minimal magnitude value. Default is 99, which means that no cut will be done.
        :type minimum: float
        :param maximum: The maximal magnitude value. Default is -99, which means that no cut will be done.
        :type maximum: float
        :param previous: True if the previous mask should be used True, else False. Default is True.
        :type previous: bool
        :return:
        """
        if band not in self._survey.all_magnitudes():
            band_new = band + 'mag'
            if band_new not in self._survey.all_magnitudes():
                raise ValueError('Magnitude column {} or {} not found.'.format(band, band_new))
            band = band_new
        cut_string = ''
        if minimum != 99:
            cut_string = 'minimum={}'.format(minimum)
        if maximum != -99:
            max_str = 'maximum={}'.format(maximum)
            if cut_string == '':
                cut_string = max_str
            else:
                cut_string = '{} and {}'.format(cut_string, max_str)

        mags = self._data[band].values
        if maximum != -99 and minimum != 99:
            mask = (mags <= minimum) & (mags >= maximum)
        elif maximum != -99:
            mask = mags >= maximum
        elif minimum != 99:
            mask = mags <= 99
        else:
            raise AttributeError('At least one of maximum or minimum must be different to the default value.')
        self.mask.add_mask(mask, 'Apply a magnitude cut on the band {} with {}'.format(band, cut_string),
                           combine=previous)

    @property
    def mag_names(self):
        return self._mag_cols

    @property
    def err_names(self):
        return self._err_cols

    @property
    def survey(self):
        return self._survey

    @survey.setter
    def survey(self, value):
        self._survey = value

    def add_survey_mags(self, mags, survey):
        """
        Adds a new set of magnitudes from another survey. If already magnitudes with the same name are
        in the current data, all the new columns will be renamed to the magnitude name plus the
        survey name.

        :param mags: The new magnitudes
        :param survey: The name of the new survey
        :type survey: str
        :return:
        """
        if type(mags) is not DataFrame:
            mags = DataFrame(mags)
        cols = guess_surveys(mags.columns)
        mags = mags[cols]
        survey = survey.lower()
        if self._data is not None:
            try:
                self._data = self._data.join(mags, how='outer')
                self._survey.add_survey(survey, cols)
            # if the magnitudes include magnitudes with the same name already
            except ValueError:
                # create a map to rename the all columns of the new survey with the survey name
                # as postfix
                name_map = {}
                for c in cols:
                    name_map[c] = '{}_{}'.format(c, survey)

                mags = mags.rename(name_map, axis='columns')

                # join the data
                self._data = self._data.join(mags, how='outer')
                self._survey.add_survey(survey, list(name_map.values()))
        # if no magnitudes are set
        else:
            self._data = mags
            self._survey = SurveyData(survey, mags)
