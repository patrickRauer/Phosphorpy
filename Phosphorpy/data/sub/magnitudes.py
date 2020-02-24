import configparser
import sqlite3
import warnings

import numba
import numpy as np
from armapy.svo import get_survey_filter_information
from astropy.table import Table
from pandas import DataFrame

from Phosphorpy.core.structure import Mask
from .colors import Colors
from .flux import FluxTable
from .plots.magnitude import MagnitudePlot
from .table import DataTable
from .tables.magnitude import Magnitude

try:
    from Phosphorpy.data.sub.interactive_plotting.mangitude import MagnitudePlot as MagnitudePlotHV
except ImportError:
    MagnitudePlotHV = None


def power_2_10(x):
    return np.power(10., -x / 2.5)


def remove_astrometry_columns(mags):
    cols = []
    for c in mags.columns:
        cl = c.lower()
        if 'ra' in cl or 'de' in cl or 'plx' in cl:
            continue
        cols.append(c)
    return cols


@numba.vectorize
def subtract(a, b):
    """
    Numba implementation of a subtraction

    :param a:
    :param b:
    :return:
    """
    return a - b


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
                    (prefix_cond and (c == s or 'mag' in c or 'Petro' in c))):
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
    skymapper = 0

    # count how often a band could be part of a survey
    for mag in mag_cols:
        # ignore error columns or J2000 coordinate columns
        if 'e_' in mag or 'J2000' in mag or 'deg' in mag or 'ra' in mag or 'dec' in mag:
            continue
        # strip mag postfix
        if 'mag' in mag:
            mag = mag.split('mag')[0]
        elif 'Petro' in mag:
            mag = mag.split('Petro')[0]

        if 'u' in mag:
            kids += 1
            sdss += 1
            skymapper += 1
        elif 'g' in mag or 'r' in mag or 'i' in mag:
            if 'p' in mag:
                apass += 1
            else:
                pan_starrs += 1
                sdss += 1
                kids += 1
                skymapper += 1

        elif 'z' in mag:
            pan_starrs += 1
            sdss += 1
            skymapper += 1
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
        elif 'v' in mag:
            skymapper += 1
    # TODO: use data from local file
    # decide which survey it is
    if skymapper == 6:
        kids -= 4
        sdss -= 5
        pan_starrs -= 4
        prefix = get_prefix(mag_cols, 'v')
        s_cols = ['u', 'v', 'g', 'r', 'i', 'z']
        mag_cols, cols = get_survey_cols(mag_cols, s_cols, prefix)
        return cols

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
    Returns only the columns names of the magnitudes by excluding columns with "e_" in it.

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
    _properties = None

    def __init__(self, name=None, mag_cols=None):
        self._survey_cols = {}
        self._properties = {}
        if name is not None:
            mag_cols = mag_cols_only(mag_cols)
            self.__add_survey_properties__(name, mag_cols)

    def __add_survey_properties__(self, survey, cols):
        survey = survey.lower()
        survey_properties = {}
        for c in cols:
            if 'mag' in c:
                c = c.split('mag')[0]
            elif 'Petro' in c:
                c = c.split('Petro')[0]
            survey_properties[c] = get_survey_filter_information(survey, c)
        self._properties[survey] = survey_properties
        self._survey_cols.update({survey: cols})

    def __getitem__(self, item):
        item = item.lower()
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

    def get_survey_id(self, survey_name):
        """
        Returns the index of the survey with the given name. If no survey with such a name is in the data,
        a ValueError raises.

        :param survey_name: The name of the survey
        :type survey_name: str
        :return: The index of the survey
        :rtype: int
        """
        survey_name = survey_name.lower()
        try:
            return np.where(
                np.array(list(self._survey_cols.keys())) == survey_name)[0][0]
        except IndexError:
            raise ValueError(f'Survey {survey_name} not known.')

    def get_corresponding_survey(self, mag_name):
        """
        Returns the first survey with a magnitude with such a name.
         If no survey with such magnitude name is found, it returns None.

        :param mag_name: The name of the magnitude
        :type mag_name: str
        :return: The survey
        :rtype:
        """
        for s in self._survey_cols:
            if mag_name in self._survey_cols[s]:
                return s

    def get_corresponding_survey_id(self, mag_name):
        return self.get_survey_id(self.get_corresponding_survey(mag_name))

    def get_survey_magnitude(self, survey_name):
        """
        Returns the names of the magnitudes of the different surveys

        :param survey_name: The name of the survey
        :type survey_name: str
        :return: A list with the names of the magnitudes
        :rtype: list
        """
        survey_name = survey_name.lower()
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
        """
        Returns the names of the error magnitude columns of  the required survey

        :param survey_name: The name of the survey
        :type survey_name: str
        :return:
            The names of the error columns, if such a survey exists.
            If no survey with such a name exists, an empty list will be returned.
        :rtype: list
        """
        survey_name = survey_name.lower()
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
        return list(self._survey_cols.keys())

    def add_survey(self, name, mag_cols):
        """
        Adds a new survey with the magnitude names to the current list.

        :param name: The name of the new survey
        :param name: str
        :param mag_cols: A list with the magnitude names
        :type mag_cols: list
        :return:
        """
        name = name.lower()
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
        survey = survey.lower()
        lambda_eff = self._properties[survey][band]['lambda_eff']
        try:
            return float(lambda_eff)
        except ValueError:
            return lambda_eff

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
        return np.float64(out)

    def get_survey_wavelengths(self, survey_name):
        """
        Returns the effective wavelengths of the bands of the survey sorted by their values.

        :param survey_name: The name of the survey
        :type survey_name: str
        :return: The wavelengths
        :rtype: numpy.ndarray
        """
        survey_name = survey_name.lower()

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
        if 'vega' in props.keys() and bool(props['vega']):
            return float(props['Vega_ergs'])
        return float(props['AB_ergs'])

    def is_vega(self, survey, band):
        """
        Returns true, if the photometric system for this survey and band is Vega, else False
        :param survey: The name of the survey
        :type survey: str
        :param band: The name of the band
        :type band: str
        :return: True, if the photometric system is vega, else False
        :rtype: bool
        """
        survey = survey.lower()
        return bool(self._properties[survey][band])

    def _to_config(self, path):
        """
        Writes the survey properties as a config file

        :param path: The path to the file
        :type path: str
        :return: None
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

    def to_dataframe(self):
        """
        Converts the survey properties to a single pandas DataFrame.

        :return: All available survey properties
        :rtype: DataFrame
        """
        surveys = self.get_surveys()
        data = {'survey': [],
                'band': []}
        for s in surveys:
            bands = self._properties[s]
            for b in bands:
                data['survey'].append(s)
                data['band'].append(b)
                props = bands[b]
                for p in list(props):
                    if p in data.keys():
                        data[p].append(props[p])
                    else:
                        data[p] = [props[p]]
        return DataFrame(data)

    def add_to_meta(self, meta, name):
        bands = self._properties[name]
        meta['name'] = name
        for b in bands:
            props = bands[b]
            for p in list(props):
                meta[f'{b}_{p}'] = props[p]

    def _to_sql(self, path):
        """
        Writes the data to a sql-database

        :param path: Path to the sql-database
        :type path: str
        :return:
        """
        con = sqlite3.Connection(path)
        data = self.to_dataframe()
        data.to_sql('survey_properties', con)

    def write(self, path, data_format='config'):
        """
        Writes the survey data to a config-file

        :param path: Path to the config file
        :type path: str
        :param data_format: Not used
        :return:
        """
        data_format = data_format.lower()
        if data_format == 'config' or data_format == 'init':
            self._to_config(path)
        elif data_format == 'sql':
            self._to_sql(path)
        else:
            raise ValueError(f'Unsupported data format: {data_format}')

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
                    key = c.split(band + '_')[-1]
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


def __check_input__(data, names):
    """
    Checks the input data if they have a proper format and covert them (if necessary) to a
    pandas DataFrame

    :param data: The input data
    :param names: The names of the columns
    :return: A DataFrame with the input data
    :rtype: pandas.DataFrame
    """
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
    return data


# class SurveyMag(DataTable):
#     pass
#
#     def __init__(self, data, names=None, survey='', mask=None):
#         DataTable.__init__(self, mask=mask)


class MagnitudeTable(DataTable):
    """
    Class to handle any interactions with magnitudes of multiple surveys
    """
    _mag_cols = None
    _err_cols = None
    _survey = None

    def __init__(self, data=None, names=None, survey='', mask=None):
        DataTable.__init__(self, mask=mask)
        self._mag_cols = []
        self._err_cols = []
        self._plot = MagnitudePlot(self)
        data = __check_input__(data, names)

        if MagnitudePlotHV.holoviews():
            self._hv_plot = MagnitudePlotHV(self)

        try:
            cols = guess_surveys(data.columns)
            if cols is not None and survey != '':
                self._survey = SurveyData(survey, cols)
                self._set_cols(cols)
                data = data[cols]
                if type(data) != Magnitude:
                    data = Magnitude(data, survey, mask=mask)
                self._data = [data]
        except ValueError:
            warnings.warn('No magnitude columns!')
        except AttributeError:
            pass

    def __str__(self):
        string = 'Available magnitudes\n'
        for d in self._data:
            string += str(d)
        return string

    def __getitem__(self, item):
        if type(item) == str:
            return self.get_survey_data(item)
        raise ValueError('No indices supported')

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

    def apply_extinction_correction_to_survey(self, correction, survey):
        """
        Applies the extinction correction to one specific survey

        :param correction:
            Table with the extinction correction.
            The number of rows must be the same as the number of detections in the survey and the columns
            must have the same names as the magnitude columns in the survey
        :type correction: Table, DataFrame
        :param survey: The name of the survey
        :type survey: str
        :return:
        """
        for d in self.data:
            if d.survey_name == survey:
                d.apply_extinction_correction(correction)
                break

    def apply_extinction_correction(self, correction):
        """
        Applies the extinction correction to the data

        :param correction: A table with the correction values
        :type correction: astropy.table.Table, pandas.DataFrame
        :return:
        """
        for d in self.data:
            d.apply_extinction_correction(correction)

    def get_names(self):
        """
        Returns the names of the magnitude columns

        :return: The name of magnitude columns
        :rtype: list
        """
        names = [d.survey_name for d in self.data]
        return names

    def get_flux(self):
        """
        Returns the flux of the columns
        :return:
        """
        if self.data is None:
            raise AttributeError('No photometric data available!\nDownload data first!')
        flux = []
        for d in self.data:
            flux.append(d.get_flux())

        # return a merged table out of the fluxes and the flux errors
        return FluxTable(flux,
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

        colors = []
        for d in self._data:
            colors.append(d.get_colors())

        return Colors(colors, mask=self._mask)

    def create_mask(self, overwrite=False):
        """
        Creates a new mask object for all data

        :param overwrite:
            True, if the old mask should be replaced, else False. Default is false.
            If a mask is set already and overwrite False, a RuntimeError will raise.
        :type overwrite: bool
        :return:
        """
        if self._mask is not None and not overwrite:
            raise RuntimeError('A mask is set and overwrite is false.')

        if self._data is None:
            raise ValueError("No data were set.")

        indices = []
        for d in self._data:
            indices.extend(d.data.index.values)
        indices = np.unique(indices)

        self._mask = Mask(len(indices), ids=indices)

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
        mag_cols = self.survey.get_survey_magnitude(survey)
        d = self._data[self.survey.get_survey_id(survey)]

        for c in mag_cols:

            k = d[c] == d[c]
            if mask is None:
                mask = k
            else:
                mask = mask & k
        try:
            self.mask.add_mask(mask, 'Mask sources with full photometry in {}'.format(survey), combine=previous)
        except AttributeError:
            warnings.warn('No mask was set. Call \'create_mask()\' to create a mask object.')

    def set_limit(self, band, survey=None, minimum=99, maximum=-99, previous=True):
        """
        Sets a magnitude limit to the magnitude columns with the name band.

        .. code-block:: python

            # select all sources with a G magnitude between 16 and 18
            # (notice that minimum and maximum are magnitudes, therefore the minimum value is larger than the maximum)
            ds.magnitudes.set_limit('G', minimum=18, maximum=16)


        :param band: The name of the band/magnitude column
        :type band: str
        :param survey: 
            The name of the specific survey or None. If None all surveys with such a magnitude name will
            used for the limiting
        :type survey: str, None
        :param minimum: The minimal magnitude value. Default is 99, which means that no cut will be done.
        :type minimum: float
        :param maximum: The maximal magnitude value. Default is -99, which means that no cut will be done.
        :type maximum: float
        :param previous: True if the previous mask should be used True, else False. Default is True.
        :type previous: bool
        :return:
        """
        for d in self._data:
            if survey is None or d.has_name(survey):
                try:
                    d.set_limit(band, minimum=minimum, maximum=maximum, previous=previous)
                except ValueError:
                    pass

    def get_survey_data(self, name):
        """
        Returns the survey data of the survey with the given name.

        :param name: The name of the survey
        :type name: str
        :return: The survey data
        :rtype: SurveyData
        """
        for d in self._data:
            if d.has_name(name):
                return d

    def get_columns(self, cols):
        """
        Creates a (sub)set of the data with the given columns

        :param cols: The name of the magnitude columns
        :type cols: list, str
        :return: DataFrame with the magnitude information or None if no columns were vound
        :rtype: pd.DataFrame
        """
        if type(cols) == str:
            cols = [cols]
        elif type(cols) != list:
            raise AttributeError('"cols" must be a string or a list.')

        out = None
        for d in self.data:
            rs = d.get_columns(cols)
            if len(rs.columns) > 0:
                if out is None:
                    out = rs
                else:
                    out = out.join(rs, rsuffix=f'_{d.survey_name}')
        return out

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

    @property
    def head(self):
        return self._survey

    def add_costum_magnitudes(self, mags, name, survey):
        """
        costum
        """

        if len(mags) == 0:
            raise ValueError('No magnitudes are inside the given data.')

        # get the columns of the given data without the astrometry columns
        cols = remove_astrometry_columns(mags)

        magnitudes = Magnitude(mags[cols], cols, name, self.mask)
        magnitudes.select_columns(cols)

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
        # get the columns of the given data without the astrometry columns
        cols = remove_astrometry_columns(mags)

        mags = Magnitude(mags, cols, survey, mask=self.mask)
        mags.select_columns(cols)

        if len(mags) == 0:
            warnings.warn(f'No data found in {survey} for the given coordinates.\n{survey} will not be listed!')
            return

        survey = survey.lower()

        if self._data is not None:
            try:
                self._data.append(mags)
                self._survey.add_survey(survey, cols)
            # if the magnitudes include magnitudes with the same name already
            except ValueError:
                # create a map to rename the all columns of the new survey with the survey name
                # as postfix
                name_map = {}
                for c in cols:
                    name_map[c] = '{}_{}'.format(c, survey)

                mags.rename(name_map, axis='columns')

                # join the data
                self._data.append(mags)
                self._survey.add_survey(survey, list(name_map.values()))
        # if no magnitudes are set
        else:
            self._data = [mags]

            self._survey = SurveyData(survey, mags.columns)
        mags.set_survey_data(self._survey)

    def write(self, path, data_format='parquet', **kwargs):
        out = []
        for d in self.data:
            out.append(
                d.write(path, data_format=data_format)
            )
        return out
