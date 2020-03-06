
from threading import Thread

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astroquery.vizier import Vizier as Viz
from astroquery.simbad import Simbad
from pandas import DataFrame
from ..config.survey_data import SURVEY_DATA
from ..config import names
from .xmatch import xmatch
from ..config.match import match_catalogs

import numpy as np
import warnings


warnings.simplefilter('ignore')

Viz.ROW_LIMIT = -1
Simbad.add_votable_fields('otype')


class Vizier:
    columns = ['*']
    catalog = ''
    radius = 2*u.arcsec
    vizier = Viz
    name = ''
    ra_name = 'RAJ2000'
    dec_name = 'DEJ2000'
    cache = None
    reference = ''
    release = ''

    def __init__(self, catalog='', name=None, columns=None):
        """
        Basic class for an interface between the Vizier class of the astroquery package

        :param catalog: The name of the catalog in the Vizier-database
        :type catalog: str
        """
        if name is not None:
            self.name = name
        self.catalog = catalog
        selected_survey = SURVEY_DATA[self.name]
        self.reference = selected_survey['reference']
        self.catalog = selected_survey['vizier']
        self.release = selected_survey['release']
        cols = []

        # check if coordinate names are set
        # if yes, take set coordinate names
        if 'coordinate' in selected_survey.keys():
            coordinate_names = selected_survey['coordinate']
        # if not take the default ones
        else:
            coordinate_names = [self.ra_name, self.dec_name]
        cols.extend(coordinate_names)

        # transform the filter names to the vizier style
        # except it is GALEX
        if self.name != names.GALEX:
            for c in selected_survey['magnitude']:
                cols.append(f'{c}mag')
                cols.append(f'e_{c}mag')
        else:
            for c in selected_survey['magnitude']:
                cols.append(c)
                cols.append(f'e_{c}')

        if columns is None:
            # if there are additional columns set
            if 'columns' in selected_survey.keys():
                cols.extend(selected_survey['columns'])
            self.columns = cols
        else:
            self.columns = columns
        self.vizier = Viz(columns=self.columns)
        self.vizier.ROW_LIMIT = -1

    def __adjust_coordinate_names__(self):
        """
        Reformat the coordinate labeling to avoid troubles with similar names
        :return:
        """
        if self.name != '' and self.name not in self.ra_name:
            self.ra_name = f'{self.name}_{self.ra_name}'
            self.dec_name = f'{self.name}_{self.dec_name}'

    def __query__(self, coordinates):
        """
        Perform the query with a split of the coordinates in 200 item
        blocks to reduce a problem with the connection

        :param coordinates: The list of coordinates of the targets
        :type coordinates: SkyCoord
        :return: A array or a list of arrays with the results of the queries
        :rtype: list, numpy.ndarray
        """

        self.__adjust_coordinate_names__()
        # set how many coordinates should be queried at the same time
        steps = 200
        out = []

        if type(self.catalog) == list:
            for i in range(len(self.catalog)):
                out.append([])
        else:
            out = [[]]

        if coordinates.isscalar:
            coordinates = SkyCoord([coordinates.ra.degree],
                                   [coordinates.dec.degree],
                                   unit=(u.deg, u.deg))
        length_coordinates = len(coordinates)

        # perform the query with maximal 200 targets at the same time
        # this is necessary because otherwise it happens from time to time that a timeout appears
        i = 0
        while i*steps < length_coordinates:
            coords = coordinates[i*steps: (i+1)*steps]
            rs = self.vizier.query_region(coords, self.radius,
                                          catalog=self.catalog)
            for j, r in enumerate(rs):
                out[j].append(r)
            i += 1

        # stack the results
        for i, out_sec in enumerate(out):
            if len(out_sec) != 0:
                o = vstack(out_sec)

                if self.name == 'GAIA':
                    ra_name = 'RA_ICRS'
                    dec_name = "DE_ICRS"
                elif self.name == 'SDSS' or self.name == 'SkyMapper':
                    ra_name = '_RAJ2000'
                    dec_name = '_DEJ2000'
                else:
                    ra_name = 'RAJ2000'
                    dec_name = 'DEJ2000'

                o.rename_column('_q', 'id')
                o.rename_column(ra_name, 'ra')
                o.rename_column(dec_name, 'dec')
                o = o.to_pandas()
                o = o.groupby('id')
                o = o.aggregate(np.nanmean)

                out[i] = o
        # print(out)
        # if only one catalog query is performed
        if len(out) == 1:
            # out = out[0]
            self.cache = out
            return out
        else:
            self.cache = out
            return out

    def __query_single__(self, coord):
        """
        Perform a query of a single catalog

        :param coord: The coordinates of the target(s)
        :type coord: astropy.coordinates.SkyCoord
        :return: The results of the query
        :rtype: astroquery.TableList
        """
        if type(self.radius) == tuple or type(self.radius) == list:
            rs = self.vizier.query_region(coord, width=self.radius[0], height=self.radius[1],
                                          catalog=self.catalog)
        else:
            rs = self.vizier.query_region(coord, self.radius,
                                          catalog=self.catalog)
        return rs

    def query_coordinate(self, coord):
        """
        Performs a query around coordinates

        :param coord: The central coordinates of the area
        :return:
        """
        # check if the coordinates are a tuple or a list
        if type(coord) == tuple or type(coord) == list:
            # if coord has 2 elements, take the first element as RA and the second as Dec
            # todo: add two additional parts with one for radius and one for box size
            if len(coord) == 2:
                ra = coord[0]
                dec = coord[1]
            # if coord has another length than 2, raise an error
            else:
                raise ValueError('If the coordinate are a tuple or list, ' +
                                 f'it must have 2 elements, not {len(coord)}!')
            coord = SkyCoord(ra*u.deg, dec*u.deg)
        elif type(coord) != SkyCoord:
            raise ValueError('"coord" must be a list or a tuple of two elements or a SkyCoord object, '
                             f'not a {type(coord)}')
        elif not coord.isscalar:
            raise ValueError('Only one coordinate set is allowed.')

        return self.__query_single__(coord)[0]

    def query(self, data, ra_name='ra', dec_name='dec', use_xmatch=False, blank=False):
        """
        Start a query around the input positions in the Vizier database

        :param data: The input data
        :type data: numpy.ndarray, astropy.table.Table
        :param ra_name: The name of the RA column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :param use_xmatch: True if xmatch should be used, even for small amount of data, else is False.
        :type use_xmatch: bool
        :param blank: True if the raw output should return, else False. Default is False.
        :type blank: bool
        :return: A list of tables with the query results
        :rtype: astropy.table.Table
        """
        if type(data) == DataFrame:
            data = Table.from_pandas(data)
        elif type(data) != Table:
            raise ValueError('Input must be a astropy Table or a pandas DataFrame.')
        # If there are more than 100 entries, use the xmatch interface to get the data
        # this will be faster
        try:
            if len(data) > 0 or use_xmatch:
                out = xmatch(data, ra_name, dec_name, self.name, blank=blank)

                return out
        except ConnectionError:
            if use_xmatch:
                raise ConnectionError('XMatch service does not respond.')
            print('No connection to XMatch, fall back to Vizier.')

        if type(data) == Table:
            # if no units are add to the coordinate columns
            if data[ra_name].unit is None:
                coordinates = SkyCoord(data[ra_name]*u.deg,
                                       data[dec_name]*u.deg)
            # if units are add to the coordinate columns
            else:
                coordinates = SkyCoord(data[ra_name],
                                       data[dec_name])
        # if the data are for example numpy.ndarray
        else:
            coordinates = SkyCoord(np.array(data[ra_name])*u.deg,
                                   np.array(data[dec_name])*u.deg)

        return self.__query__(coordinates)[0]

    def query_constrain(self, **constrains):
        """
        Do a constrain query (see :meth:`astroquery.vizier.Vizier.query_constraints` for details)

        :param constrains: The required constraints
        :type constrains: dict
        :return: The results of the query
        :rtype: astropy.table.Table
        """
        rs = self.vizier.query_constraints(catalog=self.catalog, **constrains)
        return rs[0]


class Gaia(Vizier):
    columns = ['RA_ICRS', 'DE_ICRS',
               'Plx', 'e_Plx',
               'pmRA', 'e_pmRA',
               'pmDE', 'e_pmDE',
               'Gmag', 'e_Gmag',
               'BPmag', 'e_BPmag',
               'RPmag', 'e_RPmag',
               'epsi', 'sepsi']
    name = names.GAIA
    ra_name = 'GAIA_RA_ICRS'
    dec_name = 'GAIA_DE_ICRS'

    def __init__(self):
        """
        Child class to query GAIA DR2 data
        """
        Vizier.__init__(self)


class TwoMass(Vizier):

    name = names.TWO_MASS

    def __init__(self):
        """
        Child class to query 2MASS data
        """
        Vizier.__init__(self)


class SkyMapper(Vizier):

    name = names.SkyMapper

    def __init__(self):
        """
        Child class to query 2MASS data
        """
        Vizier.__init__(self,
                        columns=[
                            '_RAJ2000', '_DEJ2000',
                            'uPetro', 'e_uPetro', 'vPetro', 'e_vPetro',
                            'gPetro', 'e_gPetro', 'rPetro', 'e_rPetro',
                            'iPetro', 'e_iPetro', 'zPetro', 'e_zPetro'
                        ]
                        )


class Wise(Vizier):
    name = names.WISE

    def __init__(self):
        """
        Child class to query AllWise data
        """
        Vizier.__init__(self)


class PanStarrs(Vizier):
    name = names.PANSTARRS

    def __init__(self):
        """
        Child class to query PanSTARRS DR1 data
        """
        Vizier.__init__(self)


class SDSS(Vizier):
    name = names.SDSS

    def __init__(self):
        """
        Child class to query SDSS DR12 data
        """
        Vizier.__init__(self)


class Apass(Vizier):
    name = names.APASS

    def __init__(self):
        """
        Child class to query APASS DR9 data
        """
        Vizier.__init__(self)


class Galex(Vizier):
    name = names.GALEX

    def __init__(self):
        """
        Child class to query GALEX DR5 data
        """
        Vizier.__init__(self)


class Kids(Vizier):

    name = names.KIDS

    def __init__(self):
        """
        Child class to query KiDS DR3 data
        """
        Vizier.__init__(self)


class Viking(Vizier):

    name = names.VIKING

    def __init__(self):
        """
        Child class to query VIKING DR2 data
        """
        Vizier.__init__(self)


class Ukidss(Vizier):

    name = names.UKIDSS

    def __init__(self):
        """
        Child class to query VIKING DR2 data
        """
        Vizier.__init__(self)


class BailerJones(Vizier):

    name = names.BAILER_JONES

    def __init__(self):
        """
        Child class to query Bailer-Jones distance estimations from the GAIA DR2
        """
        Vizier.__init__(self)


class MultiSurvey:
    surveys = []

    def __init__(self, sdss=True, panstarrs=True,
                 apass=True, two_mass=True,
                 all_wise=True, gaia=True,
                 galex=True):
        """
        This class provides a simple way to query multiple surveys at the same time.
        It also merges all the resulting data into one file with outer join (no-matches will get empty cells).

        :param sdss: True if the query should include SDSS DR12, else False
        :type sdss: bool
        :param panstarrs: True if the query should include PanSTARRS DR1, else False
        :type panstarrs: bool
        :param apass: True if the query should include APASS DR9, else False
        :type apass: bool
        :param two_mass: True if the query should include 2MASS, else False
        :type two_mass: bool
        :param all_wise: True if the query should include AllWISE, else False
        :type all_wise: bool
        :param gaia: True if the query should include GAIA DR2, else False
        :type gaia: bool
        :param galex: True if the query should include GALEX DR5, else False
        :type galex: bool
        """
        if sdss:
            self.surveys.append(SDSS())

        if panstarrs:
            self.surveys.append(PanStarrs())

        if apass:
            self.surveys.append(Apass())

        if two_mass:
            self.surveys.append(TwoMass())

        if all_wise:
            self.surveys.append(Wise())

        if gaia:
            self.surveys.append(Gaia())

        if galex:
            self.surveys.append(Galex())

    def query(self, data, ra_name='ra', dec_name='dec'):
        """
        Queries the chosen surveys

        :param data: The input data
        :type data: numpy.ndarray, astropy.table.Table
        :param ra_name: The name of the RA column
        :type ra_name: str
        :param dec_name: The name of the Dec column
        :type dec_name: str
        :return: A array with the results of the query
        :rtype: numpy.ndarray
        """
        th = []
        for s in self.surveys:
            t = Thread(target=s.query,
                       args=(data, ra_name, dec_name))
            t.start()
            th.append(t)

        # join the threads
        for t in th:
            t.join()

        data = data.copy()
        # match all the collected data in one table
        for s in self.surveys:
            if len(s.cache) > 0:
                data = match_catalogs(data, s.cache,
                                      ra_name, dec_name,
                                      s.ra_name, s.dec_name)
        return data


def get_survey(name):
    """
    Returns a survey object of the survey with the given name

    :param name: The name of the survey
    :type name: str
    :return: The survey object to perform queries
    :rtype: Vizier
    """
    name = name.lower()
    if name == 'sdss':
        survey = SDSS()
    elif name == 'panstarrs' or name == 'ps' or name == 'pan-starrs':
        survey = PanStarrs()
    elif name == '2mass':
        survey = TwoMass()
    elif name == 'gaia':
        survey = Gaia()
    elif name == 'wise':
        survey = Wise()
    elif name == 'galex':
        survey = Galex()
    elif name == 'apass':
        survey = Apass()
    elif name == 'ukidss':
        survey = Ukidss()
    elif name == 'kids':
        survey = Kids()
    elif name == 'viking':
        survey = Viking()
    elif name == 'bailer-jones' or name == 'bj':
        survey = BailerJones()
    elif name == 'skymapper':
        survey = SkyMapper()
    else:
        raise AttributeError(f'No survey with name {name} known!')
    return survey


def query_by_name(name, data, ra_name='ra', dec_name='dec'):
    """
    Queries a survey by its name

    :param name: Name of the survey
    :type name: str
    :param data: data with the coordinates, it must have at least two columns with 'ra_name' and 'dec_name'
    :type data: astropy.table.Table, pandas.DataFrame
    :param ra_name: Name of the RA column
    :type ra_name: str
    :param dec_name: Name of the Dec column
    :type dec_name: str
    :return: The data from the survey
    :rtype: numpy.ndarray
    """
    survey = get_survey(name)
    out = survey.query(data, ra_name, dec_name)
    return out


def constrain_query(name, **kwargs):
    """
    Constrain query for large surveys by their names

    .. code-block:: python

        results = constrain_query('gaia', Gmag='>16 & <18')

    :param name: The name of the survey
    :type name: str
    :param kwargs: Constrains
    :type kwargs: dict
    :return: The results of the query
    :rtype: astropy.table.Table
    """
    survey = get_survey(name)
    return survey.query_constrain(**kwargs)


def query_simbad(coordinates):
    """
    Queries the Simbad database for potential entries/classifications

    :param coordinates: The coordinates
    :type coordinates: Phosphorpy.data.sub.coordinates.CoordinateTable
    :return:
    """
    sc = coordinates.to_sky_coord()
    rs_tot = []

    # download blocks of 2000 coordinates to avoid timeout errors because of the server limitation
    steps = 2000
    i = 0
    while i*steps < len(sc):
        rs = Simbad.query_region(sc[i*steps: (i+1)*steps], radius=2*u.arcsec)
        rs['_q'] = i*steps+rs['_q']
        rs_tot.append(rs)
        i += 1
    rs_tot = vstack(rs_tot)

    if rs_tot is None:
        return rs_tot
    coord = []
    ra_str = '{}h{}m{}s'
    dec_str = ' {}d{}m{}s'
    for ra, dec in zip(rs_tot['RA'], rs_tot['DEC']):
        c = ra_str.format(*(ra.split(' ')))
        c += dec_str.format(*(dec.split(' ')))
        coord.append(c)

    s = SkyCoord(np.array(coord))
    rs_tot['ra'] = s.ra
    rs_tot['dec'] = s.dec

    rs = coordinates.match(rs_tot)

    return rs
