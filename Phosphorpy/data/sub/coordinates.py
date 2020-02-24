import numba as nb
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units.quantity import Quantity
from pandas import DataFrame
from sklearn.neighbors import NearestNeighbors

from .plots.coordinates import CoordinatePlot
from .table import DataTable

try:
    from Phosphorpy.data.sub.interactive_plotting.coordinates import CoordinatePlot as CoordinatePlotHV
except ImportError:
    CoordinatePlotHV = None

RA_NAMES = ['ra', 'Ra', 'RA', 'RAJ2000', 'RA_ICRS', '_RAJ2000']
DEC_NAMES = ['dec', 'Dec', 'DEC', 'DEJ2000', 'DECJ2000', 'DE_ICRS', '_DEJ2000']
L_NAMES = ['l', 'L']
B_NAMES = ['b', 'B']


class CoordinateTable(DataTable):
    __slots__ = ('_data', '_head', '_plot')

    def __init__(self, data, head=None, mask=None):
        """
        Table to handle coordinates

        :param data:
        :param head:
        """
        DataTable.__init__(self, mask=mask)
        self._head = head
        self._plot = CoordinatePlot(self)

        if CoordinatePlotHV.holoviews():
            self._hv_plot = CoordinatePlotHV(self)

        # create SkyCoord objects to get the galactic coordinates of the sources
        if type(data) == SkyCoord:
            ra = np.array(data.ra.degree)
            dec = np.array(data.dec.degree)
            lon = np.array(data.galactic.l)
            b = np.array(data.galactic.b)
        elif type(data) == DataFrame:
            s = SkyCoord(
                data['ra'].values*u.deg,
                data['dec'].values*u.deg
            )
            ra = s.ra.degree
            dec = s.dec.degree
            lon = s.galactic.l.degree
            b = s.galactic.b.degree
        else:
            if type(data) != np.ndarray or len(data.shape) == 1:
                columns = get_column_names(data)

                s = SkyCoord(np.array(data[find_ra_column(columns)]) * u.deg,
                             np.array(data[find_dec_column(columns)]) * u.deg)
            else:
                s = SkyCoord(data[:, 0] * u.deg, data[:, 1] * u.deg)
            ra = np.array(s.ra.degree)
            dec = np.array(s.dec.degree)
            lon = np.array(s.galactic.l)
            b = np.array(s.galactic.b)
        d = np.zeros((len(data), 4))
        d[:, 0] = ra
        d[:, 1] = dec
        d[:, 2] = lon
        d[:, 3] = b
        self._data = DataFrame(data=d, columns=['ra', 'dec', 'l', 'b'],
                               index=np.arange(len(ra))+1)
        self._data['row_id'] = np.arange(len(self._data))+1
        self._data = self._data.set_index('row_id')

    def __getitem__(self, item):
        if type(item) == str:
            # check if the item is one of the default ones
            if item in RA_NAMES:
                return self._data['ra'].values
            elif item in DEC_NAMES:
                return self._data['dec'].values
            elif item in L_NAMES:
                return self._data['l'].values
            elif item in B_NAMES:
                return self._data['b'].values
            else:
                raise KeyError('Key {} not known! Choose one of the default ones.'.format(item))
        elif type(item) == int:
            return self._data.loc[item]

    def __eq__(self, other):
        if type(other) == SkyCoord or type(other) == np.ndarray or type(other) == CoordinateTable:
            if len(self) == len(other):
                # of the matching doesn't include a -1, which means that there is no match,
                # the coordinate sets are equal
                if len(self.match(other)) == len(self):
                    return True
                else:
                    return False
            else:
                return False
        else:
            raise ValueError('Allowed types are astropy.coordinate.SkyCoord, numpy.ndarray '
                             'or SearchEngine.data.sub.coordinates.CoordinateTable, not {}'.format(type(other)))

    def __len__(self):
        return len(self._data)

    def to_table(self, galactic=False):
        """
        Returns the coordinates as pandas.DataFrame

        :param galactic: True if the galactic coordinates should be included, else False. Default is False.
        :type galactic: bool
        :return: The coordinates in a DataFrame
        :rtype: pandas.DataFrame
        """
        if galactic:
            return self.data
        else:
            return self.data[['ra', 'dec']]

    def to_sky_coord(self, source_id=-1):
        """
        Return the coordinate(s) back as a SkyCoord object

        :param source_id:
            The id of the source. Default is -1, which means that all coordinates are converted to SkyCoord.
        :type source_id: int, list, tuple
        :return: The SkyCoord object of the source(s)
        :rtype: astropy.coordinates.SkyCoord
        """
        if source_id == -1:
            s = SkyCoord(self.data['ra'].values * u.deg, self.data['dec'].values * u.deg)
        else:
            s = SkyCoord(self.data['ra'].values[source_id] * u.deg,
                         self.data['dec'].values[source_id] * u.deg)
        return s
    
    def to_astropy_table(self, category='coordinates', **kwargs):
        return super(CoordinateTable, self).to_astropy_table(category, **kwargs)

    def _combine_coordinates(self, other):

        # prepare the input coordinates for the match
        if type(other) == SkyCoord:
            d = np.zeros((len(other), 2))
            d[:, 0] = other.ra.degree
            d[:, 1] = other.dec.degree
        elif type(other) == CoordinateTable or type(other) == Table or type(other) == DataFrame:
            d = np.zeros((len(other), 2))
            d[:, 0] = other['ra']
            d[:, 1] = other['dec']
        elif type(other) == np.ndarray:
            d = other.copy()
        else:
            raise ValueError()

        x = np.concatenate((self._data[['ra', 'dec']].values, d), axis=0)
        return x

    def match(self, other, radius=2 * u.arcsec,
              ra_name='ra', dec_name='dec'):
        """
        Matches these coordinates with another coordinate set. To do that it uses a nearest neighbor algorithm
        to find the next neighbor.
        The distance is computed with the approximation of small angles

        .. math::

            d = \sqrt{\left(\Delta\\alpha \cdot \cos\left(\delta\\right)\\right)^2 + \left(\Delta \delta\\right)^2}

        as an additional approximation because performance reasons, in this algorithm its assumed that

        .. math::

            \cos \delta_1 \\approx \cos \delta_2

        and therefore we can use

        .. math::

            \Delta \\alpha \cdot \cos \delta = \\alpha_1 \cdot \cos \delta - \\alpha_2 \cdot \cos \delta \\\\
            \\approx \\alpha_1 \cos\delta_1 - \\alpha_2 \cos \delta_2

        :param other: The second coordinate set to match with it
        :type other: numpy.ndarray, astropy.coordinates.SkyCoord, Phosphorpy.data.sub.coordinates.CoordinateTable
        :param radius: The cross-match radius as float or astropy.unit. If it is a float, it will be taken as degree.
        :type radius: float, astropy.units
        :param ra_name: Name of the RA column. Default is 'ra'.
        :type ra_name: str
        :param dec_name: Name of the Dec column. Default is 'dec'.
        :type dec_name: str
        :returns: An array with the indices of the match sources. Sources without a match will have -1.
        :rtype: pandas.DataFrame
        """
        x = self._data[['ra', 'dec']].values
        # convert the radius to degree and to a float
        if type(radius) == Quantity:
            radius = radius.to(u.deg)

        if type(other) != DataFrame:
            if type(other) == Table:
                other = other.to_pandas()
            elif type(other) == np.ndarray:
                other = DataFrame(other, columns=['ra', 'dec'])
            elif type(other) == CoordinateTable:
                other = other.data
            elif type(other) == SkyCoord:
                sk = self.to_sky_coord()
                index, dis2d, dis3d = sk.match_to_catalog_sky(other)
                return other[index[dis2d < radius]]
            else:
                raise ValueError('Unsupported input type ({})'.format(str(type(other))))

        # approximate small distance
        x[:, 0] *= np.cos(np.deg2rad(x[:, 1]))
        nn = NearestNeighbors(n_neighbors=1).fit(x)

        y = other[[ra_name, dec_name]].values
        y[:, 0] *= np.cos(np.deg2rad(y[:, 1]))
        distances, indices = nn.kneighbors(y)
        m = distances[:, 0] <= radius.value

        new_indices = self.data.index.values[indices[:, 0][m]]

        # select the matched sources only and sort them in the same order
        other = other[m]
        other['row_id'] = new_indices
        other = other.set_index('row_id')
        return other

    def get_closest_source_id(self, coordinate=None, ra=None, dec=None):
        """
        Returns the ID of the closest sources and the distance to source.

        coordinate or ra and dec has to be given. If all three are given,
        the coordinate object is taken.

        :param coordinate: Coordinates (or close to) of the required ID.
        :type coordinate: SkyCoord
        :param ra: The RA component of the coordinates
        :type ra: float
        :param dec: The Dec component of the coordinates
        :type dec: float
        :return: The closest ID and the distance to the corresponding source
        :rtype: int, u.Quantity
        """
        if coordinate is None:
            if ra is None or dec is None:
                raise ValueError('If no coordinates are given, RA and Dec must be given.')
            coordinate = SkyCoord(ra*u.deg, dec*u.deg)
        distances = coordinate.separation(self.to_sky_coord())
        min_id = np.argmin(distances)
        return min_id, distances[min_id]

    def write(self, path, data_format='parquet'):
        data = self.to_table()
        if data_format == 'parquet':
            data.to_parquet(path)
        elif data_format == 'csv':
            data.to_csv(path)
        elif data_format == 'latex':
            data.to_latex(path)
        elif data_format == 'fits':
            Table.from_pandas(data).write(path, overwrite=True)


def get_column_names(d):
    """
    Returns the name of the columns
    :param d: The input data set
    :type d: numpy.ndarray, astropy.table.Table, pandas.DataFrame
    :return: A list with the column names
    :rtype: list
    """
    if type(d) == np.ndarray:
        return d.dtype.names
    elif type(d) == Table:
        return d.colnames
    elif type(d) == DataFrame:
        return d.columns


def find_ra_column(columns):
    """
    Search in the columns for a possible RA name

    :param columns: The column names of the input data set
    :type columns: list
    :return: The name of the RA column, if it is in. If the name is not included it will raise a KeyError
    :rtype: str
    """
    for r in RA_NAMES:
        if r in columns:
            return r
    raise KeyError('RA column not found!')


def find_dec_column(columns):
    """
    Search in the columns for a possible Dec name

    :param columns: The column names of the input data set
    :type columns: list
    :return: The name of the Dec column, if it is in. If the name is not included it will raise a KeyError
    :rtype: str
    """
    for d in DEC_NAMES:
        if d in columns:
            return d
    raise KeyError('Dec column not found!')


@nb.jit('f8(f8[:], f8[:])', nopython=True)
def spherical_distance(x1, x2):
    """
    Computes the spherical distance in the approximation of small distances

    .. math::

        d = \sqrt{x1\cdot \cos\left(x2\\right) + x2}

    :param x1: The first vector
    :type x1: Union
    :param x2: The second vector
    :type x2: Union
    :return: THe distance between x1 and x2
    :rtype: float
    """
    x = x1 - x2
    return np.sqrt((x[0] * np.cos(np.deg2rad(x1[1]))) ** 2 + (x[1]) ** 2)
