from astropy.coordinates import SkyCoord
from astropy.units.quantity import Quantity
from astropy import units as u
from astropy.table import Table, join
from sklearn.cluster import DBSCAN
from pandas import DataFrame
import numpy as np


from .table import DataTable
from .plots.coordinates import CoordinatePlot

RA_NAMES = ['ra', 'Ra', 'RA', 'RAJ2000', 'RA_ICRS', '_RAJ2000']
DEC_NAMES = ['dec', 'Dec', 'DEC', 'DEJ2000', 'DECJ2000', 'DE_ICRS', '_DEJ2000']
L_NAMES = ['l', 'L']
B_NAMES = ['b', 'B']


class CoordinateTable(DataTable):

    def __init__(self, data, head=None, mask=None):
        """
        Table to handle coordinates

        :param data:
        :param head:
        """
        DataTable.__init__(self, mask=mask)
        self._head = head
        self._plot = CoordinatePlot(self)

        # create SkyCoord objects to get the galactic coordinates of the sources
        if type(data) == SkyCoord:
            ra = np.array(data.ra.degree)
            dec = np.array(data.dec.degree)
            lon = np.array(data.galactic.l)
            b = np.array(data.galactic.b)
        else:
            if type(data) != np.ndarray or len(data.shape) == 1:
                columns = get_column_names(data)

                # todo: check if the data is an astropy.table.Table and if it has units
                s = SkyCoord(np.array(data[find_ra_column(columns)])*u.deg,
                             np.array(data[find_dec_column(columns)])*u.deg)
            else:
                s = SkyCoord(data[:, 0]*u.deg, data[:, 1]*u.deg)
            ra = np.array(s.ra.degree)
            dec = np.array(s.dec.degree)
            lon = np.array(s.galactic.l)
            b = np.array(s.galactic.b)
        d = np.zeros((len(data), 4))
        d[:, 0] = ra
        d[:, 1] = dec
        d[:, 2] = lon
        d[:, 3] = b
        self._data = d

    def __getitem__(self, item):
        if type(item) == str:
            # check if the item is one of the default ones
            if item in RA_NAMES:
                return self._data[:, 0]
            elif item in DEC_NAMES:
                return self._data[:, 1]
            elif item in L_NAMES:
                return self._data[:, 2]
            elif item in B_NAMES:
                return self._data[:, 3]
            else:
                raise KeyError('Key {} not known! Choose one of the default ones.'.format(item))

    def __eq__(self, other):
        if type(other) == SkyCoord or type(other) == np.ndarray or type(other) == CoordinateTable:
            if len(self) == len(other):
                # of the matching doesn't include a -1, which means that there is no match,
                # the coordinate sets are equal
                if -1 not in self.match(other):
                    return True
                else:
                    return False
            else:
                return False
        else:
            raise ValueError('Allowed types are astropy.coordinate.SkyCoord, numpy.ndarray '
                             'or SearchEngine.data.sub.coordinates.CoordinateTable, not {}'.format(type(other)))

    def __len__(self):
        return len(self._data[:, 0])

    def to_table(self, galactic=False):
        """
        Returns the coordinates as pandas.DataFrame

        :param galactic: True if the galactic coordinates should be included, else False. Default is False.
        :type galactic: bool
        :return: The coordinates in a DataFrame
        :rtype: pandas.DataFrame
        """
        if galactic:
            return DataFrame(self.data, columns=['ra', 'dec', 'l', 'b'])
            pass
        else:
            return DataFrame(self.data[:, :2], columns=['ra', 'dec'])

    def as_sky_coord(self, source_id=-1):
        """
        Return the coordinate(s) back as a SkyCoord object

        :param source_id:
            The id of the source. Default is -1 which means that all coordinates are converted to SkyCoord.
        :type source_id: int, list, tuple
        :return: The SkyCoord object of the source(s)
        :rtype: astropy.coordinates.SkyCoord
        """
        if source_id == -1:
            s = SkyCoord(self.data[:, 0]*u.deg, self.data[:, 1]*u.deg)
        else:
            s = SkyCoord(self.data[source_id, 0]*u.deg, self.data[source_id, 1]*u.deg)
        return s

    def match(self, other, radius=1*u.arcsec):
        """
        Matches these coordinates with another coordinate set

        :param other: The second coordinate set to match with it
        :type other: numpy.ndarray, astropy.coordinates.SkyCoord, SearchEngine.data.sub.coordinates.CoordinateTable
        :param radius: The cross-match radius as float or astropy.unit. If it is a float, it will be taken as degree.
        :type radius: float, astropy.units
        :returns: An array with the indices of the match sources. Sources without a match will have -1.
        :rtype: numpy.ndarray
        """
        if type(other) == SkyCoord or type(other) == np.ndarray or type(other) == CoordinateTable:
            # convert the radius to degree and to a float
            if type(radius) == Quantity:
                radius = radius.to(u.deg).value

            # prepare the input coordinates for the match
            if type(other) == SkyCoord:
                d = np.zeros((len(other), 2))
                d[:, 0] = other.ra.degree
                d[:, 1] = other.dec.degree
            elif type(other) == CoordinateTable:
                d = np.zeros((len(other), 2))
                d[:, 0] = other['ra']
                d[:, 1] = other['dec']
            elif type(other) == np.ndarray:
                d = other.copy()
            else:
                raise ValueError()

            x = np.concatenate((self._data[:, :2], d), axis=0)

            # do the cross-match with a DBSCAN
            db = DBSCAN(eps=radius, min_samples=2)
            db.fit(x)

            # todo: replace this algorithm with a more efficient one
            next_label = np.max(db.labels_)+1
            data_label = Table()
            data_label['labels'] = db.labels_[:len(self)]
            data_label['id_1'] = np.arange(0, len(self), dtype=np.int32)
            m = data_label['labels'] == -1
            data_label['labels'][m] = np.arange(next_label, len(data_label[m])+next_label,
                                                1, dtype=np.int32)

            next_label = max(next_label, np.max(data_label['labels']))
            input_label = Table()
            input_label['labels'] = db.labels_[len(self):]

            input_label['id_2'] = np.arange(0, len(d), dtype=np.int32)
            m = input_label['labels'] == -1
            input_label['labels'][m] = np.arange(next_label, len(input_label[m])+next_label,
                                                 1, dtype=np.int32)

            output = join(data_label, input_label, keys='labels', join_type='outer')
            try:
                output['id_1'][output['id_1'].mask] = -1
            except AttributeError:
                pass
            return np.array(output['id_2'])
        else:
            raise ValueError('Allowed types are astropy.coordinate.SkyCoord, numpy.ndarray '
                             'or SearchEngine.data.sub.coordinates.CoordinateTable, not {}'.format(type(other)))

    def write(self, path, data_format='parquet'):
        data = self.to_table()
        if data_format == 'parquet':
            data.to_parquet(path)
        elif data_format == 'csv':
            data.to_csv(path)
        elif data_format == 'sql':
            # todo: implement engine for sql writing
            data.to_sql(path, None)
        elif data_format == 'latex':
            data.to_latex(path)
        elif data_format == 'fits':
            Table.from_pandas(self.data).write(path)


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
