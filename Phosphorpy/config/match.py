import warnings

import numpy as np
from astropy import units as u
from astropy.table import Table, join, vstack
from astropy.units.quantity import Quantity
from sklearn.cluster import DBSCAN
import pandas

warnings.simplefilter('ignore')


def fill_missing_labels(d, max_label):
    """
    Fills the labels of unmatched sources with unique numbers higher than the highest current label

    :param d: The input data with at least the column 'label'
    :type d: astropy.table.Table
    :param max_label: The current lowest not used label number
    :type max_label: int
    :return: The input data with the new labels and the new lowest not used label number
    :rtype: astropy.table.Table, int
    """
    if 'label' not in d.colnames:
        raise ValueError('Missing column "label"!')
    mask = d['label'] == -1
    unidentified = len(d[mask])
    d['label'][mask] = np.linspace(max_label, max_label + unidentified,
                                   unidentified, dtype=np.int32)
    return d, d['label'].max() + 1


def convert_input_data(d):
    """
    Checks if the input data are a astropy.table.Table and if not it converts it to a Table.

    :param d: The input data
    :type d: numpy.ndarray, astropy.table.Table
    :return: The input data as an astropy.table.Table
    :rtype: astropy.table.Table
    """
    # TODO: reduce to coordinates and id
    if type(d) is not Table:
        if type(d) is pandas.DataFrame:
            d = Table.from_pandas(d)
        else:
            d = Table(d)
    return d


def match_catalogs(d1, d2, ra1, dec1, ra2, dec2, join_type='outer',
                   match_radius=2 * u.arcsec):
    """
    Matches two catalog into one
    :param d1: The first input catalog
    :type d1: numpy.ndarray, astropy.table.Table
    :param d2: The second input catalog
    :type d2: numpy.ndarray, astropy.table.Table
    :param ra1: The name of the RA column of the first input catalog
    :type ra1: str
    :param dec1: The name of the Dec column of the first input catalog
    :type dec1: str
    :param ra2: The name of the RA column of the second input catalog
    :type ra2: str
    :param dec2: The name of the Dec column of the second input catalog
    :type dec2: str
    :param join_type: The join mode (inner or outer), default is outer
    :type join_type: str
    :param match_radius: The matching radius in which two sources are considered as the same
    :type match_radius: int, float, astropy.units.quantity.Quantity
    :return: The matched catalog
    :rtype: numpy.ndarray
    """
    # if the match radius is not a Quantity, set the unit to arcsec
    if type(match_radius) is not Quantity:
        match_radius = match_radius * u.arcsec

    x = np.zeros((len(d1) + len(d2), 2))
    x[:len(d1), 0] = d1[ra1]
    x[:len(d1), 1] = d1[dec1]
    x[len(d1):, 0] = d2[ra2]
    x[len(d1):, 1] = d2[dec2]

    db = DBSCAN(eps=match_radius.to(u.deg).value,
                min_samples=2)
    db.fit(x)

    d1 = convert_input_data(d1)
    d2 = convert_input_data(d2)

    d1['label'] = db.labels_[:len(d1)]
    d2['label'] = db.labels_[len(d1):]

    max_label = np.max(db.labels_) + 1

    d1, max_label = fill_missing_labels(d1, max_label)
    d2, max_label = fill_missing_labels(d2, max_label)

    d = join(d1, d2, keys='label', join_type=join_type)
    d = np.array(d)
    return d


def group_by_coordinates(d, ra_name, dec_name):
    """
    Groups the input catalog sources by their coordinates in a circle of 2 arcsec and
    perform an nanmean (mean-function with ignoring nan values).

    :param d: The input catalog
    :type d: astropy.table.Table
    :param ra_name: The RA column name
    :type ra_name: str
    :param dec_name: The Dec column name
    :type dec_name: str
    :return: The mean values of sources with multiple entries or the values itself if there is one entry only
    :rtype: astropy.table.Table
    """
    print(d.colnames)
    # exclude targets with NAN values in the coordinates (in the stupid way)
    d = d[(d[ra_name] > -99) & (d[dec_name] > -99)]

    x = np.zeros((len(d), 2))
    x[:, 0] = d[ra_name]
    x[:, 1] = d[dec_name]
    db = DBSCAN(eps=3. / 3600,
                min_samples=2)
    db.fit(x)
    d['coord_label'] = db.labels_

    # make a mask of sources which have one entry only
    mask = d['coord_label'] == -1
    out = [d[mask]]

    # select sources with multiple entries
    d = d[mask == False]

    if len(d) != 0:
        # group the sources with multiple entries by their coordinate label and compute the nanmean value
        d = d.group_by('coord_label')
        d = d.groups.aggregate(np.nanmean)

        # append the results to the output list and merge both tables
        out.append(d)
        out = vstack(out)
    else:
        out = out[0]

    # delete the unused column
    del out['coord_label']
    return out
