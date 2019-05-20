import numpy as np
import pylab as pl
from astropy.table import Table
from astroquery.vizier import Vizier


class SdssSpec:
    sdss = Vizier(columns=['RA_ICRS', 'DE_ICRS',
                           'umag', 'gmag', 'rmag', 'imag', 'zmag',
                           'subCl'],
                  column_filters={"SpObjID": ">0"})
    sources = None
    stats = None
    m_stats = None

    def __init__(self, row_limit=100000, path=None):
        """
        Class to download and handle SDSS spectroscopic identificated sources

        :param row_limit: Number sources to download
        :type row_limit: int
        :param path: Path to a file with the data from SDSS DR12. Default is None.
        :type path: str
        """
        self.sdss.ROW_LIMIT = row_limit

        if path is not None:
            self.sources = Table.read(path)

    def download_classification(self, cache=True, r_min=20):
        """
        Downloads the classifications and magnitudes of sources from the SDSS DR12 catalog

        :param cache:
            If a path was set, cache=True means that it will use the data from the file.
            False means that it will download a new set. Default is True.
        :type cache: bool
        :param r_min: Minimal brightness in the r-band
        :type r_min: float
        :return: Table with the downloaded data
        :rtype: astropy.table.Table
        """
        try:
            if cache and self.sources is not None:
                return self.sources
            # query the SDSS DR12 for objects with a spectral identification
            # and a set sub-class identification with the
            # astroquery vizier class
            self.sources = self.sdss.query_constraints(catalog='V/147/sdss12',
                                                       SpObjID='>0',
                                                       subCl='!= ',
                                                       rmag='<{}'.format(r_min))[0]
            self.sources = __sub_class_rename__(self.sources)
            self.__stats__()
            return self.sources
        except IndexError:
            return None

    def __stats__(self):
        """
        Counts how often a spectral type appears in the data set

        :return: The counts of the spectral type
        :type: astropy.table.Table
        """
        sub_classes, sub_class_counts = np.unique(self.sources['subCl'], return_counts=True)
        self.stats = Table()
        self.stats['subCl'] = sub_classes
        self.stats['count'] = sub_class_counts
        self.__main_identification_stats__()
        return self.stats

    def __main_identification_stats__(self):
        """
        Returns the main identification of the sources and how often they are in the sample.
        A main identification is for a star the main spectral type.
        Galaxies and QSO's are compressed to only one type (AGN and STARFORMING).
        :return: Dict with the main identifications as keys and the values are how often they are in the sample.
        :rtype: dict
        """
        m_stats = {}
        for s in self.stats:
            sub_cl = s['subCl']

            # if the length is two, it is likely a star identification
            if len(sub_cl) == 2:
                if sub_cl[-1].isdigit():
                    # take only the spectral identification key
                    sub_cl = sub_cl[0]
                    m_stats = __add_counts__(m_stats, sub_cl, s)
                # if the second key is not a number
                else:
                    if 'OB' in sub_cl:
                        sub_cl = 'O'
                    m_stats = __add_counts__(m_stats, sub_cl, s)
            # if the length is larger than two
            else:
                # Compress galaxy and QSO's types
                if sub_cl == 'BROADLINE':
                    sub_cl = "AGN"
                elif sub_cl == 'STARBURST':
                    sub_cl = 'GALAXY'
                elif sub_cl == 'STARFORMING':
                    sub_cl = 'GALAXY'
                elif sub_cl == 'Ldwarf':
                    sub_cl = 'L'
                elif sub_cl == 'Carbon_lines':
                    sub_cl = 'Carbon'
                m_stats = __add_counts__(m_stats, sub_cl, s)
        self.m_stats = m_stats
        return m_stats

    def show_sub_classes(self):
        """
        Plots the distribution of different main spectral identifications
        as a bar plot

        :return:
        """
        pl.clf()
        sp = pl.subplot()
        x = np.zeros((len(self.m_stats), 2))
        for i, k in enumerate(self.m_stats):
            x[i, 0] = i
            x[i, 1] = self.m_stats[k]
        sp.bar(x[:, 0], x[:, 1])
        pl.xticks(x[:, 0], list(self.m_stats.keys()),
                  rotation='vertical')
        sp.set_ylabel('counts')
        sp.set_xlabel('main subclass type')
        sp.set_yscale('log')
        pl.show()


def __add_counts__(m_stats, name, s):
    if name in m_stats.keys():
        m_stats[name] += s['count']
    else:
        m_stats[name] = s['count']
    return m_stats


def __sub_class_rename__(tab):
    """
    Changes the sub-class identification of SDSS DR12 in such a way, that for stars only the main identification
    (Kx, Fy, ...) is left. WD's are also only one sub-class and it won't be distingushed between
    hot and cool ones.
    The identification of galaxies and QSO's is untouched.

    :param tab: The input table with the data from the SDSS DR12
    :type tab: numpy.ndarray, astropy.table.Table
    :return: The input data with the changed sub-class identification
    :rtype: numpy.ndarray, astropy.table.Table
    """
    for i, t in enumerate(tab):
        sub_cl = t['subCl'].split(' ')

        # if there are additional identification numbers
        if len(sub_cl) > 1:
            tab['subCl'][i] = sub_cl[0]

        sub_cl = sub_cl[0]

        # if there is sub spectral identifications. For example K5V change it to K5
        if (len(sub_cl) == 3 or len(sub_cl) == 4) and sub_cl[-1] != 'N':
            if 'sd' in sub_cl:
                sub_cl = sub_cl[2:]
            else:
                sub_cl = sub_cl[:2]
            tab['subCl'][i] = sub_cl
        # if the sub-class is not clear, use just the first one
        elif '/' in sub_cl:
            if 'Ib' in sub_cl:
                sub_cl = sub_cl.split('Ib')[0]
            else:
                sub_cl = sub_cl.split('/')[0]
            tab['subCl'][i] = sub_cl
        # if WD is part of the identification
        elif 'WD' in sub_cl:
            tab['subCl'][i] = 'WD'
        elif 'III' in sub_cl:
            sub_cl = sub_cl[:2]
            tab['subCl'][i] = sub_cl
        elif ':' in sub_cl:
            sub_cl = sub_cl.split(':')[-1]
            tab['subCl'][i] = sub_cl
        elif '.' in sub_cl:
            if 'Ib' in sub_cl:
                sub_cl = sub_cl.split('Ib')[0]
            elif 'V.' in sub_cl:
                sub_cl = sub_cl.split('V.')[0]
            sub_cl = sub_cl.split('.')[0]
            tab['subCl'][i] = sub_cl
        elif 'lb' in sub_cl:
            sub_cl = sub_cl.split('lb')[0]
            tab['subCl'][i] = sub_cl
        elif 'Ib' in sub_cl:
            sub_cl = sub_cl.split('Ib')[0]
            tab['subCl'][i] = sub_cl
        elif 'IV' in sub_cl:
            sub_cl = sub_cl.split('IV')[0]
            tab['subCl'][i] = sub_cl
        elif 'Ia' in sub_cl:
            sub_cl = sub_cl.split('Ia')[0]
            tab['subCl'][i] = sub_cl
        elif 'Vn' in sub_cl:
            sub_cl = sub_cl.split('Vn')[0]
            tab['subCl'][i] = sub_cl

    tab = tab[tab['subCl'] != '']
    return tab
