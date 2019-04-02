from Phosphorpy.external.vizier import Gaia, BailerJones
from Phosphorpy.data.sub.plots.astrometry import AstrometryPlot
from astropy.table import vstack
from .table import DataTable
import numpy as np
import pandas as pd


def _only_nearest(data):

        row_ids, row_id_count = np.unique(data['row_id'], return_counts=True)
        # find or multiple detections the closest one
        nearest = []
        for rid in row_ids[row_id_count > 1]:
            g = data[data['row_id'] == rid]
            nearest.append(g[g['angDist'] == np.min(g['angDist'])])
        nearest = vstack(nearest)
        data = vstack([nearest, data])
        return data


def _download_gaia_data(coordinates):
    """
    Downloads the Gaia data for the given coordinates

    :param coordinates: The coordinates of the targets for which Gaia data are required
    :type coordinates: Phosphorpy.data.sub.coordinates.CoordinateTable
    :return: The Gaia data
    :rtype:
    """
    g = Gaia()

    gaia = g.query(coordinates.to_table(), 'ra', 'dec', use_xmatch=True, blank=True)

    # find or multiple detections the closest one
    # gaia = _only_nearest(gaia)

    gaia_coords = gaia[['row_id', 'ra', 'dec']]
    gaia = gaia.to_pandas()
    gaia = gaia.drop_duplicates('row_id')
    gaia = gaia.set_index('row_id')

    # join the downloaded gaia data to create empty lines if gaia doesn't provide data for a specific object
    gaia = gaia[['ra', 'ra_error', 'dec', 'dec_error',
                 'parallax', 'parallax_error',
                 'pmra', 'pmra_error',
                 'pmdec', 'pmdec_error']]
    return gaia, gaia_coords


def _download_bailer_jones_data(gaia_coords):
    """
    Downloads the distances estimated by Bailer-Jones

    :param gaia_coords: The gaia coordinates of the required targets
    :return:
    """
    # download Bailer-Jones distance estimations
    bj = BailerJones()
    bj = bj.query(gaia_coords, 'ra', 'dec', use_xmatch=True, blank=True)

    bj = bj[['row_id', 'rest', 'b_rest', 'B_rest', 'rlen', 'ResFlag', 'ModFlag']].to_pandas()
    bj = bj.drop_duplicates('row_id')
    bj = bj.set_index('row_id')
    return bj


class AstrometryTable(DataTable):

    def __init__(self, mask):
        """
        AstrometryTable in the class to handle GAIA astrometry data, which includes
        coordinates, proper motion, parallax and its errors.

        :param mask: The mask object for this data
        :type mask: Phosphorpy.data.sub.table.Mask
        """
        DataTable.__init__(self, mask=mask)
        self._plot = AstrometryPlot(self)

    @staticmethod
    def load_to_dataset(ds):
        g = Gaia()
        gaia = g.query(ds.coordinates.as_sky_coord(), 'ra', 'dec', use_xmatch=True, blank=True)
        astronomy = AstrometryTable(ds.mask)
        astronomy._data = gaia[['ra', 'ra_error', 'dec', 'dec_error',
                                'parallax', 'parallax_error',
                                'pm_ra', 'pm_ra_error',
                                'pm_dec', 'pm_dec_error']]
        # bj = BailerJones()
        # bj = bj.query(ds.coordinates.as_sky_coord(), 'ra', 'dec', use_xmatch=True, blank=True)
        ds.astrometry = astronomy

    @staticmethod
    def load_astrometry(coordinates):
        """
        Loads GAIA astrometry from the server.

        :param coordinates: The coordinate table with the data.
        :type coordinates: Phosphorpy.data.sub.coordinatesCoordinateTable
        :return: An AstrometryTable with the Gaia astrometry (ra, dec, parallax, pmra, pmdec and their errors)
        :rtype: AstrometryTable
        """
        gaia, gaia_coords = _download_gaia_data(coordinates)

        # download Bailer-Jones distance estimations
        bj = _download_bailer_jones_data(gaia_coords)

        # join the original gaia data with the Bailer-Jones distance data
        gaia = gaia.join(bj, how='outer')

        df = pd.DataFrame()
        df['row_id'] = np.linspace(1, len(coordinates.data), len(coordinates.data),
                                   dtype=np.int32)
        df = df.set_index('row_id')
        gaia = df.join(gaia, how='outer')

        astronomy = AstrometryTable(coordinates.mask)

        astronomy._data = gaia
        return astronomy

    def proper_motion(self, cos_correction=False):
        """
        Return the proper motion and the errors of it.

        :param cos_correction:
            True if the RA components should be corrected by the declination, else False.
            Default is False
        :type cos_correction: bool
        :return: pmra, pmdec, pmra_error, pmdec_error
        :rtype: numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
        """

        x = self._data['pmra'].values
        y = self._data['pmdec'].values
        x_err = self._data['pmra_error'].values
        y_err = self._data['pmdec_error'].values
        if cos_correction:
            dec_rad = np.deg2rad(self._data['dec'].values)
            cos_c = np.cos(dec_rad)
            x_err = np.square(cos_c*x_err)+np.square(x*np.sin(dec_rad)*self._data['dec_error'].values)
            x_err = np.square(x_err)
            x *= cos_c
        return pd.DataFrame({'pmra': x, 'pmdec': y, 'pmra_err': x_err, 'pmdec_err': y_err})

    def total_proper_motion(self):
        """
        Return the total proper motion and its errors.

        .. math::
            \mu = \sqrt{\left(\mu_\alpha^*\right)^2+\mu_\delta^2}

        with
        .. math::
            \mu_\alpha^* = \mu_\alpha*\cos (\delta)

        :return: pm, pm_err
        :rtype: numpy.ndarray, numpy.ndarray
        """
        pm_ra, pm_dec, pm_ra_error, pm_dec_error = self.proper_motion(True)
        pm = np.hypot(pm_ra, pm_dec)
        pm_err = np.square(pm_ra*pm_ra_error)+np.square(pm_dec*pm_dec_error)
        pm_err /= pm
        return pm, pm_err

    def distance(self, kind='bailer-jones'):
        """
        Compute and return the distance in kpc

        :param kind:
            The kind of transformation. Current options are 'bailer-jones' or 'bj'
            for the distances estimated by Bailer-Jones or 'simple' for :math:`distance=1/parallax`.
            All other inputs raise a ValueError.
        :type kind: str
        :return: distance and its error in kpc
        :rtype: numpy.ndarray, numpy.ndarray
        """
        kind = kind.lower()
        # if the required distance is the distance estimated by Bailer-Jones
        if kind == ('bailer-jones' or 'bj'):
            return self.data['rest']/1000
        elif kind == 'simple':
            distance = 1/self._data['parallax'].values
            distance_error = np.abs(1/self._data['parallax'].values**2)*self._data['parallax_error'].values
            return distance, distance_error
        else:
            raise ValueError(f'Unknown kind: {kind}')

    def set_parallax_limit(self, minimal, maximal, with_errors=False):
        """
        Set a limit to the parallax

        :param minimal: The minimal parallax
        :type minimal: float
        :param maximal: The maximal parallax
        :type maximal: float
        :param with_errors: True if the errors should be consider too, else False. Default is False.
        :type with_errors: bool
        :return:
        """
        parallax = self._data['parallax'].values
        if with_errors:
            parallax_error = self._data['parallax_error'].values
            m = (parallax+parallax_error >= minimal) & (parallax-parallax_error <= maximal)
        else:
            m = (parallax >= minimal) & (parallax <= maximal)
        self.mask.add_mask(m, f'Parallax constrain: {minimal} - {maximal}')

    def proper_motion_limit(self, minimal, maximal, with_errors=False):
        """
        Set a limit to the total proper motion.
        See :meth:`total_proper_motion` fot he computing details.


        :param minimal: The minimal total proper motion in mas/yr
        :type minimal: float
        :param maximal: The maximal total proper motion in mas/yr
        :type maximal: float
        :param with_errors: True if the errors should be consider too, else False. Default is False.
        :type with_errors: bool
        :return:
        """
        pm, pm_err = self.total_proper_motion()
        if with_errors:
            m = (pm + pm_err >= minimal) & (pm-pm_err <= maximal)
        else:
            m = (pm >= minimal) & (pm <= maximal)
        self.mask.add_mask(m, f'Proper motion limit from {minimal} to {maximal}')
