import numpy as np
import pandas as pd
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, Distance
from astropy import units as u

from Phosphorpy.data.sub.plots.astrometry import AstrometryPlot
from Phosphorpy.external.vizier import get_survey
from .table import DataTable

from Phosphorpy.data.sub.interactive_plotting.astrometry import AstrometryPlot as AstrometryPlotHV


Gaia = get_survey('gaia')
BailerJones = get_survey('bailer-jones')


def _only_nearest(data):
        row_ids, row_id_count = np.unique(data['row_id'], return_counts=True)
        # find or multiple detections the closest one
        nearest = []
        for rid in row_ids[row_id_count >= 1]:
            g = data[data['row_id'] == rid]
            nearest.append(g[g['angDist'] == np.min(g['angDist'])])
        nearest = vstack(nearest)
        # data = vstack([nearest, data])
        return nearest


def _download_gaia_data(coordinates):
    """
    Downloads the Gaia data for the given coordinates

    :param coordinates: The coordinates of the targets for which Gaia data are required
    :type coordinates: Phosphorpy.data.sub.coordinates.CoordinateTable
    :return: The Gaia data
    :rtype:
    """
    g = Gaia

    gaia = g.query(coordinates.to_astropy_table(), 'ra', 'dec', use_xmatch=True, blank=True)

    # find or multiple detections the closest one
    # gaia = _only_nearest(gaia)

    if type(gaia) == Table:
        gaia = gaia.to_pandas()
    gaia = gaia.drop_duplicates('row_id')
    gaia_coords = gaia[['row_id', 'ra', 'dec']]

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
    bj = BailerJones
    if type(gaia_coords) == pd.DataFrame:
        gaia_coords = Table.from_pandas(gaia_coords)
    elif type(gaia_coords) != Table:
        gaia_coords = Table.from_pandas(gaia_coords.data)
    bj = bj.query(gaia_coords, 'ra', 'dec', use_xmatch=True, blank=True)

    bj = bj[['row_id', 'rest', 'b_rest', 'B_rest', 'rlen', 'ResFlag', 'ModFlag']]
    if type(bj) == Table:
        bj = bj.to_pandas()
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

        if AstrometryPlotHV.holoviews():
            self._hv_plot = AstrometryPlotHV(self)

    @staticmethod
    def load_to_dataset(ds):
        astrometry = AstrometryTable.load_astrometry(ds.coordinates)

        ds.set_astrometry(astrometry)

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

        # df = pd.DataFrame()
        # df['row_id'] = np.arange(coordinates.data.index.values.min(),
        #                          coordinates.data.index.values.max())
        # df = df.set_index('row_id')
        # gaia = df.join(gaia, how='outer')

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

        temp = self._data[['pmra', 'pmdec', 'pmra_error', 'pmdec_error']].copy()
        x = temp['pmra'].values
        y = temp['pmdec'].values
        x_err = temp['pmra_error'].values
        y_err = temp['pmdec_error'].values
        if cos_correction:
            dec_rad = np.deg2rad(self._data['dec'].values)
            cos_c = np.cos(dec_rad)
            x_err = np.square(cos_c*x_err)+np.square(x*np.sin(dec_rad)*self._data['dec_error'].values/1000)
            x_err = np.sqrt(x_err)
            k = x_err > 10
            if np.sum(k) > 0:
                x_err[k] = 10
            x *= cos_c
        return pd.DataFrame({'pmra': x, 'pmdec': y, 'pmra_err': x_err, 'pmdec_err': y_err}, index=self._data.index)

    def total_proper_motion(self):
        """
        Return the total proper motion and its errors.

        .. math::
            \mu = \sqrt{\left(\mu_{\\alpha^*}\\right)^2+\mu_\delta^2}

        with

        .. math::
            \mu_{\\alpha^*} = \mu_{\\alpha}*\cos (\delta)

        :return: The total proper motion and its errors in a DataFrame
        :rtype: pandas.DataFrame
        """
        temp = self.proper_motion(True)
        pm_ra, pm_dec = temp['pmra'].values, temp['pmdec'].values
        pm_ra_error, pm_dec_error = temp['pmra_err'].values, temp['pmdec_err'].values
        pm = np.hypot(pm_ra, pm_dec)
        pm_err = np.square(pm_ra*pm_ra_error)+np.square(pm_dec*pm_dec_error)
        pm_err /= pm
        return pd.DataFrame({'pm': pm, 'err': pm_err}, index=self._data.index)

    def distance(self, kind='bailer-jones'):
        """
        Compute and return the distance in kpc

        :param kind:
            The kind of transformation. Current options are 'bailer-jones' or 'bj'
            for the distances estimated by Bailer-Jones or 'simple' for :math:`distance=1/parallax`.
            All other inputs raise a ValueError.
        :type kind: str
        :return: distance and its error in kpc
        :rtype: pd.DataFrame
        """
        kind = kind.lower()
        # if the required distance is the distance estimated by Bailer-Jones
        if kind in ('bailer-jones', 'bj'):
            lower = self.data['rest'].values-self.data['b_rest'].values
            upper = self.data['B_rest'].values-self.data['rest'].values
            return pd.DataFrame({'distance': self.data['rest'].values/1000,
                                 'error': (lower+upper)/2},
                                index=self._data.index)
        elif kind in ('simple', 'parallax'):
            distance = 1/self._data['parallax'].values
            distance_error = np.abs(1/self._data['parallax'].values**2)*self._data['parallax_error'].values
            return pd.DataFrame({'distance': distance, 'error': distance_error}, index=self._data.index)
        else:
            raise ValueError(f'Unknown kind: {kind}')

    def set_distance_limit(self, minimal, maximal, kind='bailer-jones'):
        """
        Set a limit to the distances (in kpc)

        :param minimal: The minimal distance
        :type minimal: float
        :param maximal: The maximal distance
        :type maximal: float
        :param kind:
            The kind of transformation. Current options are 'bailer-jones' or 'bj'
            for the distances estimated by Bailer-Jones or 'simple' for :math:`distance=1/parallax`.
            All other inputs raise a ValueError.
        :type kind: str
        :return:
        """
        if minimal >= maximal:
            raise ValueError(f'Minimal distance must be larger than maximal distance: {minimal} >= {maximal}.')
        d = self.distance(kind)
        d = d['distance']
        m = (d >= minimal) & (d <= maximal)
        self.mask.add_mask(m, f'Distance limit: {minimal} - {maximal}')

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
        if minimal >= maximal:
            raise ValueError(f'Minimal value must be smaller than maximal: {minimal} > {maximal}')
        parallax = self._data['parallax']
        if with_errors:
            parallax_error = self._data['parallax_error']
            m = (parallax+parallax_error >= minimal) & (parallax-parallax_error <= maximal)
        else:
            m = (parallax >= minimal) & (parallax <= maximal)
        self.mask.add_mask(m, f'Parallax constrain: {minimal} - {maximal}')

    def set_total_proper_motion_limit(self, minimal, maximal, with_errors=False):
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
        if minimal >= maximal:
            raise ValueError(f'Minimal value must be smaller than maximal: {minimal} > {maximal}')

        pm = self.total_proper_motion()
        if with_errors:
            m = (pm['pm'] + pm['err'] >= minimal) & (pm['pm']-pm['err'] <= maximal)
        else:
            m = (pm['pm'] >= minimal) & (pm['pm'] <= maximal)
        self.mask.add_mask(m, f'Proper motion limit from {minimal} to {maximal}')

    def set_proper_motion_limit(self, direction, minimal, maximal, with_errors=False, cos_correction=False):
        """
        Set a limit to the proper motion

        :param direction: The direction of the proper motion (RA or Dec)
        :type direction: str
        :param minimal: The minimal proper motion
        :type minimal: float
        :param maximal: The maximal proper motion
        :type maximal: float
        :param with_errors: True if the errors should be considered in the mask, else False. Default is False.
        :type with_errors: bool
        :param cos_correction: True if the proper motion in RA direction should be declination correct, else False.
        :return:
        """
        direction = direction.lower()
        if direction != 'ra' and direction != 'dec':
            raise ValueError('Direction must be RA or Dec!')

        if minimal >= maximal:
            raise ValueError(f'Minimal value must be smaller than maximal: {minimal} > {maximal}')

        pm = self.proper_motion(cos_correction)
        pm_d = pm[f'pm{direction}']
        if with_errors:
            pm_d_err = pm[f'pm{direction}_err']
            m = (pm_d + pm_d_err >= minimal) & (pm_d-pm_d_err <= maximal)
        else:
            m = (pm_d >= minimal) & (pm_d <= maximal)
        self.mask.add_mask(m, f'Proper motion limit in {direction} direction from {minimal} to {maximal}')

    def to_sky_coord(self, bailer_jones=False):
        """
        Creates from the Gaia data a astropy.coordinates.SkyCoord object for every source.
        The SkyCoord contains beside the coordinates, the proper motions (pm ra is corrected)
        and a distance estimation.

        :param bailer_jones:
            True, if the distance estimation from Bailer-Jones et al. should be used.
            Else False for the Gaia DR2's parallax. In this case negative parallaxes are replaced with NaN.
        :type bailer_jones: bool
        :return: SkyCoord objects of the detections
        :rtype: SkyCoord
        """
        if bailer_jones:
            distance = Distance(self._data['rest']*u.pc)
        else:
            m = self._data['parallax'] == self._data['parallax']
            distance = self._data['parallax'].values
            distance[m] = np.nan
            distance = Distance(parallax=distance*u.mas)
        s = SkyCoord(self._data['ra'].values*u.deg,
                     self._data['dec'].values*u.deg,
                     pm_ra_cosdec=self._data['pmra'].values*np.cos(np.deg2rad(self._data['dec'].values))*u.mas/u.yr,
                     pm_dec=self._data['pmdec'].values*u.mas/u.yr,
                     distance=distance)
        return s
