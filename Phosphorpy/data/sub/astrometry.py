from Phosphorpy.external.vizier import Gaia
from .table import DataTable


class AstrometryTable(DataTable):

    def __init__(self, mask):
        DataTable.__init__(self, mask)

    @staticmethod
    def load_to_dataset(ds):
        g = Gaia()
        gaia = g.query(ds.coordinates, 'ra', 'dec', use_xmatch=True, blank=True)
        astronomy = AstrometryTable(ds.mask)
        astronomy._data = gaia[['ra', 'ra_error', 'dec', 'dec_error',
                                'parallax', 'parallax_error']]
        return astronomy
