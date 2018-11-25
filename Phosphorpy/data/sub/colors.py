from .table import DataTable
from .plots.color import Color


class Colors(DataTable):

    def __init__(self, data):
        DataTable.__init__(self)
        self._data = data

        self._plot = Color(self)
