import pandas as pd

from Phosphorpy.data.sub.tables.color import Color as ColorTab
from .plots.color import Color
from .table import DataTable

try:
    from Phosphorpy.data.sub.interactive_plotting.color import Color as ColorHV
except ImportError:
    ColorHV = None


class Colors(DataTable):

    def __init__(self, data=None, mask=None, survey_colors=None):
        DataTable.__init__(self, mask=mask)
        if data is None:
            data = []
        self._data = data

        self._survey_colors = survey_colors

        self._plot = Color(self)

        if ColorHV.holoviews():
            self._hv_plot = ColorHV(self)

    def __str__(self):
        s = 'Colors:\n'
        for d in self.data:
            s += str(d) + '\n'
        return s[:-2]

    @property
    def survey_colors(self):
        cols = {}
        for d in self.data:
            cols[d.survey_name] = d.columns
        return cols

    def get_survey_data(self, name):
        for d in self.data:
            if d.has_name(name):
                return d
        raise ValueError(f'No survey with {name} found.')

    def add_colors(self, data, survey_name):
        """
        Add a new ColorTab to the colors

        :param data: The new color table
        :type data: pandas.DataFrame
        :param survey_name: The name of the survey of the colors
        :type survey_name: str
        :return:
        """
        self.data.append(ColorTab(data, survey_name))

    def __get_mask_data__(self, col, minimum, maximum, previous):
        for d in self.data:
            d.set_limit(col, minimum, maximum)

    def set_limit(self, col, minimum=-99, maximum=99, previous=True, survey=None):
        """
        Sets a constrain to the colors and create a new mask of it

        :param col: The columns of the constrain.
        :type col: str, list, tuple
        :param survey:
            Name of the survey or None. If None all surveys with such a color name are used for the limiting
        :type survey: str, None
        :param minimum: The minimal value
        :type minimum: float
        :param maximum: The maximal value
        :type maximum: float
        :param previous: True if the last mask must be True too, else False to create a complete new mask.
        :type previous: bool
        :return:
        """

        if minimum >= maximum:
            raise ValueError(f'Minimum must be bigger than maximum: {minimum} >= {maximum}')
        for d in self.data:
            if survey is None or d.has_name(survey):
                d.set_limit(col, minimum=minimum, maximum=maximum, previous=previous)

    def get_columns(self, cols):
        """
        Returns all colors with the given color names

        :param cols: The required colors
        :type cols: list, str
        :return: All colors with the given names in it
        :rtype: pd.DataFrame
        """
        out = None
        for d in self.data:
            try:
                o = d.get_columns(cols)
                if len(o.columns) > 0:
                    if out is None:
                        out = o
                    else:
                        if type(cols) == pd.core.indexes.base.Index:
                            if len(o.columns) > len(out.columns):
                                out = o
                        else:
                            out = out.join(o)
            except ValueError:
                pass
        return out

    def get_column(self, col):
        """
        Returns the color with the given name

        :param col: The required color name
        :type col: str
        :return:
        """
        for d in self.data:
            if col in d.columns.values:
                return d[col]

    def outlier_detection(self, survey):
        """
        Detects outliers in the color data and mask them

        :param survey: The name of the survey where the outlier detection should be performed
        :type survey: str
        :return:
        """
        survey = survey.lower()
        for d in self.data:
            if d.has_name(survey):
                d.outlier_detection()
