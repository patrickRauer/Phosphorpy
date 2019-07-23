import pandas as pd
from .table import DataTable
from .plots.color import Color, create_color_name
from Phosphorpy.data.sub.tables.color import Color as ColorTab


class Colors(DataTable):

    def __init__(self, data=None, mask=None, survey_colors=None):
        DataTable.__init__(self, mask=mask)
        if data is None:
            data = []
        self._data = data

        self._survey_colors = survey_colors

        self._plot = Color(self)

    @property
    def survey_colors(self):
        cols = {}
        for d in self.data:
            cols[d.survey_name] = d.columns
        return cols

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
        col = create_color_name(col)
        d = self._data[col].values
        mask = (d < maximum) & (d > minimum)
        self._mask.add_mask(mask, f'Color cut (minimum={minimum}, maximum={maximum})', combine=previous)

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
        for d in self.data:
            if survey is None or survey == d.survey_name:
                d.set_limit(col, minimum=minimum, maximum=maximum, previous=previous)

    def get_columns(self, cols):
        """
        Returns all colors with the given color names

        :param cols: The required colors
        :type cols: list
        :return: All colors with the given names in it
        :rtype: pd.DataFrame
        """
        out = None
        for d in self.data:
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
            if survey == d.survey_name:
                d.outlier_detection()
