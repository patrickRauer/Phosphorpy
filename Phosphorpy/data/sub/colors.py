from .table import DataTable
from .plots.color import Color, create_color_name


class Colors(DataTable):

    def __init__(self, data, mask=None, survey_colors=None):
        DataTable.__init__(self, mask=mask)
        self._data = data

        self._survey_colors = survey_colors

        self._plot = Color(self)

    @property
    def survey_colors(self):
        return self._survey_colors

    def __get_mask_data__(self, col, minimum, maximum, previous):
        col = create_color_name(col)
        d = self._data[col].values
        mask = (d < maximum) & (d > minimum)
        self._mask.add_mask(mask, 'Color cut (minimum={}, maximum={}'.format(minimum, maximum), combine=previous)

    def set_limit(self, col, minimum=-99, maximum=99, previous=True):
        if type(col) == str:
            self.__get_mask_data__(col, minimum, maximum, previous)
        elif type(col) == list or type(col) == tuple:
            if len(col) == 2:
                for c in self._data.column:
                    if col[0] in c and col[1] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            elif len(col) == 3:
                for c in self._data.column:
                    if col[0] in c and col[1] in c and col[2] in c:
                        self.__get_mask_data__(c, minimum, maximum, previous)
                        break
            else:
                raise ValueError('The list must have 2 or 3 elements. Not more and not less.')
        else:
            raise ValueError('Sorry, I don\'t how I should handle col in the format {}. Please Use a string or '
                             'tuple/list to specify the right color')
        pass
