from astropy import units as u
from astropy.coordinates import Angle

from Phosphorpy.data.sub.interactive_plotting.interactive_plotting import HVPlot

try:
    import holoviews as hv
except ImportError:
    hv = None


def _angle_plot(x, y, sp, marker, color):
    """

    :param x: X values
    :type x: numpy.ndarray
    :param y: y values
    :type y: numpy.ndarray
    :param sp: subplot environment
    :type sp: Tuple[Any, Any]
    :param marker: Marker type
    :type marker: str
    :param color: color
    :type color: str, None
    :return:
    """
    x = Angle(x*u.deg)
    x = x.wrap_at(180*u.deg)
    y = Angle(y*u.deg)

    sp.scatter(x.radian, y.radian, marker=marker, c=color)


class CoordinatePlot(HVPlot):
    """
    Basic class to provide coordinate plots
    """

    def __init__(self, coordinate):
        self._coordinate = coordinate

    def equatorial(self, path='', marker='.', color='k',
                   mask_colors=None, mollweide=False, **hv_kwargs):
        """
        Plot the position of the entries in equatorial coordinates

        :param path: Path the the storage place. An empty string is for showing only.
        :type path: str
        :param marker: The marker type. Default is '.'
        :type marker: str
        :param color: The color of the markers. Default is black
        :type color: str
        :param mask_colors: A color or a list of colors for the mask data
        :type mask_colors: str, list
        :param mollweide:
            True, if a mollweide projection should be used, else False.
            Default is True.
        :type mollweide: bool
        :return:
        """
        return self._scatter('ra', 'dec', 'R.A. [deg]', 'Dec [deg]',
                             path, marker, color, mask_colors=mask_colors,
                             mollweide=mollweide, **hv_kwargs)

    def galactic(self, path='', marker='.', color='k',
                 mask_colors=None, mollweide=False, **hv_kwargs):
        """
        Plot the position of the entries in galactic coordinates

        :param path: Path the the storage place. An empty string is for showing only.
        :type path: str
        :param marker: The marker type. Default is '.'
        :type marker: str
        :param color: The color of the markers. Default is black
        :type color: str
        :param mask_colors: A color or a list of colors for the mask data
        :type mask_colors: str, list
        :param mollweide:
            True, if a mollweide projection should be used, else False.
            Default is True.
        :type mollweide: bool
        :return:
        """
        return self._scatter('l', 'b', 'l [deg]', 'b [deg]',
                             path, marker, color, mask_colors=mask_colors,
                             mollweide=mollweide, **hv_kwargs)

    def _scatter(self, col1, col2, x_label, y_label, path, marker, color,
                 mask_colors=None, mollweide=True, **hv_kwargs):
        """
        Makes a scatter plot of two coordinates

        :param col1: The first column name of the coordinates
        :type col1: str
        :param col2: The second column name of the coordinates
        :type col2: str
        :param x_label: The x-label
        :type x_label: str
        :param y_label: The y-label
        :type y_label: str
        :param path: The path to the save place
        :type path: str
        :param marker: The marker type
        :type marker: str
        :param color: The color of the markers
        :type color: str
        :param mask_colors: A color or a list of colors for the mask data
        :type mask_colors: str, list
        :param mollweide:
            True, if a mollweide projection should be used, else False.
            Default is True.
        :type mollweide: bool
        :return:
        """

        if mollweide:
            raise NotImplementedError('Projection is not supported by holoviews. Use the non-interactive plotting '
                                      'for the mollweide-projection instead.')

        graph = hv.Scatter(self._coordinate.data, col1, col2)
        graph = self._hover(graph)
        if self._coordinate.mask.get_mask_count() > 0:
            for mask_id in range(self._coordinate.mask.get_mask_count()):
                g = hv.Scatter(self._coordinate.data[self._coordinate.mask.get_mask(mask_id)],
                               col1, col2)
                graph *= self._hover(g)

        graph = graph.opts(xlabel=x_label,
                           ylabel=y_label,
                           **hv_kwargs)
        return graph
