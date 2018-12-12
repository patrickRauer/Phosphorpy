from astropy.coordinates import Angle
from astropy import units as u
import pylab as pl


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
    :type color: str
    :return:
    """
    x = Angle(x*u.deg)
    x = x.wrap_at(180*u.deg)
    y = Angle(y*u.deg)

    sp.scatter(x.radian, y.radian, marker=marker, c=color)


class CoordinatePlot:
    """
    Basic class to provide coordinate plots
    """

    def __init__(self, coordinate):
        self._coordinate = coordinate

    def equatorial(self, path='', marker='.', color='k'):
        """
        Plot the position of the entries in equatorial coordinates

        :param path:
        :param marker:
        :param color:
        :return:
        """
        self._scatter('ra', 'dec', 'R.A. [deg]', 'Dec [deg]',
                      path, marker, color)

    def galactic(self, path='', marker='.', color='k'):
        """
        Plot the position of the entries in galactic coordinates

        :param path:
        :param marker:
        :param color:
        :return:
        """
        self._scatter('l', 'b', 'l [deg]', 'b [deg]',
                      path, marker, color)

    def _scatter(self, col1, col2, x_label, y_label, path, marker, color):
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
        :return:
        """
        pl.clf()
        sp = pl.subplot(projection='mollweide')

        _angle_plot(self._coordinate[col1],
                    self._coordinate[col2],
                    sp, marker=marker,
                    color=color)
        if self._coordinate.mask.get_mask_count() > 0:
            # todo: implement color and/or marker iteration
            for mask_id in range(self._coordinate.mask.get_mask_count()):
                _angle_plot(self._coordinate[col1][self._coordinate.mask.get_mask(mask_id)],
                            self._coordinate[col2][self._coordinate.mask.get_mask(mask_id)],
                            sp, marker=marker,
                            color=None)
        sp.grid(True)
        sp.set_xticklabels(['14h', '16h', '18h', '20h', '22h',
                            '0h', '2h', '4h', '6h', '8h', '10h'])

        sp.set_xlabel(x_label)
        sp.set_ylabel(y_label)

        if path != '':
            pl.savefig(path)
        else:
            pl.show()
