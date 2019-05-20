class Plot:

    __plot_log = []

    def __init__(self, data_set):
        """
        This class is an interface to the plotting environments of the different sub tables.

        :param data_set:
        :type data_set: Phosphorpy.data.data.DataSet
        """
        self._data_set = data_set

    def magnitude_hist(self, *args, **kwargs):
        """
        Plots a magnitude histogram

        :param args:
        :param kwargs:
        :return:
        """
        self._data_set.magnitudes.plot.hist(*args, **kwargs)

    def sed(self, *args, **kwargs):
        """
        Plots the SED of the given target.

        :param args:
        :param kwargs:
        :return:
        """
        self._data_set.flux.plot.sed(*args, **kwargs)
        self.__plot_log.append((self._data_set.flux.plot.sed, args, kwargs))

    def color_color(self, *args, **kwargs):
        """
        Plots a color color graph of the targets for the given colors

        .. code-block:: python

            ds.plot.color_color(['BP - G', 'G - RP'])

        :param args:
        :param kwargs:
        :return:
        """
        self._data_set.colors.plot.color_color(*args, **kwargs)
        self.__plot_log.append((self._data_set.colors.plot.color_color, args, kwargs))

    def color_hist(self, *args, **kwargs):
        """
        Plot the distribution of colors for the given colors

        :param args:
        :param kwargs:
        :return:
        """
        self._data_set.colors.plot.color_hist(*args, **kwargs)
        self.__plot_log.append((self._data_set.colors.plot.color_hist, args, kwargs))

    def equatorial_coordinates(self, **kwargs):
        """
        Plot the coordinates of the targets in a equatorial coordinate system.
        The plot will be a mollweide projection of the coordinates

        .. code-block:: python

            ds.plot.equatorial_coordinates()

        :param kwargs:
            Coordinate plot arguments (path, marker and color).
            See :meth:`Phosphorpy.data.sub.plot.coordinates.equatorial` for details.
        :return:
        """
        self._data_set.coordinates.plot.equatorial(**kwargs)
        self.__plot_log.append((self._data_set.coordinates.plot.equatorial, None, None))

    def galactic_coordinates(self, **kwargs):
        """
        Plot the coordinates of the targets in a galactic coordinate system.
        The plot will be a mollweide projection of the coordinates

        :param kwargs:
            Coordinate plot arguments (path, marker and color).
            See :meth:`Phosphorpy.data.sub.plot.coordinates.galactic` for details.
        :return:
        """
        self._data_set.coordinates.plot.galactic(**kwargs)
        self.__plot_log.append((self._data_set.coordinates.plot.galactic, None, None))
