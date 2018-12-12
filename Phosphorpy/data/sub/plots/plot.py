class Plot:

    __plot_log = []

    def __init__(self, data_set):
        """

        :param data_set:
        :type data_set: Phosphorpy.data.data.DataSet
        """
        self._data_set = data_set

    def sed(self, *args, **kwargs):
        self._data_set.flux.plot.sed(*args, **kwargs)
        self.__plot_log.append((self._data_set.flux.plot.sed, args, kwargs))

    def color_color(self, *args, **kwargs):
        self._data_set.colors.plot.color_color(*args, **kwargs)
        self.__plot_log.append((self._data_set.colors.plot.color_color, args, kwargs))

    def color_hist(self, *args, **kwargs):
        self._data_set.colors.plot.color_hist(*args, **kwargs)
        self.__plot_log.append((self._data_set.colors.plot.color_hist, args, kwargs))

    def equatorial_coordinates(self):
        self._data_set.coordinates.plot.equatorial()
        self.__plot_log.append((self._data_set.coordinates.plot.equatorial, None, None))

    def galactic_coordinates(self):
        self._data_set.coordinates.plot.galactic()
        self.__plot_log.append((self._data_set.coordinates.plot.galactic, None, None))
