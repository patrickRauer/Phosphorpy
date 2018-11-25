class Plot:

    def __init__(self, data_set):
        """

        :param data_set:
        :type data_set: Phosphorpy.data.data.DataSet
        """
        self._data_set = data_set

    def sed(self, *args, **kwargs):
        self._data_set.flux.plot.sed(*args, **kwargs)

    def color_color(self, *args, **kwargs):
        self._data_set.colors.plot.color_color(*args, **kwargs)

    def color_hist(self, *args, **kwargs):
        self._data_set.colors.plot.color_hist(*args, **kwargs)
