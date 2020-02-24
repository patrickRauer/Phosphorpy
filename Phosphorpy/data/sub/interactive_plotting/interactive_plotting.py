try:
    import holoviews as hv
except ImportError:
    hv = None


class HVPlot:
    BOKEH = 'bokeh'
    PLOTLY = 'plotly'
    MATPLOTLIB = 'matplotlib'

    def _hover(self, g):
        """
        Add the parameter tools with the value 'hover' to the graph options, if the backend is 'bokeh'.

        :param g: The graph
        :return:
        """
        if hv.Store.current_backend == self.BOKEH:
            g = g.opts(tools=['hover'])
        return g

    @staticmethod
    def holoviews():
        if hv is not None:
            return True
        else:
            return False
