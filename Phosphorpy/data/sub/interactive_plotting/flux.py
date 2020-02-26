from Phosphorpy.data.sub.interactive_plotting.interactive_plotting import HVPlot

try:
    import holoviews as hv
except ImportError:
    hv = None


class FluxPlot(HVPlot):
    _data = None

    def __init__(self, data):
        """

        :param data:
        :type data: Phosphorpy.data.sub.flux.FluxTable
        """
        self._data = data

    def sed(self, index, fit=None, path='', x_log=False, y_log=False, legend=False, **hv_kwargs):
        """
        Plot the SED of the source with the given index.

        :param index: The index of the source
        :type index: int
        :param fit:
            A function to fit, a string with the name of a predefined function or an integer for a polynomial fit.
        :type fit: function, str, int
        :param path:
            Path to the place where the figure should be saved. Default is an empty string, which means that
            the figure will be shown only.
        :type path: str
        :param x_log: True if the x-axis should be in log scale, else False. Default is False.
        :type x_log: bool
        :param y_log: True if the y-axis should be in log-scale, else False. Default is False
        :type y_log: bool
        :param legend: True if the legend should be shown, else False. Default is False.
        :type legend: bool
        :return:
        """
        graph = None

        # iterate over all surveys
        for s in self._data.data:

            errs = hv.ErrorBars(
                (
                    self._data.survey_head.get_survey_wavelengths(s.survey_name),
                    s.get_flux(index),
                    s.get_error(index)
                )
            )

            g = hv.Scatter(
                (
                    self._data.survey_head.get_survey_wavelengths(s.survey_name),
                    s.get_flux(index)
                ),
                label=s.survey_name)
            g = errs * self._hover(g)

            if graph is None:
                graph = g
            else:
                graph *= g

        graph = graph.opts(
            xlabel='wavelength [$\\AA$]',
            ylabel='flux',
            logx=x_log,
            logy=y_log, **hv_kwargs
        )
        return graph
