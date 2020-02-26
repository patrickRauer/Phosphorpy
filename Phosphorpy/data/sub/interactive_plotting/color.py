from collections.abc import Iterable
import pylab as pl
import seaborn
import numpy as np

from Phosphorpy.data.sub.interactive_plotting.interactive_plotting import HVPlot

try:
    import holoviews as hv
    import numba
    from holoviews.operation.datashader import datashade
except ImportError:
    hv = None


def replace_labels(axes, cols, labels):
    """
    Sets the the labels of the graph in a proper format. For example removes 'mag' from the old labels

    :param axes: The axes object of the subplot

    :param cols: The names of the cols
    :type cols: list
    :param labels: The new labels
    :type labels: list, dict
    :return:
    """
    if type(labels) == list:
        axes.set_xlabel(labels[0])
        axes.set_ylabel(labels[1])
    elif type(labels) == dict:
        axes.set_xlabel(labels[cols[0]])
        axes.set_ylabel(labels[cols[1]])
    else:
        axes.set_xlabel(cols[0].replace('mag', ''))
        axes.set_ylabel(cols[1].replace('mag', ''))


def create_color_name(c):
    """
    Checks if the color name is a proper color name. If not it will convert the name to a proper style.
    :param c: The input color name
    :type c: str
    :return: The input color name in a proper style
    :rtype: str
    """
    if 'mag' not in c:
        c = c.split('-')
        c = '{}mag - {}mag'.format(c[0].replace(' ', ''),
                                   c[1].replace(' ', ''))
    return c


class Color(HVPlot):
    """
    Class to handle different plot with the colors
    """

    _color = None
    _magnitude = None

    def __init__(self, color, magnitude=None):
        self._color = color
        self._magnitude = magnitude

    def __color_color_multi__(self, cols, labels, legend=False):
        """
        Creates a color-color plot with multiple colors (more than 2) in a multiple grid

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :param legend: True if a legend with the selection labels should be shown, else False. Default is False.
        :type legend: bool
        :return:
        """
        # exclude data with nan values
        d = self._color.get_columns(cols)

        m = d[cols[0]] > -999
        for i in range(1, len(cols)):
            m = m & (d[cols[i]] > -999)

        hue = None
        if self._color.mask.get_mask_count() > 0:
            d['selection'] = '          '
            hue = 'selection'
            for i in range(1, self._color.mask.get_mask_count()):
                d.loc[self._color.mask.get_mask(i)]['selection'] = self._color.mask.get_description(i)
        d = d[m]
        pp = seaborn.PairGrid(d, hue=hue)
        pp.map_diag(pl.hist)
        pp.map_lower(pl.scatter, marker='.')

        pp.fig.subplots_adjust(wspace=0.02, hspace=0.02)

        # change the current labels, which are the default column names to a proper style
        # like removing 'mag' from the labels or replacing the labels with given labels
        for i in range(len(cols)):
            for j in range(len(cols)):
                axes = pp.axes[i][j]
                xlabel = axes.get_xlabel()
                ylabel = axes.get_ylabel()
                if type(labels) == list:
                    if xlabel != '':
                        axes.set_xlabel(labels[i])
                    if ylabel != '':
                        axes.set_ylabel(labels[j])
                elif type(labels) == dict:
                    if xlabel != '':
                        axes.set_xlabel(labels[xlabel])
                    if ylabel != '':
                        axes.set_ylabel(labels[ylabel])
                else:
                    if xlabel != '':
                        axes.set_xlabel(xlabel.replace('mag', ''))
                    if ylabel != '':
                        axes.set_ylabel(ylabel.replace('mag', ''))
        if hue is not None and legend:
            pp.add_legend()
            # pp.fig.tight_layout()

    def __color_color_survey__(self, survey, labels, legend=False):
        d = self._color.get_survey_data(survey)
        m = np.array(len(d)*[True])
        for c in d.columns:
            m = m & (d[c] > -999)

        cols = d.columns.values

        hue = None
        if self._color.mask.get_mask_count() > 0:
            d['selection'] = '          '
            hue = 'selection'
            for i in range(1, self._color.mask.get_mask_count()):
                d.loc[self._color.mask.get_mask(i)]['selection'] = self._color.mask.get_description(i)
        d = d[m]
        pp = seaborn.PairGrid(d, hue=hue)
        pp.map_diag(pl.hist)
        pp.map_lower(pl.scatter, marker='.')

        pp.fig.subplots_adjust(wspace=0.02, hspace=0.02)

        # change the current labels, which are the default column names to a proper style
        # like removing 'mag' from the labels or replacing the labels with given labels
        for i in range(len(cols)):
            for j in range(len(cols)):
                axes = pp.axes[i][j]
                xlabel = axes.get_xlabel()
                ylabel = axes.get_ylabel()
                if type(labels) == list:
                    if xlabel != '':
                        axes.set_xlabel(labels[i])
                    if ylabel != '':
                        axes.set_ylabel(labels[j])
                elif type(labels) == dict:
                    if xlabel != '':
                        axes.set_xlabel(labels[xlabel])
                    if ylabel != '':
                        axes.set_ylabel(labels[ylabel])
                else:
                    if xlabel != '':
                        axes.set_xlabel(xlabel.replace('mag', ''))
                    if ylabel != '':
                        axes.set_ylabel(ylabel.replace('mag', ''))
        if hue is not None and legend:
            pp.add_legend()
            # pp.fig.tight_layout()
        pass

    def __color_color_single__(self, cols, labels, legend=False, use_datashade=True, **hv_kwargs):
        """
        Creates a color-color plot with two magnitudes

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :param legend: True if a legend with the selection labels should be shown, else False. Default is False.
        :type legend: bool
        :return:
        """
        # todo: make sure that one survey is used
        color1 = self._color.get_column(cols[0])
        color2 = self._color.get_column(cols[1])

        color1, color2 = color1.align(color2)

        graph = hv.Points(
            (color1,
             color2)
        )

        if len(color1) > 1000 and use_datashade:
            graph = datashade(graph)
        graph = self._hover(graph)

        # iterate over all masks
        for i in range(self._color.mask.get_mask_count()):
            m = self._color.mask.get_mask(i)
            g = hv.Points(
                (
                    color1[m],
                    color2[m]
                ),
                label=self._color.mask.get_description(i)
            )
            if len(color1) > 1000 and use_datashade:
                g = datashade(g)
            g = self._hover(g)
            graph *= g

        graph = graph.opts(
            xlabel=cols[0].replace('mag', ''),
            ylabel=cols[1].replace('mag', ''),
            **hv_kwargs
        )

        return graph

    def color_color(self, survey=None, cols=None, path='', labels=None, legend=False, use_datashade=True, **hv_kwargs):
        """
        Plots a color color diagram. If their are more than two columns, it will
        plot the color-color diagrams in a grid

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list, dict
        :param path:
            Path to the save place. Default is an empty string, which means that the figure will be shown, only.
        :type path: str
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :param survey: A name of a loaded survey to make a color-color plot with the colors of this survey.
        :type survey: str
        :param legend: True if a legend with the selection labels should be shown, else False. Default is False.
        :type legend: bool
        :return:
        """
        # todo: if a survey and cols are set, check if the cols are in the survey colors and then use this colors only
        if cols is None:
            if survey is None:
                cols = self._color.data.columns.values
            else:
                self.__color_color_survey__(survey, labels, legend=legend)
                return
                # cols = self._color.survey_colors[survey].values
        else:
            use_cols = []
            if type(cols) is list:
                for c in cols:
                    use_cols.append(create_color_name(c))
            elif type(cols) == dict:
                # todo: implement a col selection with a dict (maybe a list of dicts?)
                raise NotImplementedError()
            else:
                raise ValueError(f'{type(cols)} is not supported format for color names.')
            cols = use_cols

        if len(cols) > 2:
            raise NotImplementedError('Multi color-color plots are not supported in the interactive plotting '
                                      'environment.\n Choose the standard plotting environment to create these plots.')
        else:
            return self.__color_color_single__(cols, labels, legend=legend, use_datashade=use_datashade, **hv_kwargs)

    def color_hist(self, cols=None, bins='auto', histtype='step', xlimit=None, density=False, path='',
                   labels=None, **hv_kwargs):
        """
        Plots a histogram of the color(s).

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param bins: The number of bins in the histogram. Default is 'auto'.
        :type bins: int, str
        :param histtype: The type of the histogram. Default is 'step'.
        :type histtype: str
        :param xlimit: The range of the x_axis. Default is None.
        :type xlimit: list, tuple
        :param density: True if the histogram should be a density histogram, else False. Default is False.
        :type density: bool
        :param path:
            Path to the save place. Default is an empty string, which means that the figure will be shown, only.
        :type path: str
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :return:
        """

        if type(cols) == str:
            cols = [cols]
        # exclude data with nan values
        # m = self._color.get_columns(cols[0]) > -999
        # for i in range(1, len(cols)):
        #     m = m & (self._color.get_columns(cols[i]) > -999)

        if density:
            ylabel = 'density'
        else:
            ylabel = '#'

        graph = None

        if isinstance(cols, Iterable):
            for c in cols:
                color = self._color.get_columns(c)
                m = color == color
                # hist, edge = np.histogram(color[m.values], bins=bins,
                #                           range=xlimit, density=density)
                g = hv.Distribution(color[m.values], label=c.replace('mag', ''))
                # g = hv.Histogram((hist, edge), label=c.replace('mag', ''))

                g = self._hover(g)

                if graph is None:
                    graph = g
                else:
                    graph *= g
        else:
            color = self._color.get_columns(cols)
            m = color == color
            hist, edge = np.histogram(color[m.values], bins=bins,
                                      range=xlimit, density=density)
            graph = hv.Histogram((hist, edge)).opts(tools=['hover'])
            graph = self._hover(graph)

        graph = graph.opts(
            xlabel='color',
            ylabel=ylabel,
            **hv_kwargs
        )

        return graph
