import pylab as pl
import seaborn


def replace_labels(axes, cols, labels):
    if type(labels) == list:
        axes.set_xlabel(labels[0])
        axes.set_ylabel(labels[1])
    elif type(labels) == dict:
        axes.set_xlabel(labels[cols[0]])
        axes.set_ylabel(labels[cols[1]])
    else:
        axes.set_xlabel(cols[0].replace('mag', ''))
        axes.set_ylabel(cols[1].replace('mag', ''))


class Color:
    """
    Class to handle different plot with the colors
    """

    _color = None
    _magnitude = None

    def __init__(self, color, magnitude=None):
        self._color = color
        self._magnitude = magnitude

    def __color_color_multi__(self, cols, labels):
        """
        Creates a color-color plot with multiple colors (more than 2) in a multiple grid

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :return:
        """
        pp = seaborn.PairGrid(self._color.data[cols])
        pp.map_diag(pl.hist)
        pp.map_offdiag(pl.scatter)
        pp.fig.subplots_adjust(wspace=0, hspace=0)

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

            pp.fig.tight_layout()

    def __color_color_single__(self, cols, labels):
        """]Creates a color-color plot with two magnitudes

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :return:
        """
        sp = pl.subplot()
        sp.scatter(self._color.data[cols[0]],
                   self._color.data[cols[1]],
                   marker='.', c='k')
        replace_labels(sp, cols, labels)

    def color_color(self, cols=None, path='', labels=None):
        """
        Plots a color color diagram. If their are more than two columns, it will
        plot the color-color diagrams in a grid

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param path:
            Path to the save place. Default is an empty string, which means that the figure will be shown, only.
        :type path: str
        :param labels:
            Replacement of the default labels in a list or a dict. Default is None, which means that default labels
            are used.
        :type labels: list, dict
        :return:
        """
        if cols is None:
            cols = self._color.data.columns

        pl.clf()
        if len(cols) > 2:
            self.__color_color_multi__(cols, labels)
            # todo: implement color color plot (grid if there are more than two colors)
            pass
        else:
            self.__color_color_single__(cols, labels)
        if path != '':
            pl.savefig(path)
        pl.show()

    def color_hist(self, cols=None, bins='auto', histtype='step', range=None, density=False, path='', labels=None):
        """
        Plots a histogram of the color(s).

        :param cols:
            The name of color columns. Default is None, which means that all colors are taken.
        :type cols: list
        :param bins: The number of bins in the histogram. Default is 'auto'.
        :type bins: int, str
        :param histtype: The type of the histogram. Default is 'step'.
        :type histtype: str
        :param range: The range of the x_axis. Default is None.
        :type range: list, tuple
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
        pl.clf()
        sp = pl.subplot()
        if type(cols) == list:
            for c in cols:
                sp.hist(self._color.data[c], bins=bins, histtype=histtype, range=range, density=density)
        else:
            sp.hist(self._color.data[cols], bins=bins, histtype=histtype, range=range, density=density)

        replace_labels(sp, cols, labels)
        if path != '':
            pl.savefig(path)
        pl.show()
