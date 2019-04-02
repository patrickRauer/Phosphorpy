from dominate import tags, document
from dominate.tags import div

from astropy.coordinates import Angle
from astropy import units as u
import pylab as pl
import os


class Report:
    root_path = './'
    static_path = '_static/'
    alabaster_link = 'alabaster.css'
    parts = []

    def add_section(self, sec):
        self.parts.append(sec)

    def html(self, path=''):
        doc = document(title='source sheet')

        with doc.head:
            tags.script(src='https://code.jquery.com/jquery-3.3.1.js',
                        type='text/javascript')
            tags.script(src='https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js',
                        type='text/javascript')
            tags.script(src="https://ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/jquery-ui.min.js",
                        type='text/javascript')
            tags.script(src='coordinates.js',
                        type='text/javascript')
            tags.script(src="https://cdn.datatables.net/colreorder/1.1.2/js/dataTables.colReorder.min.js",
                        type='text/javascript')
            tags.script(src="https://cdn.datatables.net/plug-ins/725b2a2115b/api/fnFilterClear.js",
                        type='text/javascript')
            tags.script(src="https://cdn.jsdelivr.net/jquery.ui-contextmenu/1.7.0/jquery.ui-contextmenu.min.js",
                        type='text/javascript')
            tags.link(rel='stylesheet', href=self.alabaster_link)
            tags.link(rel='stylesheet', href='default.css')
            tags.script(src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/' +
                            'MathJax.js?config=TeX-AMS-MML_HTMLorMML',
                        type='text/javascript')
            tags.link(rel='stylesheet', href="https://cdn.datatables.net/v/dt/dt-1.10.18/datatables.min.css",
                      type='text/css')
        with doc:
            tags.nav()
            with div(_class='document'):
                with div(_class='documentwrapper'):
                    with div(_class='bodywrapper'):
                        with div(_class='body', role='main'):
                            for p in self.parts:
                                p.html()
        render = doc.render()

        if path != '':
            with open(path, 'w') as f:
                f.write(render)
        return render


class SourceList(Report):
    """
    Creates a HTML-table with the sources
    """
    __ra_string__ = '{:02.0f}:{:02.0f}:{:04.1f}'
    __dec_string__ = '{:+03.0f}:{:02.0f}:{:04.1f}'
    cds_link = 'http://cdsportal.u-strasbg.fr/?target={} {}'
    ned_link = 'https://ned.ipac.caltech.edu/?q=nearposn&lon={}&lat={}&sr=2&incsrcs=0&coordsys=Equatorial&equinox=J2000'
    sdss_link = 'http://skyserver.sdss.org/dr12/en/tools/explore/summary.aspx?ra={}&dec={}'
    coordinates = None

    def __init__(self, coordinates):
        self.coordinates = coordinates

    @div(_class='dropdown')
    def __id_div__(self, id_number, coordinate):
        tags.label(str(id_number), id='sc_{}'.format(id_number))
        with div(_class='dropdown-content'):
            # link to the details of the source
            tags.a('details',
                   href='#{}'.format(coordinate.to_string('hmsdms').replace(' ', '_')))
            # link to CDS
            tags.a('CDS',
                   href=self.cds_link.format(coordinate.ra.degree,
                                             coordinate.dec.degree),
                   target='_blank')
            # link to SDSS
            tags.a('SDSS',
                   href=self.sdss_link.format(coordinate.ra.degree,
                                              coordinate.dec.degree),
                   target='_blank')
            # link to NED
            tags.a('NED',
                   href=self.ned_link.format(*(coordinate.to_string('hmsdms').split(' '))),
                   target='_blank')
            # link to CSS
            # TODO: check if this is possible
            tags.a('CSS',
                   href='',
                   target='_blank')

    @div(id='source_list')
    def html(self, path=''):

        with tags.table(_class='display', _id='example'):
            # table header

            # table body
            with open('coordinates.js', 'w') as f:
                f.write('var data = [\n')
                for i, c in enumerate(self.coordinates):
                    ra = self.__ra_string__.format(*c.ra.hms)
                    dec = self.__dec_string__.format(c.dec.dms[0],
                                                     abs(c.dec.dms[1]),
                                                     abs(c.dec.dms[2]))
                    f.write(f'\t\t["{i}", "{ra}", "{dec}"],\n')
                f.write('\t]')

            tags.script("$(document).ready(function () {\n" +
                        "\t$('#example').DataTable( {\n" +
                        "\t\tdata: data,\n" +
                        "\t\tcolumns: [{ title: 'id' },\n" +
                        "\t\t\t{ title: 'ra' },\n" +
                        "\t\t\t{ title: 'dec' }]\n" +
                        "\t} );\n" +
                        "\tconsole.log('create table');\n" +
                        "\t});\n" +
                        "$(document).contextmenu({\n" +
                        "\tdelegate: '.dataTable td',\n" +
                        "\tmenu: [\n" +
                        "\t\t{title: 'Filter', uiIcon: 'ui-icon-volume-off ui-icon-filter'},\n" +
                        "\t\t{title: 'Remove filter', uiIcon: 'ui-icon-volume-off ui-icon-filter'}\n" +
                        "]});\n")


class ReportPlot(Report):
    """
    Basic class to create static plots for the Report-HTML file
    """
    image_name = 'dummy.png'
    mask = None

    def __init__(self, mask=None):
        self.mask = mask

    def _scatter_(self, sp, x, y):
        """
        Does a scatter plot. If data are masked, they will be plotted in a different color.

        :param sp: subplot environment
        :type sp: Tuple[Any, Any]
        :param x: The x-values
        :type x: list
        :param y: The y-values
        :type y: list
        :return:
        """

        sp.scatter(x, y,
                   marker='.', c='k')

        if self.mask is not None:
            sp.scatter(x[self.mask],
                       y[self.mask],
                       marker='.', c='r')

    def _check_dir_(self):
        path = ''.join([self.root_path, self.static_path])
        if not os.path.exists(path):
            os.makedirs(path)


class SourcePositionPlot(ReportPlot):
    """
    Plot environment for position plots (equatorial and galactic
    """
    image_name = 'position.png'
    projection = None
    galactic = False

    def __init__(self, coordinates, projection=None, galactic=False, mask=None):
        ReportPlot.__init__(self, mask=mask)
        self.coordinates = coordinates
        self.projection = projection
        self.galactic = galactic
        if galactic:
            self.image_name = 'position_galactic.png'

    def __create_image__(self):
        """
        Creates an image of the positions of the sources and save it in the static image directory
        :return:
        """
        pl.clf()
        # choose the coordinate type
        if self.galactic:
            x = self.coordinates.galactic.l
            y = self.coordinates.galactic.b
        else:
            x = self.coordinates.ra
            y = self.coordinates.dec

        if self.projection is not None:
            sp = pl.subplot(projection=self.projection)
            x = Angle(x)
            x = x.wrap_at(180 * u.degree).radian
            y = Angle(y).radian
        else:
            sp = pl.subplot()

        # do a scatter plot
        self._scatter_(sp, x, y)

        pl.grid(True)

        # choose the right labeling for the different coordinate types
        if self.galactic:
            sp.set_xlabel('l [deg]')
            sp.set_ylabel('b [deg]')
        else:
            sp.set_xlabel('RA [deg]')
            sp.set_ylabel('Dec [deg]')

        path = ''.join([self.root_path, self.static_path])

        # check if the static dir exists
        self._check_dir_()
        pl.savefig(path+self.image_name)

    @div(id='position_plot')
    def html(self, path=''):
        self.__create_image__()
        path = ''.join([self.root_path, self.static_path, self.image_name])
        with tags.a(href=path):
            tags.img(src=path)


def split_magnitude_name(c):
    c = c.split('_')[1]
    c = c.split('mag')[0]
    return c


class MagnitudePlot(ReportPlot):

    data = None
    col1 = None
    col2 = None
    col3 = None
    col4 = None

    def __init__(self, data, col1, col2, col3, col4, mask=None):
        ReportPlot.__init__(self, mask)
        self.data = data
        self.col1 = col1
        self.col2 = col2
        self.col3 = col3
        self.col4 = col4

    def __create_image__(self):
        pl.clf()
        sp = pl.subplot()
        x = self.data[self.col1]-self.data[self.col2]
        y = self.data[self.col3]-self.data[self.col4]

        self._scatter_(sp, x, y)

        sp.set_xlabel('{} - {}'.format(split_magnitude_name(self.col1),
                                       split_magnitude_name(self.col2)))
        sp.set_ylabel('{} - {}'.format(split_magnitude_name(self.col3),
                                       split_magnitude_name(self.col4)))

        self._check_dir_()

        pl.savefig(self.root_path+self.static_path+'{}-{}_{}-{}.png'.format(self.col1,
                                                                            self.col2,
                                                                            self.col3,
                                                                            self.col4))

    @div
    def html(self, path=''):
        self.__create_image__()
        tags.img(src=self.root_path+self.static_path+'{}-{}_{}-{}.png'.format(self.col1,
                                                                              self.col2,
                                                                              self.col3,
                                                                              self.col4))


class DataSetReport(Report):

    _ds = None

    def __init__(self, dataset):
        """
        Creates a html report of the given dataset
        :param dataset: The dataset for the report
        :type dataset: Phosphorpy.DataSet
        """
        self._ds = dataset

    def html(self, path=''):
        self.parts.append(SourceList(self._ds.coordinates.as_sky_coord()))

        Report.html(self, path)
