from SearchEngine.report.Report import Report

from SearchEngine.external.image import PanstarrsImage

from astropy.coordinates import SkyCoord
from astropy import units as u
from dominate import tags
from dominate.tags import div

import warnings
warnings.simplefilter('ignore')


def check_survey(data, survey):
    cols = data.colnames
    for c in cols:
        if survey in c and '_e_' not in c:
            if data[c] == data[c]:
                return True
    return False


@div
def mag_section(data):
    """
    Creates a HTML table with the magnitudes of a source
    :param data: The data of the source
    :return:
    """
    surveys = ['APASS', 'KIDS', 'SDSS12', '2MASS', 'VIKING', 'PS', 'UKIDSS']
    mag = ['u', 'B', 'V', 'g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
    cols = data.colnames
    with tags.table(_class='docutils', border='1'):
        with tags.tr():
            tags.td()
            for m in mag:
                tags.td(m)
        for s in surveys:
            if check_survey(data, s):
                with tags.tr():
                    tags.td(s)
                    for c in mag:
                        if '{}_{}mag'.format(s, c) in cols:
                            tags.td('{: 5.2f}'.format(data['{}_{}mag'.format(s, c)]))
                        elif '{}_{}mag'.format(s, c.lower()) in cols:
                            tags.td('{: 5.2f}'.format(data['{}_{}mag'.format(s, c.lower())]))
                        elif '{}_{}mag'.format(s, c.lower()+'_1') in cols:
                            tags.td('{: 5.2f}'.format(data['{}_{}mag'.format(s, c.lower()+'_1')]))

                        else:
                            tags.td()


@div
def image_section(coord):
    path = './_static/source_ps_image/source_{}.png'.format(coord.to_string('hmsdms').replace(' ', '_'))
    ps_img = PanstarrsImage()
    # ps_img.get_color_image(coord, path, smooth=5)
    return path


class Detail(Report):
    """
    Main class to create a detail view for a source
    """
    source_data = None

    def __init__(self, source_data):
        self.source_data = source_data

    @div
    def html(self, path=''):
        coords = SkyCoord(self.source_data['ra']*u.deg,
                          self.source_data['dec']*u.deg)
        s = coords.to_string('hmsdms')
        with tags.table():
            with tags.tr():

                tags.th('source {}'.format(s),
                        id=s.replace(' ', '_'))
            with tags.tr():
                with tags.td():
                    mag_section(self.source_data)
            with tags.tr():
                with tags.td():
                    img_path = image_section(coords)
                    with tags.a(href=img_path):
                        tags.img(src=img_path)
