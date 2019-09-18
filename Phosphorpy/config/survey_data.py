import configparser
import pkg_resources
ADS_LINK = 'https://ui.adsabs.harvard.edu/?#abs/{}'


def create_dict(line):
    """
    Creates a dict from a line of the config file.

    :param line: A line from the config file
    :type line: str
    :return: A dict with one element
    :rtype: dict
    """
    line = line.split(' ')
    temp = [line_part for line_part in line if line_part != '']
    items = temp[1:]
    if len(items) == 1:
        items = items[0]
    return {temp[0]: items}


def read_survey_data():
    """
    Read the survey description from the local survey file

    :return: Dict with all local surveys
    """

    resource_package = 'Phosphorpy'
    resource_path = '/'.join(('local', 'survey.conf'))
    shortcut_path = pkg_resources.resource_filename(resource_package, resource_path)

    conf = configparser.ConfigParser()
    conf.read(shortcut_path)

    surveys = {}
    for s in conf.sections():
        c = conf[s]
        temp = {}
        for k in c:
            line = c[k].split(', ')
            if len(line) > 1:
                temp[k] = line
            else:
                temp[k] = c[k]
        surveys[c['name']] = temp
    return surveys


def add_survey(name, vizier_path, release, reference, magnitudes, coordinates=None, xmatch=None):

    resource_package = 'Phosphorpy'
    resource_path = '/'.join(('local', 'survey.conf'))
    shortcut_path = pkg_resources.resource_filename(resource_package, resource_path)

    conf = configparser.ConfigParser()
    conf.read(shortcut_path)

    sections = list(conf.sections())

    sid = 'survey{}'.format(int(sections[-1].split('survey')[-1])+1)
    conf.add_section(sid)

    conf.set(sid, 'name', value=name)
    conf.set(sid, 'release', release)
    conf.set(sid, 'reference', reference)
    conf.set(sid, 'magnitude', ', '.join(magnitudes))
    conf.set(sid, 'vizier', vizier_path)

    if coordinates is not None:
        conf.set(sid, 'coordinate', ', '.join(coordinates))

    if xmatch is not None:
        conf.set(sid, 'xmatch', xmatch)

    # conf.write()


SURVEY_DATA = read_survey_data()
