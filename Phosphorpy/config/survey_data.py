import configparser
import pkg_resources
import os
import urllib
import warnings


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
    
    # if the config file was not included in the package (for some reasons)
    if not os.path.exists(shortcut_path):
        warnings.warn('Survey configuration file is missing. Download it from GitHub.')
        url = 'https://raw.githubusercontent.com/patrickRauer/Phosphorpy/master/Phosphorpy/local/survey.conf'
        urllib.request.urlretrieve(url, shortcut_path)

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
    """
    Adds new survey data to the config file

    :param name: The name of the survey
    :type name: str
    :param vizier_path: The Vizier code to the survey
    :type vizier_path: str
    :param release: The current release of the data
    :type release: str
    :param reference: The link the publication
    :type reference: st
    :param magnitudes: The names of the magnitudes
    :type magnitudes: Union
    :param coordinates: The names of the coordinates
    :type coordinates: Union
    :param xmatch: XMatch specific column name sections
    :type xmatch: str
    :return:
    """
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


def get_survey_data(survey_name):
    """
    Returns the survey data of the survey with the specific name (lower or upper letters do not matter)

    :param survey_name: The name of the survey
    :type survey_name: str
    :return: The data of the survey, if one was found.
    :rtype: dict
    """
    survey_name = survey_name.lower()
    for k in SURVEY_DATA:
        if k.lower() == survey_name:
            return SURVEY_DATA[k]
    raise ValueError('Survey not found!')


SURVEY_DATA = read_survey_data()
