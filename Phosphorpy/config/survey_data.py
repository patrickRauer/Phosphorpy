import configparser
import pkg_resources
ADS_LINK = 'https://ui.adsabs.harvard.edu/?#abs/{}'


def create_dict(line):
    line = line.split(' ')
    temp = []
    for line_part in line:
        if line_part != '':
            temp.append(line_part.strip(','))
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


SURVEY_DATA = read_survey_data()
