from .survey_data import SURVEY_DATA

# Set the names of the surveys
# TODO: use a file instead of hard coding
PANSTARRS = 'PAN-STARRS'
SDSS = 'SDSS'
WISE = 'WISE'
TWO_MASS = '2MASS'
APASS = 'APASS'
KIDS = 'KiDS'
VIKING = 'VIKING'
GALEX = 'GALEX'
GAIA = 'GAIA'
UKIDSS = 'UKIDSS'
BAILER_JONES = 'Bailer-Jones'
SkyMapper = 'SkyMapper'

# Dict where the keys are the survey names and
# the items are a list of the corresponding filter names
MAG_COLUMNS = {}
for k in SURVEY_DATA:
    MAG_COLUMNS[k] = SURVEY_DATA[k]['magnitude']


def get_available_surveys():
    """
    Returns a list with the available surveys

    :return: The name of the surveys
    :rtype: list
    """
    surveys = list(MAG_COLUMNS.keys())
    return surveys


def get_survey_magnitude_names(survey_name):
    """
    Returns a list with the filter names of the survey

    :param survey_name: The name of the survey
    :type survey_name: str
    :return: List with the filter names
    :rtype: list
    """
    if survey_name in MAG_COLUMNS.keys():
        return MAG_COLUMNS[survey_name]
    else:
        raise AttributeError('No survey with name \'{}\' is in the list!'.format(survey_name))
