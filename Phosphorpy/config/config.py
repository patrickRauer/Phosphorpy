from collections import OrderedDict
from astropy import units as u
from threading import Thread
import configparser
from SearchEngine.external.vizier import Vizier
from SearchEngine.report.Report import Report
import os


class MultiDict(OrderedDict):
    _unique = 0   # class variable

    def __setitem__(self, key, val):
        if isinstance(val, dict):
            self._unique += 1
            key += str(self._unique)
        OrderedDict.__setitem__(self, key, val)


class SearchConfig(Thread):
    data = None
    commands = []

    def __init__(self, path=None):
        """
        Class to represent the configuration settings of the search

        :param path: Path to the config-file. Default is None, which means default settings are used.
        :type path: str
        """
        Thread.__init__(self)
        self.commands = configparser.ConfigParser(dict_type=MultiDict, strict=False)
        self.commands.read(path)
        self.commands.sections()

    def __do_query__(self, query_data):
        mode = query_data['mode'].lower()
        if 'circle' in mode or 'box' in mode:
            viz = Vizier(name=query_data['survey'].upper())

            # convert central coordinates
            ra, dec = query_data['center'].split(', ')
            ra = float(ra)
            dec = float(dec)

            # circle mode
            if 'circle' in mode:
                # set the vizier radius with the chosen radius
                viz.radius = float(query_data['radius'])*__get_units__(query_data)
                print(viz.radius)
            # box mode
            elif 'box' in mode:
                r = query_data['radius'].split(', ')
                unit = __get_units__(query_data)
                viz.radius = [r[0]*unit, r[1]*unit]
            else:
                raise AttributeError('{} is not allowed as query mode type.\n' +
                                     'Allowed mode types are \'circle\' or \'box\'.')

            # query around coordinates
            self.data = viz.query_coordinate([ra, dec])
        print('query mode', mode)

    def __create_report__(self, report_data):
        r = Report()
        pass

    def __do_command__(self, command,
                       command_data):
        com = command.lower()
        if 'query' in com:
            print(command_data)
            self.__do_query__(command_data)
        elif 'match' in com:
            pass
        elif 'image' in com:
            pass
        elif 'load' in com:
            pass
        elif 'write' in com or 'save' in com:
            pass
        elif 'report' in com:
            self.__create_report__(command_data)

    def run(self):
        for command in self.commands.sections():
            self.__do_command__(command, dict(self.commands.items(command)))


def __get_units__(query_data):
    """
    Returns the chosen radius unit.
    If no unit is given, it will return degree.

    :param query_data: The query configuration data
    :type query_data: dict
    :return: The chosen unit
    :rtype:
    """
    if 'radius_unit' in query_data.keys():
        unit_string = query_data['radius_unit']
        if 'arcsec' in unit_string:
            return u.arcsec
        elif 'arcmin' in unit_string:
            return u.arcmin
        else:
            return u.deg
    else:
        return u.deg
