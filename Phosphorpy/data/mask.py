from datetime import datetime


class MaskLog:
    _mask = None
    _label = ''
    _date = None

    def __init__(self, mask, label):
        self._mask = mask
        self._label = label
        self._date = datetime.now()

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        raise AttributeError('Attribute setting is not allowed')

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value


class Mask:
    _mask_log = None
    _mask = None

    def __init__(self):
        """
        The Mask-class provides the basic mask system.
        """
        self._mask_log = []

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask_log.append(MaskLog(value, 'new {}'.format(len(self._mask_log)+1)))
        self._mask = value

    def __getitem__(self, item):
        for log in self._mask_log:
            if item == log.label:
                return log
        raise KeyError('Key not found!')

    def __setitem__(self, key, value):
        self._mask_log.append(MaskLog(key, value))
