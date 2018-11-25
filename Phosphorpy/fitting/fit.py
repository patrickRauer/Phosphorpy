import numpy as np


class Fit:
    _x_log = False
    _y_log = False
    _func = None
    _fit_value = None

    def __init__(self, func, x_log=False, y_log=False):
        self.func = func
        self.x_log = x_log
        self.y_log = y_log

    def __call__(self, *args, **kwargs):
        x = args[0]
        y = args[1]
        if self.x_log:
            x = np.log10(x)

        if self.y_log:
            y = np.log10(y)

        fit = self.func(x, y)
        self._fit_value = fit
        return fit

    @property
    def x_log(self):
        return self._x_log

    @x_log.setter
    def x_log(self, value):
        if type(value) == bool:
            self._x_log = value
        else:
            raise ValueError('x_log must be a boolean')

    @property
    def y_log(self):
        return self._y_log

    @y_log.setter
    def y_log(self, value):
        if type(value) == bool:
            self._y_log = value
        else:
            raise ValueError('y_log must be a boolean')

    @property
    def func(self):
        return self._func

    @func.setter
    def func(self, f):
        if callable(f):
            self._func = f
        else:
            raise ValueError('Value must be callable.')

    @property
    def fit_value(self):
        return self._fit_value

    @fit_value.setter
    def fit_value(self, value):
        raise AttributeError('Attribute can not be set')


def poly(degree):
    def poly_fit(x, y):
        return np.polyfit(x, y, degree)
    return poly_fit


class PolyFit(Fit):

    def __init__(self, degree, **kwargs):
        Fit.__init__(self, poly(degree), **kwargs)
