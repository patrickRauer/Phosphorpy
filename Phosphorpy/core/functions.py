import numpy as np
import numba


def power_2_10(x):
    return np.power(10., -x/2.5)


@numba.vectorize
def subtract(a, b):
    """
    Numba implementation of a subtraction

    :param a:
    :param b:
    :return:
    """
    return a-b


def gaus(x, a, b, c, d, e):
    """
    Fitting function with a gaussian component, a linear component and a static component


    .. math:

        f(\lambda) = a*e^{-\frac{(\lambda-lambda_0)**2}{2*\sigma**2}}+b*\lambda+c


    :param x: The variable values
    :param a: The strength of the gaussian
    :param b: The shift of the center
    :param c: The sigma of the gaussian
    :param d: The steepness rate of the linear component
    :param e: The value of the static component
    :return: The corresponding function values
    :type: np.ndarray
    """
    return a * np.exp(-np.square(x-b)/(2*c**2))+e+d*x