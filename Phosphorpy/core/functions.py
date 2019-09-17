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

    :math:`f(x) = a*\exp\left(-\\frac{(x-x_0)^2}{2*\sigma^2}\\right)+b*x+c`

    :param x: The variable values
    :param a: The strength of the gaussian
    :param b: The shift of the center
    :param c: The sigma of the gaussian
    :param d: The steepness rate of the linear component
    :param e: The value of the static component
    :return: The corresponding function values
    :rtype: np.ndarray
    """
    return a * np.exp(-np.square(x-b)/(2*c**2))+e+d*x


def smooth2d(mat, c=5):
    """
    Smooths a 2d-array

    :param mat: The input data
    :type mat: np.ndarray
    :param c: The number of smooths
    :param c: int
    :return: The c-times smoothed input data
    """
    if c == 0:
        return mat
    out = np.zeros(mat.shape)

    out[0, 0] = (2*mat[0, 0]+mat[1, 0]+mat[0, 1])/4
    out[0, -1] = (2*mat[0, -1]+mat[1, -1]+mat[0, -2])/4
    out[-1, 0] = (2*mat[-1, 0]+mat[-2, 0]+mat[-1, 1])/4
    out[-1, -1] = (2*mat[-1, -1]+mat[-2, -1]+mat[-1, -2])/4

    out[1:-1, 0] = (3*mat[1:-1, 0]+mat[1:-1, 1]+mat[:-2, 0]+mat[2:, 0])/6
    out[1:-1, -1] = (3*mat[1:-1, -1]+mat[1:-1, -2]+mat[:-2, -1]+mat[2:, -1])/6
    out[0, 1:-1] = (3*mat[0, 1:-1]+mat[1, 1:-1]+mat[0, :-2]+mat[0, 2:])/6
    out[-1, 1:-1] = (3*mat[-1, 1:-1]+mat[-2, 1:-1]+mat[-1, :-2]+mat[-1, 2:])/6

    out[1:-1, 1:-1] = (4*mat[1:-1, 1:-1]+mat[:-2, 1:-1]+mat[2:, 1:-1]+mat[1:-1, :-2]+mat[1:-1, 2:])/8
    return smooth2d(out, c=c-1)
