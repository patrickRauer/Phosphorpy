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
