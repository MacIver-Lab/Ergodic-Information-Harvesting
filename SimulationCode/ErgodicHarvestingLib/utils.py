import numpy as np
from scipy.stats import norm
from functools import reduce


def matmult(*x):
    """
    Shortcut for standard matrix multiplication.
    matmult(A,B,C) returns A*B*C.
    """
    return reduce(np.dot, x)


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like 1-D arrays to form the cartesian product of.
    out : ndarray Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)
    m = n / arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j * m : (j + 1) * m, 1:] = out[0:m, 1:]
    return out


def bimodal_gaussian_pdf(X):
    # define a Gaussian PDF over grid X
    nrmPDF = norm.pdf(X, 0.8, 0.2)
    nrmPDF /= nrmPDF.sum()
    return nrmPDF


def print_cyan(msg):
    print(f"\033[96m{msg}\033[00m")


def print_green(msg):
    print(f"\033[92m{msg}\033[00m")


def print_yellow(msg):
    print(f"\033[93m{msg}\033[00m")


def print_color(msg, color=None):
    if not isinstance(color, str):
        print(msg)
        return
    if "cyan" in color:
        print_cyan(msg)
    elif "green" in color:
        print_green(msg)
    elif "yellow" in color:
        print_yellow(msg)
