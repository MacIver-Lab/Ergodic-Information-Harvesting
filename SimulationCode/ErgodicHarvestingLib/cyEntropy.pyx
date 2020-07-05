# -*- coding: utf-8 -*-

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt, log, pi


@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping
def evalGaussianConst(double sigma):
    gConst1 = ( 1.0 / (sqrt(2.0*pi) * sigma) )
    gConst2 = -1.0 / (2.0 * sigma * sigma)
    return gConst1, gConst2


@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping
def cyNormPDF(np.ndarray[np.double_t, ndim=1] x, double m, double const1, double const2, int nrm=0):
    x = const1 * np.exp( (x - m)**2 * const2 )
    if nrm == 0:
        x /= x.sum()
    elif nrm == 1:
        x /= x.max()
    return x


@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping
def cyNormPDF2(np.ndarray[np.double_t, ndim=2] x, double m, double const1, double const2):
    x = const1 * np.exp( (x - m)**2 * const2 )
    x = (x.T / x.sum(axis=0)).T
    return x