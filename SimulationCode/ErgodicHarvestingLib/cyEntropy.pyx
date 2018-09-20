# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 10:57:38 2018

@author: chenc
"""

cimport cython
import numpy as np
cimport numpy as np
from scipy.stats import norm, entropy
from libc.math cimport exp, sqrt, log

@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def cyNormPDF(np.ndarray[np.double_t, ndim=1] x, double m, double sigma, int nrm=0):
    cdef double leftPart = 1.0 / sqrt(2.0*np.pi*sigma*sigma)
    x = leftPart * np.exp(-(x - m)**2 / (2*sigma*sigma))
    if nrm == 0:
        x /= x.sum()
    elif nrm == 1:
        x /= x.max()
    return x