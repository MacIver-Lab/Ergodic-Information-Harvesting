# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 10:50:25 2018

@author: chenc
"""
cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport exp, sin, expm1

@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def rk4c(f, np.ndarray t_in, double y0, int tScale=1):
    cdef int n, idx
    cdef double dt0, dt, t0, t1
    cdef np.ndarray vt, vy
    
    n = len(t_in)
    dt0 = t_in[1] - t_in[0]
    dt = dt0 / float(tScale)
    t0 = t_in[0]
    t1 = t_in[n-1]
    vt = vy = np.empty(n)
    vt[0] = t = t0
    vy[0] = y = y0
    t += dt
    idx = 1
    while t < t1:
        k1 = dt * f(t, y)
        k2 = dt * f(t + 0.5 * dt, y + 0.5 * k1)
        k3 = dt * f(t + 0.5 * dt, y + 0.5 * k2)
        k4 = dt * f(t + dt, y + k3)
        if idx % tScale == 0:
            i = int(idx / tScale)
            vt[i] = t
            vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6.0
        t += dt
        idx += 1
    return vt, vy

@cython.profile(False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def rk2(f, np.ndarray t_in, double y0, int tScale=1):
    cdef int n, idx
    cdef double dt0, dt, t0, t1, coef
    cdef np.ndarray vt, vy
    
    n = len(t_in)
    dt0 = t_in[1] - t_in[0]
    dt = dt0 / float(tScale)
    t0 = t_in[0]
    t1 = t_in[n-1]
    vt = vy = np.empty(n)
    vt[0] = t = t0
    vy[0] = y = y0
    idx = 1
    coef = 2.0 / 3.0
    t += dt
    while t < t1:
        k1 = dt * f(t, y)
        k2 = dt * f(t + coef * dt, y + coef * k1) * 3.0
        if idx % tScale == 0:
            i = int(idx / tScale)
            vt[i] = t
            vy[i] = y = y + (k1 + k2) / 4.0
        t += dt
        idx += 1
    return vt, vy
