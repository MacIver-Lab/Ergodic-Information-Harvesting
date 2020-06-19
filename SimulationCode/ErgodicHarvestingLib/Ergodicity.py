# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import trapz, quad


class Ergodicity(object):
    def __init__(self, ergParam):
        self.Nfourier = ergParam.nFourier
        self.dimw = 1  # workspace dimension
        wlimit = float(1)
        self.wlimit = float(1)
        self.res = ergParam.res  # Spatial Resolution
        self.time = np.linspace(0.0, 1.0, self.res)
        self.xlist = np.linspace(0.0, 1.0, self.res)
        # set up a grid over the frequency
        klist = np.arange(self.Nfourier)
        # set up hk
        s = (float(self.dimw) + 1.0) / 2.0
        self.Lambdak = 1.0 / (1.0 + klist ** 2) ** s
        self.klist = klist / self.wlimit * np.pi
        self.hk = np.zeros(self.Nfourier).flatten()
        for index in range(self.Nfourier):
            integ = quad(
                lambda x: (np.cos(x * self.klist[index])) ** 2, 0, float(wlimit)
            )
            self.hk[index] = np.sqrt(integ[0])

    def normalize_pdf(self):
        # function to normalize a pdf
        self.pdf /= np.sum(self.pdf) / np.product(self.pdf.shape)

    def calculate_uk(self, pdf):
        # calculate Fourier coefficients of the distribution
        self.uk = np.zeros(self.Nfourier).flatten()
        for index in range(len(self.uk)):
            uk_interior = pdf / self.hk[index]
            basis_part = np.cos(self.klist[index] * self.xlist)
            uk_interior *= self.wlimit / self.res * basis_part
            self.uk[index] = np.sum(uk_interior)

    def setPDF(self, pdf):
        # input pdf
        self.pdf = pdf
        self.normalize_pdf()
        self.calculate_uk(self.pdf)
        pass

    def ergodicEval(self):
        # evaluate the ergodic metric (ck, uk, need to be calculated)
        self.erg = np.sum(self.Lambdak * (self.ck - self.uk) ** 2)
        return self.erg

    def ckEval(self):
        X = self.X_current
        time = self.time
        T = time[-1]
        # change coordinates from configuration to ergodic workspace
        W = X.flatten()
        self.ck = np.zeros(self.Nfourier).flatten()
        for index in range(len(self.ck)):
            ck_interior = 1.0 / (self.hk[index] * T)
            basis_part = np.cos(self.klist[index] * W)
            ck_interior = ck_interior * basis_part
            self.ck[index] = trapz(ck_interior, time)

    def computeErgMeasure(self, x, pdf):
        self.X_current = x
        self.setPDF(pdf)
        self.ckEval()
        return self.ergodicEval()
