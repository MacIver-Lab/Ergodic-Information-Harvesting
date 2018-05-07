import numpy as np

def DistErgodicity(pdf):
    # set up a grid over the frequency
    klist = np.arange(Nfourier)

    # do some ergodic stuff
    s = (float(dimw) + 1.0) / 2.0
    Lambdak = 1.0 / (1.0 + klist**2)**s
    klist = klist / wlimit * np.pi
    hk = np.zeros(Nfourier).flatten()
    for index in range(Nfourier):
        integ = quad(lambda x: (np.cos(x*klist[index]))**2, 0.0, float(wlimit))
        hk[index] = np.sqrt(integ[0])

    def calculate_uk(self, pdf):
        # calculate Fourier coefficients of the distribution
        #pdf = pdf.T
        uk = np.zeros(Nfourier).flatten()
        for index in range(len(uk)):
            uk_interior = pdf / hk[index]
            basis_part = np.cos(klist[index] * xlist)
            uk_interior *= wlimit / res * basis_part
            uk[index] = np.sum(uk_interior)
        
        #print("uk", uk[0])
    
    def ckeval(self):
        X = X_current
        time = time
        T = time[-1]
        # change coordinates from configuration to ergodic workspace
        W = X.flatten()
        ck = np.zeros(Nfourier).flatten()
        #xlist = tj.T
        for index in range(len(ck)):
            ck_interior = 1.0 / hk[index] * 1.0 / (float(T))
            basis_part = np.cos(klist[index] * W)
            ck_interior = ck_interior * basis_part
            ck[index] = trapz(ck_interior,time) #np.sum(dt*ck_interior)
