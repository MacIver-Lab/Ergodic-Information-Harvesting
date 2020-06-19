# -*- coding: utf-8 -*-

# Import required libraries
import numpy as np
from scipy.stats import entropy

from ErgodicHarvestingLib.cyEntropy import cyNormPDF, cyNormPDF2, evalGaussianConst


def EntropyEID(prior, wsRes, sigmaM, sigmaL):
    """Computes Entropy EID map
    
    EntropyEID function will output predicted Entropy-based EID map given 
    input parameters. The Entropy-based EID is the negative of entropy gain 
    over the entropy of prior. This functions assumes a workspace of [0, 1]
    
    Args:
        mmLambda (lambda): Lambda function of measurement model
            Measurement model in the form of a Lambda function, 
            e.g. mmLambda = lambda x: scipy.stats.norm.pdf(samples, x, sigma) 
            for a symmetric Gaussian measurement model
        
        prior (numpy.array): Prior distribution
        wsRes (int): Spatial resolution of the workspace
        sigmaM (float): Sigma (variance) of the measurement model
        sigmaL (float): Sigma (variance) of the likelihood function
        
            
    """
    ## Figure out workspace spatial resolution and generate spaial samples
    wsSamples = np.linspace(0.0, 1.0, wsRes)

    ## Measurement samples
    mmSamples = wsSamples

    ## Make sure prior is a probability distribution (integral = 1)
    prior /= prior.sum()
    prior = prior.flatten()

    # Compute base entropy to compare
    sBase = entropy(prior)

    # Hashtable for norm pdf
    mmDict = {}

    # measurement modal constant
    mmConst1, mmConst2 = evalGaussianConst(sigmaM)
    pLConst1, pLConst2 = evalGaussianConst(sigmaL)

    ## Create Measurement model object
    def upsilon(theta=0.5):
        if type(theta) is float:
            mm = mmDict.get(theta, None)
            if mm is None:
                mm = cyNormPDF(wsSamples, theta, mmConst1, mmConst2, nrm=1)
                mmDict[theta] = mm
        else:
            # Queries multiple possible target thetas
            # Check theta level key hit
            key = (theta[0], theta[-1])
            mm = mmDict.get(key, None)
            if mm is None:
                mm = []
                for t in theta:
                    # Check t level key hit
                    mm_ = mmDict.get(t, None)
                    if mm_ is None:
                        mm_ = cyNormPDF(wsSamples, t, mmConst1, mmConst2, nrm=1)
                        mmDict[t] = mm_
                    mm.append(mm_)
                mm = np.array(mm)
                mmDict[key] = mm
        return mm

    def Pvtheta(v, theta):
        # Compute likelihood
        pvt = cyNormPDF2(upsilon(theta), v, pLConst1, pLConst2)
        return pvt

    def Pvx(v):
        # Likelihood function
        pL = Pvtheta(v, wsSamples)
        # Measurement probability
        pvx = pL
        # Loop through theta samples (all possible target location, weighted by
        # prior belief)
        for idx in range(pvx.shape[0]):
            pvx[idx, :] = pvx[idx, :] * prior[idx]
        # Total probability of measurement
        pvx = pvx.sum(axis=0)

        # OUTPUT:
        #   pL : Individual likelihood function of all the possible source location
        #   pvx: Overall likehood of measurement at every location x
        return pvx, pL

    """ Loop through all of the possible measurement states for
        posterior belief and entropy
    """
    deltaS = np.zeros(wsRes)
    for mms in mmSamples:
        # Posterior
        pvx, pL = Pvx(mms)
        pp = pL * prior
        # Entropy (already negated)
        dS = sBase - entropy(pp.T)
        # Remove nan values
        dS[(np.isnan(dS) | np.isinf(dS))] = 0
        # Update result
        deltaS += pvx * dS

    # Convert deltaS to Entropy EID
    eid = deltaS
    eid[np.where(eid < 0)] = 0
    eid /= eid.sum()
    return eid
