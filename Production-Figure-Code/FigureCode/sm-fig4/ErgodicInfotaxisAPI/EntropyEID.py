# -*- coding: utf-8 -*-
"""
Entropy EID Algorithm

Computes Entropy-based EID map given measurement model and posterior belief

Created on Wed Jan 17 13:32:35 2018

@author: Chen Chen
"""

# Import required libraries
import numpy as np
from scipy.stats import norm, entropy

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
    
    ## Create Measurement model object
    def upsilon(theta=0.5):
        if type(theta) is float:
            mm = norm.pdf(wsSamples,theta,sigmaM)
            mm /= max(mm)
        else:
            # Queries multiple possible target thetas
            mm = []
            for t in theta:
                mm_ = norm.pdf(wsSamples,t,sigmaM)
                mm_ /= max(mm_)
                mm.append(mm_)
            mm = np.array(mm)
        return mm
    
    # Measurement likelihood lambda function
    #Pvtheta = lambda v, theta: ( 1.0 / (np.sqrt(2.0*np.pi) * sigmaL) ) * \
    #          np.exp( -(v - upsilon(theta))**2 / (2.0 * sigmaL**2) )
    pLConst1 = ( 1.0 / (np.sqrt(2.0*np.pi) * sigmaL) )
    pLConst2 = -1 / (2.0 * sigmaL**2) 
    def Pvtheta(v, theta):
        # Compute likelihood
        pvt = pLConst1 * np.exp( (v - upsilon(theta))**2 * pLConst2 )
          
        # Normalize
        if pvt.ndim == 1:
            return pvt / pvt.sum()
        else:
            for idx in range(pvt.shape[0]):
                pvt[idx, :] /= pvt[idx, :].sum()
          
        return pvt
        
    # Measurement probability Pvx over all possible sensor location    
    def Pvx(v):
        # Likelihood function
        pL = Pvtheta(v, wsSamples)
        # Measurement probability
        pvx = pL
          
        # Loop through theta samples (all possible target location, weighted by 
        # prior belief)
        for idx in range(pvx.shape[0]):
            pL[idx, :] /= pL[idx, :].sum()
            pvx[idx, :] = pvx[idx, :] * prior[idx]

        # Average over all the weighted target location samples to get the 
        # overall probability of getting measurement vx
        pvx = pvx.mean(axis=0)
        pvx /= pvx.sum()
          
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

        # Entropy
        dS = entropy(pp.T) - sBase
        
        # Remove nan values
        dS[dS == np.nan] = 0
        
        # Update result
        deltaS += pvx * dS

    # Convert deltaS to Entropy EID
    eid = -deltaS
    
    # Remove invalid (infinity) numbers
    eid[eid == np.inf] = 0
    
    # Remove negative numbers
    eid[eid < 0] = 0
    #if eid.min() < 0:
    #    eid -= eid.min()
    
    #print("Belief Entropy = {0}, ".format(sBase))
    eid /= eid.sum()
    return eid

""" Unit test
wsSamples = np.linspace(0,1,101)
priorP = norm.pdf(wsSamples,0.3,0.1) #+ norm.pdf(wsSamples,0.6,0.2)
priorP /= priorP.sum()
EID_Entropy = EntropyEID(priorP, 101, 0.1, 0.1)
plt.plot(wsSamples, priorP)
plt.plot(wsSamples, EID_Entropy)
#"""