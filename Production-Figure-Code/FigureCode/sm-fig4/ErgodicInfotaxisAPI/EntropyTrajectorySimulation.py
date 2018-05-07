# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:53:25 2017

@author: Chen Chen @ MacIver Lab
"""

# Import basic packages
import numpy as np
import matplotlib.pyplot as plt
# Import Ergodic packages
from ErgodicInfotaxisAPI.EID import EID
from ErgodicInfotaxisAPI.Ergodicity import Ergodicity
# Time
from time import strftime
# Save to .mat file
from ErgodicInfotaxisAPI.save2mat import save2mat
from IPython import display
from time import sleep
# Entropy
from scipy.stats import entropy


def TrajEIDSim(ergParam, eidParam, rawTraj, showMsg=True):
    # Reseed the numpy random module for deterministic behavior cross machines
    np.random.seed(eidParam.randSeed)
    
    # Initialize
    eid = EID(eidParam)
    erg = Ergodicity(ergParam)

    eidList = np.ones([ergParam.res, eidParam.maxIter])
    eidList[:, 1] /= eidList[:, 1].sum()
    sTrajList = np.array([eidParam.sInitPos])
    oTrajList = np.empty(0)
    phiList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
    pBList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
    ergList = np.empty(eidParam.maxIter)
    enpList = np.empty(eidParam.maxIter)
    pLast = np.ones([ergParam.res, 1])

    # Main Loop
    for it in range(eidParam.maxIter-1):
        # Read new trajectory segments
        traj = rawTraj[it*ergParam.res:(it+1)*ergParam.res]

        # Update measurement and EID with current trajectory
        eid.UpdateTraj(traj, pLast)
        phi, pB, objTraj = eid.Simulate(it)

        # Update belief and EID map
        eidList[:, it+1] = phi[:, -1]
        phiList[:, :, it+1] = phi
        pBList[:, :, it+1] = pB
        pLast = pB[:, -1]
        # Evaluate ergodicity
        ergList[it] = erg.computeErgMeasure(traj, phi[:, -1])
        # Evaluate Entropy
        enpList[it] = entropy(pLast)
        # Update trajectory list
        sTrajList = np.append(sTrajList, traj)
        oTrajList = np.append(oTrajList, objTraj)

        if showMsg:
            print('*** Iteration #{0} *** Erg = {1}, Enp = {2}'.format(
                it, ergList[it], enpList[it]))


    # Save data to MATLAB
    matDict = {}
    matDict['oTrajList'] = oTrajList
    matDict['sTrajList'] = sTrajList
    matDict['eidList'] = eidList
    matDict['pB'] = pBList
    matDict['phi'] = phiList
    matDict['ergList'] = ergList
    matDict['enpList'] = enpList
    matDict['AttenuateMetrics'] = eidParam.attenMetrics
    matDict['gAtten'] = eidParam.attenMetrics['attenGain'].flatten()[0][0][0]
    matDict['SNR'] = eidParam.SNR
    matDict['dt'] = eidParam.dt
    matDict['maxTime'] = eidParam.maxT
    matDict['wControl'] = ergParam.wControl
    matDict['sInitPos'] = eidParam.sInitPos
    matDict['objCenter'] = eidParam.objCenter
    matDict['simType'] = eidParam.simType
    matDict['eidType'] = eidParam.eidType
    matDict['randSeed'] = eidParam.randSeed
    matDict['pLSigmaAmp'] =  eidParam.pLSigmaAmp
    matDict['pLSigmaAmpEID'] = eidParam.pLSigmaAmpEID
    matDict['pLSigmaAmpBayesian'] = eidParam.pLSigmaAmpBayesian
    matDict['procNoiseSigma'] = eidParam.procNoiseSigma
    matDict['tRes'] = ergParam.tRes
    matDict['objAmp'] = eidParam.objAmp
    matDict['Sigma'] = eidParam.Sigma
    matDict['timeDone'] = strftime("%b-%d-%H-%M-%S")

    save2mat(eidParam.fileName, matDict, path=eidParam.saveDir, appendIdx=False)
