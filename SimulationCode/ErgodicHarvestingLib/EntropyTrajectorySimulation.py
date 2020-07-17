# -*- coding: utf-8 -*-

# Import basic packages
import numpy as np
from numpy.random import default_rng

# Time
from time import strftime

# Entropy
from scipy.stats import entropy

# Import Ergodic packages
from ErgodicHarvestingLib.EID import EID
from ErgodicHarvestingLib.Ergodicity import Ergodicity
from ErgodicHarvestingLib.save2mat import save2mat


def TrajEIDSim(ergParam, eidParam, rawTraj, showMsg=True):
    # Initialize the RNG with the provided seed
    rng = default_rng(eidParam.randSeed)
    # Initialize
    eid = EID(eidParam, rng)
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
    for it in range(eidParam.maxIter - 1):
        # Read new trajectory segments
        traj = rawTraj[it * ergParam.res : (it + 1) * ergParam.res]
        # Update measurement and EID with current trajectory
        eid.UpdateTraj(traj, pLast)
        phi, pB, objTraj = eid.Simulate(it)
        # Update belief and EID map
        eidList[:, it + 1] = phi[:, -1]
        phiList[:, :, it + 1] = phi
        pBList[:, :, it + 1] = pB
        pLast = pB[:, -1]
        # Evaluate ergodicity
        if phi.size >= ergParam.res and phi.shape[1] > 0:
            ergList[it] = erg.computeErgMeasure(traj, phi[:, -1])
        # Evaluate Entropy
        if len(pLast) > 1:
            enpList[it] = entropy(pLast)
        # Update trajectory list
        sTrajList = np.append(sTrajList, traj)
        oTrajList = np.append(oTrajList, objTraj)
        if showMsg:
            print(
                "*** Iteration #{0} *** Erg = {1}, Enp = {2}".format(
                    it, ergList[it], enpList[it]
                )
            )

    # Save data to MATLAB
    matDict = {}
    matDict["oTrajList"] = oTrajList
    matDict["sTrajList"] = sTrajList
    matDict["eidList"] = eidList
    matDict["pB"] = pBList
    matDict["phi"] = phiList
    matDict["ergList"] = ergList
    matDict["enpList"] = enpList
    matDict["AttenuateMetrics"] = eidParam.attenMetrics
    matDict["gAtten"] = eidParam.attenMetrics["attenGain"].flatten()[0][0][0]
    matDict["SNR"] = eidParam.SNR
    matDict["dt"] = eidParam.dt
    matDict["maxTime"] = eidParam.maxT
    matDict["wControl"] = ergParam.wControl
    matDict["sInitPos"] = eidParam.sInitPos
    matDict["objCenter"] = eidParam.objCenter
    matDict["simType"] = eidParam.simType
    matDict["eidType"] = eidParam.eidType
    matDict["randSeed"] = eidParam.randSeed
    matDict["pLSigmaAmp"] = eidParam.pLSigmaAmp
    matDict["pLSigmaAmpEID"] = eidParam.pLSigmaAmpEID
    matDict["pLSigmaAmpBayesian"] = eidParam.pLSigmaAmpBayesian
    matDict["procNoiseSigma"] = eidParam.procNoiseSigma
    matDict["tRes"] = ergParam.tRes
    matDict["objAmp"] = eidParam.objAmp
    matDict["Sigma"] = eidParam.Sigma
    matDict["timeDone"] = strftime("%b-%d-%H-%M-%S")

    save2mat(eidParam.fileName, matDict, path=eidParam.saveDir, appendIdx=False)
