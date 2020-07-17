# -*- coding: utf-8 -*-

# Time
from time import strftime

# Entropy
from scipy.stats import entropy
import numpy as np
from numpy.random import default_rng
from scipy.interpolate import interp1d

# Import Ergodic packages
from ErgodicHarvestingLib.EID import EID
from ErgodicHarvestingLib.Ergodicity import Ergodicity
from ErgodicHarvestingLib.EID_Opt import ergoptimize
from ErgodicHarvestingLib.save2mat import save2mat


def EIDSim(ergParam, eidParam, showMsg=True):
    # Initialize the RNG with the provided seed
    rng = default_rng(eidParam.randSeed)
    # Initialize
    eid = EID(eidParam, rng)
    erg = Ergodicity(ergParam)
    # Initial control input is zero
    uinit = 0

    # Initialize
    if eidParam.simType == "EH":
        eidList = np.ones([ergParam.res, eidParam.maxIter])
        eidList[:, 1] /= eidList[:, 1].sum()
        phi = eidList[:, 1].reshape(len(eidList[:, 1]), 1)
        sTrajList = np.array([eidParam.sInitPos])
        oTrajList = np.empty(0)
        phiList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
        pBList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
        ergList = np.empty(eidParam.maxIter)
        enpList = np.empty(eidParam.maxIter)
        pLast = np.ones([ergParam.res, 1])
    elif eidParam.simType == "IF":
        eidParam.maxIter = np.int(np.round(eidParam.maxT / ergParam.dt))
        eidList = np.ones([eidParam.res, eidParam.maxIter])
        sTrajList = np.array([eidParam.sInitPos])
        oTrajList = np.empty(0)
        phiList = np.ones([eidParam.res, 1, eidParam.maxIter])
        pBList = np.ones([eidParam.res, 1, eidParam.maxIter])
        ergList = np.empty(0)
        enpList = np.empty(eidParam.maxIter)
        pLast = np.ones([eidParam.res, 1])
        # Workspace lattice
        wsLattice = np.linspace(0, 1, eidParam.res)
        pLast /= pLast.sum()
        phiList[:, :, 0] = eid.CalcEntropyEID(pLast)

    # Main Loop
    for it in range(eidParam.maxIter - 1):
        if eidParam.simType == "EH":
            # Use ergodic optimization to derive a ergodic trajectory
            [rawTraj, uinit] = ergoptimize(
                phi[:, -1],
                sTrajList[-1],
                control_init=uinit,
                ergParam=ergParam,
                showStats=False,
                showMsg=showMsg,
            )
            # [Optional] Filter the wiggle to evaluate causation
            traj = interp1d(ergParam.time, rawTraj)(ergParam.eidTime)
            # Update measurement and EID with current trajectory
            eid.UpdateTraj(traj, pLast)
            phi, pB, objTraj = eid.Simulate(it)
            # Update belief and EID map
            eidList[:, it] = phi[:, -1]
            phiList[:, :, it] = phi
            pBList[:, :, it] = pB
            pLast = pB[:, -1]
            # Evaluate ergodicity
            ergList[it] = erg.computeErgMeasure(traj, pLast)
        elif eidParam.simType == "IF":
            # Plan for next move
            traj = eid.PlanNextMove(
                phiList[:, :, it],
                sTrajList[it],
                wsLattice,
                stepSize=eidParam.stepSize,
                planDepth=eidParam.planDepth,
            )
            # Update measurement and EID with current move
            eid.UpdateTraj(traj, pLast)
            phi, pB, objTraj = eid.Simulate(it)
            # Update belief and EID map
            eidList[:, it + 1] = phi[:, -1]
            phiList[:, :, it + 1] = phi
            pBList[:, :, it + 1] = pB
            pLast = pB[:, -1]
            # Evaluate Entropy
            enpList[it] = entropy(pLast)

        # Update trajectory list
        sTrajList = np.append(sTrajList, traj)
        oTrajList = np.append(oTrajList, objTraj)

    # Save data in MATLAB mat format
    matDict = {}
    matDict["oTrajList"] = oTrajList
    matDict["sTrajList"] = sTrajList
    matDict["gAtten"] = 0
    matDict["timeHorizon"] = ergParam.timeHorizon
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
    matDict["eidList"] = eidList
    matDict["pB"] = pBList
    matDict["phi"] = phiList
    matDict["blindIdx"] = eidParam.blindIdx
    matDict["blindType"] = eidParam.blind
    matDict["timeDone"] = strftime("%b-%d-%H-%M-%S")
    matDict["enpList"] = enpList
    if eidParam.simType == "EH":
        matDict["ergList"] = ergList
    elif eidParam.simType == "IF":
        matDict["StepSize"] = eidParam.stepSize
        matDict["planDepth"] = eidParam.planDepth

    save2mat(eidParam.filename, matDict, path=eidParam.saveDir, appendIdx=False)
