# -*- coding: utf-8 -*-

from scipy.io import loadmat
from numpy import linspace
from os import getcwd, makedirs
from os.path import exists

from ErgodicHarvestingLib.SimParameters import ErgodicParameters, EIDParameters
from ErgodicHarvestingLib.EntropyTrajectorySimulation import TrajEIDSim


def EIH_Sim(*args):
    import time
    from numpy.random import rand
    time.sleep(0.5 * rand())
    return

    if len(args) > 0:
        # Parameter Objects
        ergParam = ErgodicParameters()
        eidParam = EIDParameters()
        simType = args[0]
        eidParam.randSeed = int(args[1])
        eidParam.rawTraj = False
        eidParam.objCenter = 0.50
        eidParam.sInitPos = 0.40
        # Conditions to simulate
        ergParam.timeHorizon = 1.0
        eidParam.res = 101  # Spatial resolution
        ergParam.tRes = 101  # Time resolution
        ergParam.time = linspace(0.0, ergParam.timeHorizon, ergParam.tRes)
        ergParam.eidTime = linspace(0.0, ergParam.timeHorizon, eidParam.res)
        eidParam.maxT = 50
        ergParam.wControl = 10
        eidParam.pLSigmaAmp = 200
        eidParam.pLSigmaAmpBayesian = 200
        eidParam.pLSigmaAmpEID = 200
        eidParam.procNoiseSigma = 0.01
        eidParam.pLHistDepth = 1
        eidParam.Sigma = 0.06
        eidParam.objAmp = 0.20
        eidParam.blindIdx = False
        eidParam.blind = "NA"
        eidParam.eidType = "Entropy"
        eidParam.simType = "EH"

        filename = f"EIH-SNR-{eidParam.SNR}"
        if simType == "trackingOnly":
            # Parse input parameters
            sourceFile = args[2]
            eidParam.fileName = args[3]
            eidParam.saveDir = args[4]
            if not exists(eidParam.saveDir):
                print(
                    "Save folder {0} does not exist, creating...".format(
                        eidParam.saveDir
                    )
                )
                makedirs(eidParam.saveDir)
            # Load source file
            print("Loading source file {0}".format(sourceFile))
            sourceData = loadmat(sourceFile)
            eidParam.trajSim = True
            eidParam.sTraj = sourceData["sTrajList"].flatten()
            eidParam.rawTraj = sourceData["oTrajList"].flatten()
            dt = sourceData["dt"].flatten()[0]
            maxT = sourceData["maxTime"].flatten()[0]
            SNR = sourceData["SNR"].flatten()[0]
            wControl = sourceData["wControl"].flatten()[0]
            objAmp = sourceData["objAmp"].flatten()[0]
            Sigma = sourceData["Sigma"].flatten()[0]
            eidParam.attenMetrics = sourceData["AttenuateMetrics"]
            # Conditions to simulate
            eidParam.maxT = float(maxT)
            eidParam.sInitPos = eidParam.sTraj[0]
            eidParam.SNR = float(SNR)
            eidParam.Sigma = float(Sigma)
            eidParam.objAmp = float(objAmp)
            ergParam.wControl = float(wControl)
            ergParam.dt = float(dt)
            eidParam.UpdateDeltaT(ergParam.dt)
            filename += "-gAttn-" + str(
                eidParam.attenMetrics["attenGain"].flatten()[0][0][0]
            )
            print("Processing file {0} ...".format(sourceFile))
            TrajEIDSim(ergParam, eidParam, eidParam.sTraj, showMsg=True)
            print(
                "Attenuated trial simulated, output file - {0}".format(
                    eidParam.fileName
                )
            )
    else:
        print(getcwd())
        print("Invalid input parameters! len={0}, param={1}".format(len(args), args))
