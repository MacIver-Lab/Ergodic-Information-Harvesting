# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:54:01 2017

@author: Chen Chen

"""
from ErgodicInfotaxisAPI.SimParameters import ErgodicParameters, EIDParameters
from ErgodicInfotaxisAPI.EntropyTrajectorySimulation import TrajEIDSim
from scipy.io import loadmat
from sys import argv
from numpy import linspace
from os import getcwd
from os import makedirs
from os.path import exists

def EIH_Sim(*argv):
    if len(argv) > 0:
        print(getcwd())
        # Parameter Objects
        ergParam = ErgodicParameters()
        eidParam = EIDParameters()

        #argv = argv[0]
        simType = argv[0]
        eidParam.randSeed = int(argv[1])

        eidParam.rawTraj = False
        eidParam.objCenter = 0.50
        eidParam.sInitPos = 0.40

        # Conditions to simulate
        ergParam.timeHorizon = 1.0
        eidParam.res = 101 # Spatial resolution
        ergParam.tRes = 101 # Time resolution
        ergParam.time = linspace(0.0, ergParam.timeHorizon, ergParam.tRes)
        ergParam.eidTime = linspace(0.0, ergParam.timeHorizon, eidParam.res)
        eidParam.maxT = 50
        ergParam.wControl = 20
        eidParam.pLSigmaAmp = 200
        eidParam.pLSigmaAmpBayesian = 200
        eidParam.pLSigmaAmpEID = 200
        eidParam.procNoiseSigma = 0.02
        eidParam.pLHistDepth = 1
        eidParam.Sigma = 0.06
        eidParam.objAmp = 0.20
        eidParam.blindIdx = False
        eidParam.blind = 'NA'
        eidParam.eidType = 'Entropy'
        eidParam.simType = 'EH'

        filename = 'EIH' + \
                '-SNR-' + str(eidParam.SNR)
                
        if simType == 'trackingOnly':
            # Parse input parameters
            sourceFile = argv[2]
            eidParam.fileName = argv[3]
            eidParam.saveDir = argv[4]
            if not exists(eidParam.saveDir):
                print("Save folder {0} does not exist, creating...".format(eidParam.saveDir))
                makedirs(eidParam.saveDir)

            # Load source file
            print('Loading source file {0}'.format(sourceFile))
            sourceData = loadmat(sourceFile)

            eidParam.sTraj = sourceData['sTrajList'].flatten()
            eidParam.rawTraj = sourceData['oTrajList'].flatten()
            dt = sourceData['dt'].flatten()[0]
            maxT = sourceData['maxTime'].flatten()[0]
            SNR = sourceData['SNR'].flatten()[0]
            wControl = sourceData['wControl'].flatten()[0]
            objAmp = sourceData['objAmp'].flatten()[0]
            Sigma = sourceData['Sigma'].flatten()[0]
            eidParam.attenMetrics = sourceData['AttenuateMetrics']

            # Conditions to simulate
            eidParam.maxT = float(maxT)
            eidParam.sInitPos = eidParam.sTraj[0]
            eidParam.SNR = float(SNR)
            eidParam.Sigma = float(Sigma)
            eidParam.objAmp = float(objAmp)
            ergParam.wControl = float(wControl)
            ergParam.dt = float(dt)
            eidParam.UpdateDeltaT(ergParam.dt)

            filename += '-gAttn-' + str(eidParam.attenMetrics['attenGain'].flatten()[0][0][0])
            print('Processing file {0} ...'.format(sourceFile))

            TrajEIDSim(ergParam, eidParam, eidParam.sTraj, showMsg=True)
            print('Attenuated trial simulated, output file - {0}'.format(eidParam.fileName))
    else:
        print(getcwd())
        print('Invalid input parameters! len={0}, param={1}'.format(len(argv),argv))

#param = 'trackingOnly 0 ../../SimData/SimTrial-001/AttenTraj/EIH-SNR-25-gAttn-5-RawTraj.mat EIH-SNR-25-gAttn-5-RandSeed-0.mat ../../SimData/SimTrial-001/SimTraj/'
#param = param.split()
#EIH_Sim(param)