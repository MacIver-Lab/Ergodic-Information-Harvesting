from itertools import product
from multiprocessing import Pool, cpu_count
from ErgodicHarvestingLib.SimParameters import ErgodicParameters, EIDParameters
from ErgodicHarvestingLib.Simulation import EIDSim
from ErgodicHarvestingLib.ParameterIO import loadParams
from os import makedirs
from os.path import exists
from numpy import linspace
import numpy as np
from scipy.io import loadmat


# LoadMothData() load and normalize moth tracking data
# for EIH simulations
# Returns [sensor, target] trajectory
def flatten(x):
    out = []
    for r in x[0]:
        out.append(r[0])
    return np.array(out)

def normalize(s, t, mid=0.5, gain=0.1):
    # Remove offset
    sMean = np.mean(s)
    s -= np.mean(s)
    t -= sMean
    # Scale
    sScale = np.max(s) * gain**-1
    s /= sScale
    t /= sScale
    # Offset
    s += mid
    t += mid
    return s, t

def loadMothData(target="M300lux", trialID=0, nrmMid=0.5, nrmGain=0.1):
    traj = flatten(loadmat('./MothData.mat')[f'trial_{target}'])[trialID,:,:]
    s, t = normalize(traj[:, 0], traj[:, 1], nrmMid, nrmGain)
    return [s, t]

def SimulationMain(dataFiles, nThread=1):
    if type(dataFiles) is not list:
        dataFiles = [dataFiles]
    paramList = []
    simParamList = []
    nSimJobsList = []
    for dat in dataFiles:
        # Load parameters
        paramList.append(loadParams(dat))
        # Permutate conditions
        simParamList.append(
            list(
                product(
                    paramList[-1]['SNR'],
                    paramList[-1]['procNoiseSigma'],
                    paramList[-1]['pLSigmaAmp'],
                    paramList[-1]['Sigma'],
                    paramList[-1]['objAmp'],
                    paramList[-1]['dt'],
                    paramList[-1]['wControl'],
                    paramList[-1]['randSeed']
                )
            )
        )
        nSimJobsList.append(len(simParamList[-1]))

    nTrials = len(nSimJobsList)
    nTotalJobs = sum(nSimJobsList)
    # Limit pool size when job size is smaller than total available threads
    if nTotalJobs < nThread:
        nThread = nTotalJobs
    # Start a new parallel pool
    print("Starting parallel pool with {0} threads".format(nThread))
    pool = Pool(processes=nThread)
    jobs = []
    simCnt = 0
    for trial in range(nTrials):
        # Parse parameters
        param = paramList[trial]
        simParam = simParamList[trial]
        nJobs = nSimJobsList[trial]
        eidParam = param['eidParam']
        ergParam = param['ergParam']
        filename = param['filename']
        figName = param['figName']
        # Check if saveDir exists, create the folder if not
        if not exists(eidParam.saveDir):
            print("Save folder {0} does not exist, creating...".format(
                eidParam.saveDir))
            makedirs(eidParam.saveDir)

        for it in range(nJobs):
            eidParam.SNR = simParam[it][0]
            eidParam.procNoiseSigma = simParam[it][1]
            eidParam.pLSigmaAmp = simParam[it][2]
            eidParam.pLSigmaAmpBayesian = simParam[it][2]
            eidParam.pLSigmaAmpEID = simParam[it][2]
            eidParam.Sigma = simParam[it][3]
            eidParam.objAmp = simParam[it][4]
            if type(eidParam.rawTraj) is str and 'moth' in eidParam.rawTraj:
                eidParam.rawTraj = np.array(loadMothData(
                    target="M300lux", trialID=0, nrmGain=eidParam.objAmp)[1])
            ergParam.dt = simParam[it][5]
            eidParam.UpdateDeltaT(simParam[it][5])
            if eidParam.simType == 'IF':
                eidParam.maxIter = round(eidParam.maxT / ergParam.dt)
            ergParam.wControl = simParam[it][6]
            eidParam.randSeed = simParam[it][7]
            ergParam.time = linspace(0.0, ergParam.timeHorizon, ergParam.tRes)
            ergParam.eidTime = linspace(
                0.0, ergParam.timeHorizon, eidParam.res)
            eidParam.filename = (
                filename
                .replace('SNR', 'SNR-' + str(eidParam.SNR))
                .replace('RandSeed', 'RandSeed-' + str(eidParam.randSeed))
            )

            # Wait for all the thread to finish if the pool is full
            if len(jobs) == nThread:
                print("All thread used, waiting until finish... ...")
                # Wait until all the active thread to finish
                for job in jobs:
                    job.join()
                # Clear the pool
                jobs = []

            # Start a new job thread
            p = pool.Process(target=EIDSim, args=(
                ergParam, eidParam, False, False, False, True))
            p.start()
            jobs.append(p)
            simCnt += 1
            print("Starting New Simulation Thread ({0:d}/{1:d}), progress {2:.2f}% ({3:d}/{4:d})".format(
                len(jobs), nThread,
                100.0*(simCnt)/nTotalJobs, simCnt, nTotalJobs))

    print("Waiting for remaining jobs to finish... ...")
    # Wait until all the active thread to finish
    for job in jobs:
        job.join()
