from itertools import product 
from multiprocessing import Pool, cpu_count
from ErgodicHarvestingLib.SimParameters import ErgodicParameters, EIDParameters
from ErgodicHarvestingLib.Simulation import EIDSim
from ErgodicHarvestingLib.ParameterIO import loadParams
from os import makedirs
from os.path import exists
from numpy import linspace

def SimulationMain(dataFile, nThread=1):
    # Load parameters 
    param = loadParams(dataFile)
    for key in param:
        globals()[key] = param[key]
    
    eidParam.UpdateDeltaT(ergParam.dt)
    # Check if saveDir exists, create the folder if not
    if not exists(eidParam.saveDir):
        print("Save folder {0} does not exist, creating...".format(eidParam.saveDir))
        makedirs(eidParam.saveDir)

    # Permutate conditions
    simParam = list(product(SNR,procNoiseSigma,pLSigmaAmp,Sigma,objAmp,dt,wControl,randSeed))
    nSimTrials = len(simParam)
    
    # Limit pool size when job size is smaller than total available threads 
    if nSimTrials < nThread:
        nThread = nSimTrials
    # Start a new parallel pool  
    print("Starting parallel pool with {0} threads".format(nThread))
    pool = Pool(processes=nThread)
    jobs = []
    
    for it in range(nSimTrials):
        # Determine if thread is already done by any previous session
        eidParam.SNR = simParam[it][0]
        eidParam.procNoiseSigma = simParam[it][1]
        eidParam.pLSigmaAmp = simParam[it][2]
        eidParam.pLSigmaAmpBayesian = simParam[it][2]
        eidParam.pLSigmaAmpEID = simParam[it][2]
        eidParam.Sigma = simParam[it][3]
        eidParam.objAmp = simParam[it][4]
        ergParam.dt = simParam[it][5]
        eidParam.UpdateDeltaT(simParam[it][5])
        if eidParam.simType == 'IF':
            eidParam.maxIter = round(eidParam.maxT / ergParam.dt)
        ergParam.wControl = simParam[it][6]
        eidParam.randSeed = simParam[it][7]
        ergParam.time = linspace(0.0, ergParam.timeHorizon, ergParam.tRes)
        ergParam.eidTime = linspace(0.0, ergParam.timeHorizon, eidParam.res)
        
        eidParam.filename = filename.replace('SNR', 'SNR-'+str(eidParam.SNR))
        eidParam.filename = eidParam.filename.replace('RandSeed', 'RandSeed-'+str(eidParam.randSeed))

        # Wait for all the thread to finish if the pool is full
        if len(jobs) == nThread:
            print("All thread used, waiting until finish... ...")
            # Wait until all the active thread to finish
            for job in jobs:
                job.join()
            
            print("Simulation Batch done!")
            # Clear the pool
            jobs = []
        
        # Start a new job thread
        #eidParam.maxIter = 3
        p = pool.Process(target=EIDSim, args=(ergParam, eidParam, False,False,False,True))
        p.start()
        jobs.append(p)
        
        print("Starting New Simulation Thread ({0:d}/{1:d}), progress {2:.2f}% ({3:d}/{4:d})".format(
            len(jobs), nThread,
            100.0*(it+1)/nSimTrials, it+1, nSimTrials))
    
    print("Waiting for remaining jobs to finish... ...")
    # Wait until all the active thread to finish
    for job in jobs:
        job.join()
        
        
#SimulationMain('./FigParameters/fig1/fig1-ErgodicHarvest-SNR.json', 1)