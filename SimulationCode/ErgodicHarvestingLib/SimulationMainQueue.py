from itertools import product
from multiprocessing import Pool, JoinableQueue, get_context
from queue import Empty
from os import makedirs, getpid, remove
from os.path import exists
from numpy import linspace
import numpy as np
from scipy.io import loadmat
import time
import tempfile
sim_jobs_path = tempfile.mkdtemp()
import pickle as pkl

from ErgodicHarvestingLib.Simulation import EIDSim
from ErgodicHarvestingLib.EIH_API_Sim_Entropy import EIH_Sim
from ErgodicHarvestingLib.ParameterIO import loadParams
from ErgodicHarvestingLib.Targets import RealTarget, Distractor
from ErgodicHarvestingLib.utils import print_color


# Load and normalize moth tracking data for EIH simulations
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
    sScale = np.max(s) * gain ** -1
    s /= sScale
    t /= sScale
    # Offset
    s += mid
    t += mid
    return s, t


def loadMothData(target="M300lux", trialID=0, nrmMid=0.5, nrmGain=0.1):
    traj = flatten(
        loadmat(
            "../Production-Figure-Code/PublishedData/animal_behavior_data/Moth/MothData.mat"
        )[f"trial_{target}"]
    )[trialID, :, :]
    s, t = normalize(traj[:, 0], traj[:, 1], nrmMid, nrmGain)
    return [s, t]


def QueueWorker(mp_queue):
    time.sleep(2.5)
    while True:
        try:
            job_path = mp_queue.get(block=True, timeout=30.0)
            remaining_jobs = mp_queue.qsize()
            with open(job_path, 'rb') as fp:
                args = pkl.load(fp)
            if args is None:
                continue
            remove(job_path)
            mp_queue.task_done()
            if isinstance(args, list):
                # wiggle attenuation sim
                print_color(
                    f"[WorkerNode-{getpid()}] Received new job {args[3]}, remaining jobs {remaining_jobs}",
                    color="yellow",
                )
                EIH_Sim(*args)
            else:
                # other sims
                print_color(
                    f"[WorkerNode-{getpid()}] Received new job {args[1].filename}, remaining jobs {remaining_jobs}",
                    color="yellow",
                )
                EIDSim(*args)
        except Empty:
            print_color(
                f"[WorkerNode-{getpid()}] no more work to be done, existing",
                color="yellow",
            )
            return
        except Exception:
            raise


def SimulationMainQueue(dataFiles, nThread=1):
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
                    paramList[-1]["SNR"],
                    paramList[-1]["procNoiseSigma"],
                    paramList[-1]["pLSigmaAmp"],
                    paramList[-1]["Sigma"],
                    paramList[-1]["objAmp"],
                    paramList[-1]["dt"],
                    paramList[-1]["wControl"],
                    paramList[-1]["randSeed"],
                )
            )
        )
        nSimJobsList.append(len(simParamList[-1]))

    nTrials = len(nSimJobsList)
    # additional wiggle attenuation sims
    with open("./SimParameters/SimJobList.txt", "r") as fp:
        attenuation_sim_trials = fp.readlines()
        attenuation_sim_trials.sort()
    nAttenuationSimTrials = len(attenuation_sim_trials)
    nRegularSimJobs = sum(nSimJobsList)
    nTotalJobs = nRegularSimJobs + nAttenuationSimTrials
    # Limit pool size when job size is smaller than total available threads
    if nTotalJobs < nThread:
        nThread = nTotalJobs
    # Start a new parallel pool
    print("Starting parallel pool with {0} threads".format(nThread))
    ctx = get_context("fork")
    pool = Pool(processes=nThread)
    work_queue = JoinableQueue(maxsize=nTotalJobs)
    jobs = []
    remaining_jobs = nTotalJobs
    job_id = 1

    for it in range(nAttenuationSimTrials):
        job_path = f'{sim_jobs_path}/eih_sim_job_{job_id}.pkl'
        with open(job_path, 'wb') as fp:
            pkl.dump(attenuation_sim_trials[it].split(), fp, pkl.HIGHEST_PROTOCOL)
        # Fill in work queue
        work_queue.put(job_path, block=True, timeout=None)
        job_id += 1
        remaining_jobs -= 1
        print_color(
            f"[MasterNode-{getpid()}]: Adding new job {attenuation_sim_trials[it].split()[3]}",
            color="green",
        )

    # Kick off worker threads
    for _ in range(nThread):
        # Start a new job thread
        try:
            p = pool.Process(target=QueueWorker, args=(work_queue,))
        except Exception:
            if ctx is not None:
                # Fallback to use context
                p = ctx.Process(target=QueueWorker, args=(work_queue,))
        p.start()
        jobs.append(p)

    for trial in range(nTrials):
        # Parse parameters
        param = paramList[trial]
        simParam = simParamList[trial]
        nJobs = nSimJobsList[trial]
        eidParam = param["eidParam"]
        ergParam = param["ergParam"]
        filename = param["filename"]
        # Check if saveDir exists, create the folder if not
        if not exists(eidParam.saveDir):
            print(f"Save folder {eidParam.saveDir} does not exist, creating...")
            makedirs(eidParam.saveDir, exist_ok=True)
        ergParam.time = None
        ergParam.eidTime = None
        for it in range(nJobs):
            eidParam.SNR = simParam[it][0]
            eidParam.procNoiseSigma = simParam[it][1]
            eidParam.pLSigmaAmp = simParam[it][2]
            eidParam.pLSigmaAmpBayesian = simParam[it][2]
            eidParam.pLSigmaAmpEID = simParam[it][2]
            eidParam.Sigma = simParam[it][3]
            eidParam.objAmp = simParam[it][4]
            ergParam.dt = simParam[it][5]
            eidParam.UpdateDeltaT(simParam[it][5])
            if eidParam.simType == "IF":
                eidParam.maxIter = round(eidParam.maxT / ergParam.dt)
            ergParam.wControl = simParam[it][6]
            eidParam.randSeed = simParam[it][7]
            eidParam.filename = (
                filename.replace("SNR", "SNR-" + str(eidParam.SNR))
                .replace("wC", "wC-" + str(ergParam.wControl))
                .replace("RandSeed", "RandSeed-" + str(eidParam.randSeed))
            )
            # Do the extra initialization here to speed up.
            if isinstance(eidParam.rawTraj, str) and "moth" in eidParam.rawTraj:
                eidParam.rawTraj = np.array(
                    loadMothData(
                        target="M300lux", trialID=0, nrmGain=eidParam.objAmp
                    )[1]
                )
            ergParam.time = linspace(0.0, ergParam.timeHorizon, ergParam.tRes)
            ergParam.eidTime = linspace(0.0, ergParam.timeHorizon, eidParam.res)
            # initialize multiple targets tracking
            if eidParam.multiTargetTracking == "dual":
                eidParam.multiTargetTracking = True
                eidParam.otherTargets = [RealTarget(eidParam.multiTargetInitialPos)]
            elif eidParam.multiTargetTracking == "distractor":
                eidParam.multiTargetTracking = True
                eidParam.otherTargets = [Distractor(eidParam.multiTargetInitialPos)]

            # Fill in work queue if it's not full
            job_path = f'{sim_jobs_path}/eih_sim_job_{job_id}.pkl'
            with open(job_path, 'wb') as fp:
                pkl.dump((ergParam, eidParam, False), fp, pkl.HIGHEST_PROTOCOL)
            work_queue.put(job_path, block=True, timeout=None)
            job_id += 1
            remaining_jobs -= 1
            print_color(
                f"[MasterNode-{getpid()}]: Adding new job {eidParam.filename}",
                color="green",
            )

    # Wait until all the active thread to finish
    work_queue.join()
    for job in jobs:
        job.join()
