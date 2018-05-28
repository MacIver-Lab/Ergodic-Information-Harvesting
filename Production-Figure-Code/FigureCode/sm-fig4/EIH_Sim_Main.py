from EIH_API_Sim_Entropy import EIH_Sim
from scipy.io import loadmat
from sys import argv
from numpy import linspace
from os import getcwd
from itertools import product
from multiprocessing import Pool, cpu_count
# UTC Timestamp
from time import time

if __name__ == '__main__': 
    # Read simulation jobs
    if len(argv) > 1:
        # Read input arguments
        simDataFile = argv[1]
        nThread = int(argv[2])
        f = open(simDataFile, 'r')
        trials = f.readlines()
        nSimTrials = len(trials)

        # Submit simulations
        if nThread > cpu_count():
            nThread = cpu_count()
        if nSimTrials < nThread:
            nThread = nSimTrials
        # Start a new parallel pool  
        print("Starting parallel pool with {0} threads".format(nThread))
        pool = Pool(processes=nThread)
        jobs = []
        for it in range(nSimTrials):
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
            p = pool.Process(target=EIH_Sim, args=trials[it].split())
            p.start()
            jobs.append(p)
            
            print("Starting New Simulation Thread ({0:d}/{1:d}), progress {2:.2f}% ({3:d}/{4:d})".format(
                len(jobs), nThread,
                100.0*it/nSimTrials, it+1, nSimTrials))
        print("Waiting for remaining jobs to finish... ...")
        # Wait until all the active thread to finish
        for job in jobs:
            job.join()
            
        timeStamp = "{0:.0f}".format(time())
        print("All done! EOF at time stamp {0}".format(timeStamp))
    else:
        print(getcwd())
        print('Invalid input parameters! len={0}, param={1}'.format(len(argv),argv))