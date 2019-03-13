import sys
from itertools import product
from multiprocessing import Pool, cpu_count
from scipy.io import loadmat, savemat
from numpy import linspace, arange, concatenate
from os import makedirs, scandir
from os.path import exists
from shutil import copy2
import numpy as np
# Suppress all warnings
np.seterr(all='ignore')
import warnings
warnings.filterwarnings("ignore")
# UTC Timestamp
import time
from datetime import datetime
# Ergodic Harvesting Code
from ErgodicHarvestingLib.SimParameters import ErgodicParameters, EIDParameters
from ErgodicHarvestingLib.Simulation import EIDSim
from ErgodicHarvestingLib.EIH_Attenuation_Sim_Main import EID_Attenuation_Sim
from ErgodicHarvestingLib.SimulationMain import SimulationMain

# get number of parallel threads, EAFP way
try:
    print(sys.argv[1])
    nThread = int(sys.argv[1])
except:
    nThread = cpu_count() # default

print(f'using {nThread} threads')

timeStampStart0 = datetime.fromtimestamp(time.time())
params = []
timeStampStart = datetime.fromtimestamp(time.time())

# load parameter files
paramPath = './SimParameters/'
if len(sys.argv) > 3:
    params.append(paramPath + sys.argv)
else:
    for f in scandir(paramPath):
        if f.is_file() and '.json' in f.name:
            params.append(paramPath + f.name)

nSimJobs = len(params)
print(f'Submitting {nSimJobs} total simulation jobs...')
print('---------------------------------------------------------------')
SimulationMain(params, nThread=nThread)

# wiggle attenuation simulations 
EID_Attenuation_Sim('./SimParameters/SimJobList.txt', nThread)

timeStampEnd = datetime.fromtimestamp(time.time())
timeString = timeStampEnd.strftime("%b-%d-%Y %T")
durationMinutes = (timeStampEnd - timeStampStart0).total_seconds() / 60.0
print("All done! EOF at {0}, total time taken for all simulation(s) {1:.2f} minutes".format(
    timeString, durationMinutes))