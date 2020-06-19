from ErgodicHarvestingLib.SimulationMainQueue import SimulationMainQueue
from datetime import datetime
import time
import warnings
import sys
from multiprocessing import cpu_count
from os import scandir
import numpy as np

# Suppress all warnings
np.seterr(all="ignore")
warnings.filterwarnings("ignore")

if __name__ == "__main__":
    # get number of parallel threads, EAFP way
    print(sys.argv)
    try:
        nThread = int(sys.argv[1])
    except BaseException:
        nThread = cpu_count()  # default

    print(f"using {nThread} threads")
    timeStampStart0 = datetime.fromtimestamp(time.time())
    timeStampStart = datetime.fromtimestamp(time.time())
    params = []

    # load parameter files
    paramPath = "./SimParameters/"
    for f in scandir(paramPath):
        if f.is_file() and ".json" in f.name:
            params.append(paramPath + f.name)
    # sort the file list so we have deterministic ordering
    params.sort()
    nSimJobs = len(params)
    print(f"Submitting {nSimJobs} total simulation jobs...")
    print("---------------------------------------------------------------")
    SimulationMainQueue(params, nThread=nThread)

    timeStampEnd = datetime.fromtimestamp(time.time())
    timeString = timeStampEnd.strftime("%b-%d-%Y %T")
    durationMinutes = (timeStampEnd - timeStampStart0).total_seconds() / 60.0
    print(
        "All done! EOF at {0}, total time taken for all simulation(s) {1:.2f} minutes".format(
            timeString, durationMinutes
        )
    )
