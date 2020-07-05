from os import getpid
from multiprocessing import Pool, Queue, cpu_count, get_context
from queue import Empty
import time

from ErgodicHarvestingLib.EIH_API_Sim_Entropy import EIH_Sim
from ErgodicHarvestingLib.utils import print_color


def QueueWorker(mp_queue):
    while True:
        try:
            args, remaining_jobs = mp_queue.get(block=True, timeout=5.0)
            print_color(
                f"[WorkerNode-{getpid()}] Submitting new job {args[3]}", color="cyan"
            )
            EIH_Sim(*args)
        except Empty:
            print_color(
                f"[WorkerNode-{getpid()}] no more work to be done, existing",
                color="yellow",
            )
            return
        except Exception:
            raise


def EID_Attenuation_Sim(simDataFile, nThread):
    # Read simulation jobs
    f = open(simDataFile, "r")
    trials = f.readlines()
    nSimTrials = len(trials)
    # Submit simulations
    if nThread > cpu_count():
        nThread = cpu_count()
    if nSimTrials < nThread:
        nThread = nSimTrials
    # Start a new parallel pool
    print("Starting parallel pool with {0} threads".format(nThread))
    ctx = get_context("fork")
    pool = Pool(processes=nThread)
    max_queue_size = min(2 * nThread, nSimTrials)
    work_queue = Queue(maxsize=max_queue_size)
    jobs = []
    remaining_jobs = nSimTrials
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

    for it in range(nSimTrials):
        # Fill in work queue
        work_queue.put((trials[it].split(), remaining_jobs), block=True, timeout=None)
        remaining_jobs -= 1
        print_color(
            f"[MasterNode-{getpid()}]: Adding new job {trials[it].split()[3]}, "
            f"remaining jobs {remaining_jobs}",
            color="green",
        )
        # Unfortunately we need to wait briefly before adding new data into the queue.
        # This is because it takes some time for the object to get properly ingested.
        time.sleep(0.1)

    # Wait until all the active thread to finish
    for job in jobs:
        job.join()
