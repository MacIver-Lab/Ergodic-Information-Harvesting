import numpy as np
import sys

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as fp:
        log = fp.read()
    logs = [line.split(' ')[-1][:-5] for line in log.split('\n') if 'Received' in line]
    unique_logs = np.unique(logs)
    print(f'Found {len(unique_logs)} total unique sims, expected 938, diff is {938-len(unique_logs)}')
    assert len(logs) == len(unique_logs)