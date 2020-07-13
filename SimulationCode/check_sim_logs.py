import numpy as np
import sys

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as fp:
        log = fp.read()
    logs = np.unique([line.split(' ')[-1][:-5] for line in log.split('\n') if 'Received' in line])
    print(f'Found {len(logs)} total unique sims, expected 938, diff is {938-len(logs)}')