#! /usr/bin/env python
import os
import glob
import datetime as dt
import numpy as np

_DATA_DIR = os.path.join(os.getenv('HOME'), 'docker_io/RS2_beaufort_2013/')
_DAYS_IN_SEC = 24*3600
_THRESH = 3.1 # days: 3 -> 249 pairs (max=2.93, 2d 22h 19min); 3.1 -> 271 pairs (max=3.08, 3d 1h 55min)

def get_time(f):
    i = f.index('_HH')
    datestr = f[i-15:i]
    return dt.datetime.strptime(datestr, '%Y%m%d_%H%M%S')

def filtered_pair(f1, f2):
    dto1 = get_time(f1)
    dto2 = get_time(f2)
    diff = (dto2-dto1).total_seconds()/_DAYS_IN_SEC
    if diff > _THRESH:
        return []
    return [(os.path.basename(f1), os.path.basename(f2), diff)]

# get the pairs that are close enough in time
filelist = sorted(glob.glob(os.path.join(_DATA_DIR, 'RS2*')), key=get_time)
pairs = []
for i, f1 in enumerate(filelist[:-1]):
    others = filelist[i+1:]
    for f2 in others:
        pairs += filtered_pair(f1, f2)

# save the pairs to file
outfile = 'RS2_pairs.csv'
print(f'Saving {len(pairs)} pairs to {outfile}')
maxdiff = 0
with open(outfile, 'w') as f:
    f.write('File1,File2\n')
    for f1, f2, diff in pairs:
        f.write(f'{f1},{f2},{diff}\n')
        maxdiff = np.max([maxdiff, diff])
print(f'Max interval is {maxdiff} days.')
