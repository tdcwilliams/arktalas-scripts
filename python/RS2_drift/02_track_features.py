#! /usr/bin/env python
import os
import glob
import datetime as dt
import numpy as np
import pandas as pd

from sea_ice_drift import get_n
from sea_ice_drift.ftlib import feature_tracking

_DATA_DIR = os.getenv('RS2_dir')
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

df = pd.read_csv('out/RS2_pairs.csv')
df1 = pd.DataFrame(columns=df.columns)
for i in df.index:
    df_ = df.loc[df.index==i]
    f1 = os.path.join(os.getenv('RS2_dir'), df_['File1'].values[0])
    f2 = os.path.join(os.getenv('RS2_dir'), df_['File2'].values[0])
    tdiff = dt.timedelta(df_['Time interval, days'].values[0])
    print(f1, f2, sep='\n')
    print(tdiff)

    # create Nansat objects with one band only. 
    n1 = get_n(f1, bandName='sigma0_HH', remove_spatial_mean=True)
    n2 = get_n(f2, bandName='sigma0_HH', remove_spatial_mean=True)

    # Run Feature Tracking
    # get start/end coordinates in the image coordinate system (colums/rows)  
    c1, r1, c2, r2 = feature_tracking(n1, n2, nFeatures=20000, ratio_test=0.6)
    if len(c1) > 0:
        print('No tracked features')
        continue
    df1 = df1.append(df_, ignore_index=True)
    print(f'Number of pairs with features tracked = {len(df1)}')

os.makedirs('out', exist_ok=True)
outfile = 'out/RS2_pairs_ft.csv'
print(f'Saving {len(df1)} pairs to {outfile}')
df1.to_csv(outfile)
hi


# get the pairs that are close enough in time
filelist = sorted(glob.glob(os.path.join(_DATA_DIR, 'RS2*')), key=get_time)
pairs = []
for i, f1 in enumerate(filelist[:-1]):
    others = filelist[i+1:]
    for f2 in others:
        pairs += filtered_pair(f1, f2)

# save the pairs to file
os.makedirs('out', exist_ok=True)
outfile = 'out/RS2_pairs.csv'
print(f'Saving {len(pairs)} pairs to {outfile}')
maxdiff = 0
with open(outfile, 'w') as f:
    f.write('File1,File2,"Time interval, days"\n')
    for f1, f2, diff in pairs:
        f.write(f'{f1},{f2},{diff}\n')
        maxdiff = np.max([maxdiff, diff])
print(f'Max interval is {maxdiff} days.')
