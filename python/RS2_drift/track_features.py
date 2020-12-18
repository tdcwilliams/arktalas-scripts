#! /usr/bin/env python
import os
import glob
import datetime as dt
import numpy as np
import pandas as pd

from sea_ice_drift import get_n
from sea_ice_drift.ftlib import feature_tracking

df = pd.read_csv('out/RS2_pairs.csv')
df1 = pd.DataFrame(columns=df.columns)
os.makedirs('out', exist_ok=True)
outfile = 'out/RS2_pairs_ft.csv'
num_ft = []
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
    if len(c1) == 0:
        print('No tracked features')
        continue
    num_ft += [len(c1)]
    df1 = df1.append(df_, ignore_index=True)
    print(f'Saving {len(df1)} pairs to {outfile}')
    df1.assign(NumFT=num_ft).to_csv(outfile)
