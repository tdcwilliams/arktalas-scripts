#! /usr/bin/env python
import os
import glob
import datetime as dt
import numpy as np
import pyproj
from shapely.geometry import Polygon

from sea_ice_drift import get_n


_DATA_DIR = os.getenv('RS2_dir')
_DAYS_IN_SEC = 24*3600
_THRESH = 3.1 # days: 3 -> 249 pairs (max=2.93, 2d 22h 19min); 3.1 -> 271 pairs (max=3.08, 3d 1h 55min)

def get_time(f):
    i = f.index('_HH')
    datestr = f[i-15:i]
    return dt.datetime.strptime(datestr, '%Y%m%d_%H%M%S')

def get_proj():
    if 0:
        # nextsim proj
        srs = '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-45 +x_0=0 +y_0=0 +a=6378273 +b=6356889.44891059 +units=m +no_defs' #nextsim proj
        return pyproj.Proj(srs)
    # Beaufort sea - good for breakup incident
    lon, lat = np.array([
        (-109.4, 74.7), # TR
        (-118.8, 67.5), # BR
        (-161.3, 70.5), # BL: Barrow
        (-153.9, 75.8), # TL
        ]).T
    lonc, latc = lon.mean() -2, lat.mean()
    srs = f'+proj=stere +lat_0={latc} +lat_ts={latc} +lon_0={lonc} +x_0=0 +y_0=0 +a=6378273 +b=6356889.44891059 +units=m +no_defs'
    return pyproj.Proj(srs)

def get_border(f, proj):
    """
    get border as a shapely polygon

    Parameters:
    -----------
    filename : str
    proj : pyproj.Proj

    Returns:
    --------
    p : shapely.geometry.Polygon
    """
    n = get_n(f, bandName='sigma0_HH', remove_spatial_mean=True)
    xb, yb = proj(*n.get_border())
    return Polygon([(x, y) for x, y in zip(xb, yb)])

def get_overlap(f1, f2, proj):
    """
    Parameters:
    -----------
    filename1 : str
    filename2 : str
    proj : pyproj.Proj

    Returns:
    --------
    overlap : float
        fractional overlap
    """
    p1 = get_border(f1, proj)
    p2 = get_border(f2, proj)
    intersection = p1.intersection(p2)
    overlaps = [intersection.area/p.area for p in [p1, p2]]
    return np.max(overlaps)

def filtered_pair_info(f1, f2, proj):
    """
    get pairs that are close enough in time
    and with > 30% overlap

    Parameters:
    -----------
    filename1 : str
    filename2 : str
    proj : pyproj.Proj

    Returns:
    --------
    info: tuple (if images are close enough in time and have enough overlap)
          or None
        if not None it is (basename1, basename2, interval, overlap):
            basename1 : str
                basename of 1st file
            basename2 : str
                basename of 2nd file
            interval : float
                time interval between the images
            overlap : float
                fractional overlap
    """
    dto1 = get_time(f1)
    dto2 = get_time(f2)
    diff = (dto2-dto1).total_seconds()/_DAYS_IN_SEC
    if diff > _THRESH:
        return
    overlap = get_overlap(f1, f2, proj)
    if overlap < .3:
        return
    info = (os.path.basename(f1), os.path.basename(f2), diff, overlap)
    print(info)
    return info

def save_pairs(pair_info, outfile):
    """
    Parameters:
    -----------
    pair_info : list(tuple)
        each tuple is (basename1, basename2, interval, overlap):
            basename1 : str
                basename of 1st file
            basename2 : str
                basename of 2nd file
            interval : float
                time interval between the images
            overlap : float
                fractional overlap
    outfile : str
        name of output file
    """
    print(f'Saving {len(pair_info)} pairs to {outfile}')
    maxdiff = 0
    with open(outfile, 'w') as f:
        f.write('File1,File2,"Time interval, days",Overlap\n')
        for f1, f2, diff, overlap in pair_info:
            f.write(f'{f1},{f2},{diff},{overlap}\n')
            maxdiff = np.max([maxdiff, diff])
    print(f'Max interval is {maxdiff} days.')

filelist = sorted(glob.glob(os.path.join(_DATA_DIR, 'RS2*')), key=get_time)
proj = get_proj() #pyproj to calc area of overlap
pair_info = []
os.makedirs('out', exist_ok=True)
outfile = 'out/RS2_pairs.csv'
for i, f1 in enumerate(filelist[:-1]):
    others = filelist[i+1:]
    for f2 in others:
        info = filtered_pair_info(f1, f2, proj)
        if info is not None:
            pair_info += [info]
            save_pairs(pair_info, outfile)
