#! /usr/bin/env python
import os
from argparse import ArgumentParser
import numpy as np
import datetime as dt

from scipy.interpolate import RegularGridInterpolator

from pynextsim.projection_info import ProjectionInfo
import mod_netcdf_utils as mnu

NS_PROJ = ProjectionInfo()

def parse_args():
    p = ArgumentParser('Script to compare moorings with RS2-derived drift')
    p.add_argument('rs2_file', help='path to npz file with pattern matching')
    p.add_argument('moorings_file', help='path to moorings file')
    p.add_argument('outdir', help='Where to save results')
    return p.parse_args()

def read_rs2_file(rs2_file):
    print(f'Reading                 : {rs2_file}')
    # get times
    basename, _ = os.path.splitext(os.path.basename(rs2_file))
    dto1, dto2 = [
            dt.datetime.strptime(s, '%Y%m%dT%H%M%SZ')
            for s in basename.split('_')[1].split('-')]
    print(f'Start                   : {dto1}')
    print(f'End                     : {dto2}')
    print(f'Interval                : {dto2 - dto1}')

    # get spatial info from pattern-matching results
    pm_results = dict(np.load(rs2_file))
    gpi = np.isfinite(pm_results['upm_clean']*pm_results['vpm_clean'])
    xy1 = NS_PROJ.pyproj(
            pm_results['lon1pm'][gpi], pm_results['lat1pm'][gpi])
    xy2 = NS_PROJ.pyproj(
            pm_results['lon2pm'][gpi], pm_results['lat2pm'][gpi])

    print(f'Number of drift vectors : {np.sum(gpi)}')
    return dto1, dto2, pm_results, xy1, xy2

def match_files(nci, xy1, xy2):
    x, y = NS_PROJ.pyproj(*nci.get_lonlat())
    x = x[0] 
    y = y[:,0].flatten()
    xmin = np.min([xy1[0].min(), xy2[0].min()])
    xmax = np.max([xy1[0].max(), xy2[0].max()])
    ymin = np.min([xy1[1].min(), xy2[1].min()])
    ymax = np.max([xy1[1].max(), xy2[1].max()])
    def reduce_range(a, amin, amax):
        gpi = (a>=amin)*(a<=amax)
        a_ = a[gpi]
        amin_ = a_.min()
        amax_ = a_.max()
        lst = list(a)
        return a_, lst.index(amin_), lst.index(amax_)
    x, j0, j1 = reduce_range(x, xmin, xmax)
    y, i0, i1 = reduce_range(y, ymin, ymax)
    return x, y, [i0, i1+1, j0, j1+1]

def time_iterator(nci, dto1, dto2):

    # default time resolution
    dt_ref = nci.datetimes[1] - nci.datetimes[0]

    # 1st time step info
    dto1_, time_index1 = nci.nearestDate(dto1)
    if dto1_ - .5*dt_ref < dto1:
        time_index1 -= 1
        dto1_ = nci.datetimes[time_index1]
    dt1 = (dto1_ + .5*dt_ref - dto1).total_seconds()

    # last time step info
    dto2_, time_index2 = nci.nearestDate(dto2)
    if dto2_ + .5*dt_ref > dto2:
        time_index2 += 1
        dto2_ = nci.datetimes[time_index2]
    dt2 = (dto2 - dto2_ + .5*dt_ref).total_seconds()

    # Yield the 1st time index and integration time
    yield time_index1, dt1
    for time_index in range(time_index1+1, time_index2):
        # Yield the intermediary time indices and integration times
        yield time_index, dt_ref.total_seconds()
    # Yield the last time index and integration time
    yield time_index2, dt2

def get_interpolators(time_index, nci, x, y, ij_range):
    siu = nci.get_var(
            'siu', ij_range=ij_range, time_index=time_index).values
    siv = nci.get_var(
            'siv', ij_range=ij_range, time_index=time_index).values
    return [RegularGridInterpolator([y, x], a,
            bounds_error=False, fill_value=np.nan)
            for a in [siu, siv]]

def integrate_one_time_step(x0, y0, dt, *args):
    dst_p = np.array([y0, x0]).T
    rgi_uv = get_interpolators(*args)
    u0, v0 = np.array([rgi(dst_p) for rgi in rgi_uv])
    return x0 + u0*dt, y0 + v0*dt

def integrate_velocities(x0, y0, dto1, dto2,
        nci, xy_info):

    x = np.array(x0)
    y = np.array(y0)
    dt_tot = 0
    for i, dt_i in time_iterator(nci, dto1, dto2):
        x, y = integrate_one_time_step(
                x, y, dt_i, i, nci, *xy_info)
        dt_tot += dt_i
        print(i, np.any(np.isnan(x*y)))
    assert(dt_tot == (dto2-dto1).total_seconds())
    return (x, y), dt_tot

def compare(xy1, xy2, xy2_mod, dt_tot):
    fac = 24*3600*1e-3 #m/s to km/day
    dx_obs, dy_obs = [xy2[i] - xy1[i] for i in range(2)]
    dx_mod, dy_mod = [xy2_mod[i] - xy1[i] for i in range(2)]
    # bias in speed
    spd_obs = fac*np.hypot(dx_obs, dy_obs)/dt_tot
    spd_mod = fac*np.hypot(dx_obs, dy_obs)/dt_tot
    bias_speed = np.mean(spd_mod - spd_obs)
    print(f'Bias in speed = {bias_speed} km/day')
    # RMSE in speed
    rmse_speed = np.sqrt(np.mean((spd_mod - spd_obs)**2))
    print(f'RMSE in speed = {rmse_speed} km/day')
    # Vector RMSE
    vdiff = fac*np.hypot(dx_mod - dx_obs, dy_mod - dy_obs)/dt_tot
    vrmse = np.sqrt(np.mean(vdiff**2))
    print(f'VMRSE = {vrmse} km/day')

def run():
    args = parse_args()
    dto1, dto2, pm_results, xy1, xy2 = read_rs2_file(args.rs2_file) 
    nci = mnu.nc_getinfo(args.moorings_file)
    xy_info = match_files(nci, xy1, xy2)
    xy2_ns, dt_tot = integrate_velocities(
            *xy1, dto1, dto2, nci, xy_info)
    compare(xy1, xy2, xy2_ns, dt_tot)

if __name__ == '__main__':
    run()
