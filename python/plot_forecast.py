#! /usr/bin/env python
import os
import numpy as np
from string import Template
import datetime as dt
from matplotlib import pyplot as plt

from pynextsim.gmshlib import GmshMesh
from pynextsim.gridding import Grid
from pynextsim.openers import OpenerVariable, Opener
import pynextsim.lib as nsl

import mod_netcdf_utils as mnu

_DATA_DIR = '/cluster/projects/nn9624k/wrf_init/march2017/'
_REFDATE = dt.datetime(1900, 1, 1)
#_IJ_RANGE = []

#should be inputs
varname = 't2m'
dto = dt.datetime(2013,2,22)
meshfile = os.path.join(os.getenv('NEXTSIM_MESH_DIR'), 'wrf_arctic_10km.msh')
plot_res = 10 #km

# get target grid
gmsh = GmshMesh(meshfile)
grid = gmsh.boundary.get_grid(resolution=plot_res*1000)

def floor_time(dto, avg_period_hours=12):
    h = avg_period_hours*(dto.hour//avg_period_hours)
    return dt.datetime(dto.year, dto.month, dto.day, h)

def average_data(nci, varname, dto0, dto1):
    dtimes = np.array(nci.datetimes)
    dtimes = dtimes[(dtimes>=dto0)*(dtimes<=dto1)]
    n = len(dtimes)
    if n == 0:
        return
    indices = [nci.datetimes.index(dto) for dto in dtimes]
    data = 0.
    for i in indices:
        data += nci.get_var(varname, time_index=i).values/n
    return data.filled(np.nan)

def get_averaged_data(nci, varname, avg_period_hours=12):
    delt = dt.timedelta(hours=avg_period_hours)
    dto0 = floor_time(nci.datetimes[0] , avg_period_hours=avg_period_hours)
    dto2 = floor_time(nci.datetimes[-1], avg_period_hours=avg_period_hours)
    while dto0 <= dto2:
        dto1 = dto0 + delt
        av_data = average_data(nci, varname, dto0, dto1)
        yield dto0, dto1, average_data(nci, varname, dto0, dto1)
        dto0 = dto1

def plot_wind(grid, uv, dints, outdir, step=10):
    spd = np.hypot(*uv)
    fig, ax = grid.plot(spd, add_landmask=False,
            cmap='viridis', clim=[0,10], clabel='Wind speed, m/s')
    ax.coastlines(resolution='50m')
    x = grid.xy[0][::step,::step] 
    y = grid.xy[1][::step,::step]
    u = uv[0][::step,::step]
    v = uv[1][::step,::step]
    u, v = nsl.rotate_lonlat2xy(grid.projection, x, y, u, v)
    ax.quiver(x, y, u, v, units='xy', angles='xy', color='r')

    # tidy up
    ttl = 'ECMWF FC\n' + ' - '.join([
        dto.strftime('%Y-%m-%d %H:%M') for dto in dints])
    figname = os.path.join(outdir, 'ecmwf_fc_%s.png' %(
        '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])))
    ax.set_title(ttl)

    #fig.show()
    print(f'Saving {figname}')
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(figname)
    plt.close()

if 0:
    #scalar field

    # get source file
    op = OpenerEra5(varname)
    f = op.find(dto) # filename for date if it exists
    nci = mnu.nc_getinfo(f)
    tind = nci.datetimes.index(dto)
    data = grid.get_netcdf_data(nci, vlist=[varname], time_index=tind)[varname]
    fig, ax = grid.plot(data, cmap='viridis', add_landmask=False)
    ax.coastlines(resolution='50m')
    fig.show()
else:
    #wind
    t = Template(os.path.join(_DATA_DIR, 'od.ans.201302-201303.sfc.${varname}.nc'))
    pairs = []
    av_per_h = 12
    step = 10
    outdir = 'figs/ecmwf_fc'
    for varname in ['u10', 'v10']:
        f = t.safe_substitute(dict(varname=varname))
        pairs += [(mnu.nc_getinfo(f), varname)]
    igi = grid.get_interpolator(
            pairs[0][0].get_lonlat(), interp_from=False, latlon=True)
    for u10, v10 in zip(
            get_averaged_data(*pairs[0], avg_period_hours=av_per_h),
            get_averaged_data(*pairs[1], avg_period_hours=av_per_h)
            ):
        print(u10[:2])
        print(u10[2].shape)
        print(v10[2].shape)
        data = []
        if u10[2] is None or v10[2] is None:
            continue
        for tup in [u10, v10]:
            data += [igi.interp_field(tup[2])]
        plot_wind(grid, data, u10[:2], outdir, step=step)
    hi

    for varname in ['u10', 'v10']:
        # get source file
        f = t.safe_substitute(dict(varname=varname))
        nci = mnu.nc_getinfo(f)
        tind = nci.datetimes.index(dto)
        data += [grid.get_netcdf_data(nci, vlist=[varname], time_index=tind)[varname]]

    spd = np.hypot(*data) 
    fig, ax = grid.plot(spd, cmap='viridis', add_landmask=False)
    ax.coastlines(resolution='50m')
    x = grid.xy[0][::10,::10] 
    y = grid.xy[1][::10,::10]
    u = data[0][::10,::10]
    v = data[1][::10,::10]
    u, v = nsl.rotate_lonlat2xy(grid.projection, x, y, u, v)
    ax.quiver(x, y, u, v, units='xy', angles='xy', color='r')
    ax.set_title(f'ECMWF FC\n{dto.strftime("%Y-%m-%d %H:%M")}')
    #fig.show()
    fig.savefig('test.png')
