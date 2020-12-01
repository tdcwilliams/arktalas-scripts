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
_IJ_RANGE = [0, 160, 0, 1440]

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
        data += nci.get_var(varname,
                time_index=i, ij_range=_IJ_RANGE).values/n
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
    datestr = '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
    figname = os.path.join(outdir, 'ecmwf_fc_wind_%s.png' %datestr)
    ax.set_title(ttl)

    #fig.show()
    print(f'Saving {figname}')
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(figname)
    plt.close()

def plot_scalar(grid, data, varname, outdir, dints=None):
    fig, ax = grid.plot(data, add_landmask=False, #clim=[-10,10],
            cmap='viridis', clabel='2-m air temperature, $^\circ$C')
    ax.coastlines(resolution='50m')
    print(data.min(), data.max())

    # tidy up
    datestr1 = ''
    datestr2 = ''
    if dints is not None:
        datestr1 = '\n' + ' - '.join([dto.strftime('%Y-%m-%d %H:M') for dto in dints])
        datestr2 = '_' + '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
    ax.set_title(f'ECMWF FC{datestr1}')

    #fig.show()
    figname = os.path.join(outdir, f'ecmwf_fc_{varname}{datestr2}.png')
    print(f'Saving {figname}')
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(figname)
    plt.close()

def test_plot(grid, igi, lonlat, outdir):
    ''' test interpolation by plotting lon, lat '''
    #print(lonlat)
    for varname, v in zip(['lon', 'lat'], lonlat):
        data = igi.interp_field(v)
        plot_scalar(grid, data, varname, outdir)

#should be inputs
dto = dt.datetime(2013,2,22)
meshfile = os.path.join(os.getenv('NEXTSIM_MESH_DIR'), 'wrf_arctic_10km.msh')
plot_res = 10 #km
av_per_h = 12
step = 10
outdir = 'figs/ecmwf_fc'
#varnames = ['t2m']
varnames = ['u10', 'v10']

# get target grid
print('Getting target grid')
gmsh = GmshMesh(meshfile)
grid = gmsh.boundary.get_grid(resolution=plot_res*1000)

# get interpolator
print('Getting interpolator')
t = Template(os.path.join(_DATA_DIR, 'od.ans.201302-201303.sfc.${varname}.nc'))
f = t.safe_substitute(dict(varname=varnames[0]))
lonlat = mnu.nc_getinfo(f).get_lonlat(ij_range=_IJ_RANGE)
igi = grid.get_interpolator(lonlat, interp_from=False, latlon=True)
#test_plot(grid, igi, lonlat, outdir)

if len(varnames) == 1:
    #scalar field
    varname = varnames[0]
    print(f'Making plots for {varname}')
    f = t.safe_substitute(dict(varname=varnames[0]))
    for dto0, dto1, v in get_averaged_data(
            mnu.nc_getinfo(f), varname, avg_period_hours=av_per_h):
        if v is None:
            continue
        data = igi.interp_field(v)
        if varname in ['t2m', 'd2m']:
            data -= 273.15 #kelvin to deg C
        os.makedirs(outdir, exist_ok=True)
        plot_scalar(grid, data, varname, outdir, dints=(dto0, dto1))
else:
    #wind
    print(f'Making plots for wind speed')
    pairs = []
    for varname in varnames:
        f = t.safe_substitute(dict(varname=varname))
        pairs += [(mnu.nc_getinfo(f), varname)]
    for (dto0, dto1, u10), (_, _, v10) in zip(
            get_averaged_data(*pairs[0], avg_period_hours=av_per_h),
            get_averaged_data(*pairs[1], avg_period_hours=av_per_h)
            ):
        dints = (dto0, dto1)
        if u10 is None or v10 is None:
            continue
        data = [igi.interp_field(v) for v in [u10, v10]]
        plot_wind(grid, data, dints, outdir, step=step)
