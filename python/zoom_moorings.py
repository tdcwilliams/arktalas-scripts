#! /usr/bin/env python
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset
import datetime as dt
from matplotlib import pyplot as plt

from pynextsim.gridding import Grid
from pynextsim.gmshlib import GmshMesh
from pynextsim.projection_info import ProjectionInfo
from pynextsim.irregular_grid_interpolator import IrregularGridInterpolator

_PLOT_INFO = dict(
        #rls = ('Upward longwave flux, W/m$^2$', 'rls_%Y%m%dT%H%M%SZ.png', 'viridis', [0,150]),
        sit = ('Thickness, m', 'sit_%Y%m%dT%H%M%SZ.png', 'viridis', [0,3]),
        #sic = ('Concentration', 'sic_%Y%m%dT%H%M%SZ.png', cmocean.cm.ice, [0,1]),
        #wspeed = ('Wind speed, m/s', 'wspeed_%Y%m%dT%H%M%SZ.png', 'viridis', [0,10]),
        )

class Gridx(Grid):
    ''' extend Grid class '''
    @property
    def coord_vectors_regular_xy(self):
        """
        For a regular grid, return the coordinate vectors y, x

        Returns:
        --------
        y: numpy.ndarray
        x: numpy.ndarray
        """
        assert(self.regular_xy)
        return grid.xy[1][:,0], grid.xy[0][0]

    def get_regular_grid_interpolator(self, field):
        """
        Parameters:
        -----------
        field : numpy.ndarray
            2d field to be interpolatated from the regular grid,
            shape ny,nx

        Returns:
        --------
        rgi : scipy.interpolate.RegularGridInterpolator
        """
        return RegularGridInterpolator(
                self.coord_vectors_regular_xy, field, **kwargs)

    def get_irregular_grid_interpolator(self, xout, yout, latlon=False):
        """
        Parameters:
        -----------
        xout : numpy.ndarray
        yout : numpy.ndarray
        """
        if latlon:
            x, y = self.projection.pyproj(xout, yout)
        else:
            x, y = xout, yout
        return IrregularGridInterpolator(*self.xy, x, y)

    def interpolate(self, data, xout, yout, latlon=False, igi=None):
        if not igi:
            if latlon:
                x, y = self.projection.pyproj(xout, yout)
            else:
                x, y = xout, yout
            igi = self.get_irregular_grid_interpolator(x, y)
        out = dict()
        for vname, arr in data.items():
            out[vname] = igi.interpolate(arr)
        return out, igi

def get_target_grid():
    res = 2000 # resolution in meters of target grid
    factors = dict(xmin=1.8, xmax=1, ymin=1.5, ymax=1.3) # for tuning extent of grid
    lon, lat = np.array([
        (-109.4, 74.7), # TR
        (-118.8, 67.5), # BR
        (-161.3, 70.5), # BL: Utqiagvik (trad. name of Barrow)
        (-153.9, 75.8), # TL
        ]).T
    lonc, latc = lon.mean() - 2, lat.mean()
    print(lonc, latc)
    proj = ProjectionInfo(lat_0=latc, lon_0=lonc, lat_ts=latc)
    x, y = proj.pyproj(lon, lat)
    xc, yc = np.mean(x), np.mean(y)
    dx = x.max() - xc
    dy = y.max() - yc
    grid_params = dict(
            xmin = xc - factors['xmin']*dx,
            xmax = xc + factors['xmax']*dx,
            ymin = yc - factors['ymin']*dy,
            ymax = yc + factors['ymax']*dy,
            )
    grid_params['nx'] = int((grid_params['xmax'] - grid_params['xmin'])/res)
    grid_params['ny'] = int((grid_params['ymax'] - grid_params['ymin'])/res)
    return Grid.init_from_grid_params(grid_params, projection=proj)

def plot(dst_grid, gmo, array, vname, dto, outdir):
    clabel, figname, cmap, clim = _PLOT_INFO[vname]
    figname = os.path.join(outdir, dto.strftime(figname))
    print(f'Saving {figname}')
    fig, ax = dst_grid.plot(array, clim=clim, cmap=cmap, clabel=clabel)#, land_zorder=1)
    ax.set_title(dto.strftime(f'%Y-%m-%d %H:%M (12h average)'))
    #ax.coastlines()
    gmo.boundary.plot(ax=ax)
    fig.savefig(figname)
    plt.close()

def make_plots(ds, dst_grid, gmfil, outdir):
    os.makedirs(outdir, exist_ok=True)
    src_grid = Gridx(
            ds.variables['longitude'][:],
            ds.variables['latitude'][:],
            latlon=True
            )
    gmo = GmshMesh(gmfil, projection=dst_grid.projection)
    igi = src_grid.get_irregular_grid_interpolator(*dst_grid.lonlat, latlon=True)
    print(src_grid.shape, igi.src_shape)
    #i = -1; time = float(ds.variables['time'][i])
    for i, time in enumerate(ds.variables['time'][:]):
        dto = dt.datetime(1900,1,1) + dt.timedelta(time)
        for vname in _PLOT_INFO:
            arr = ds.variables[vname][i].filled(np.nan)
            print(arr.shape)
            arr = igi.interp_field(arr)
            print(dst_grid.shape, arr.shape)
            plot(dst_grid, gmo, arr, vname, dto, outdir)
    
dst_grid = get_target_grid()
rootdir_5km = '/cluster/work/users/timill/nextsim-stand-alone/wrf_arctic_5km/breakup_2013'
rootdir_10km = '/cluster/work/users/timill/nextsim-stand-alone/wrf_arctic_10km/breakup_2013'
if 0:
    moorings = os.path.join(rootdir_10km, 'expt_00_wrf10km', 'outputs', 'Moorings.nc')
    outdir = 'figs/zoom_10km'
    gmfil = '/cluster/projects/nn2993k/sim/mesh/wrf_arctic_10km.msh'
else:
    moorings = os.path.join(rootdir_5km, 'expt_01_wrf_10km-C1.5', 'outputs', 'Moorings.nc')
    outdir = 'figs/zoom_5km'
    gmfil = '/cluster/projects/nn2993k/sim/mesh/wrf_arctic_5km.msh'
with Dataset(moorings, 'r') as ds:
    make_plots(ds, dst_grid, gmfil, outdir)
