#! /usr/bin/env python
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset
import datetime as dt
from matplotlib import pyplot as plt
import subprocess
import argparse

from pynextsim.gridding import Grid
from pynextsim.gmshlib import GmshMesh
from pynextsim.projection_info import ProjectionInfo
from pynextsim.irregular_grid_interpolator import IrregularGridInterpolator

_PLOT_INFO = dict(
        rls = ('Upward longwave flux, W/m$^2$', 'rls_%Y%m%dT%H%M%SZ.png', 'viridis', [0,150]),
        sit = ('Thickness, m', 'sit_%Y%m%dT%H%M%SZ.png', 'viridis', [0,3]),
        sic = ('Concentration', 'sic_%Y%m%dT%H%M%SZ.png', cmocean.cm.ice, [0,1]),
        wspeed = ('Wind speed, m/s', 'wspeed_%Y%m%dT%H%M%SZ.png', 'viridis', [0,10]),
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

def parse_args():
    parser = argparse.ArgumentParser("script to plot moorings file with zoom in on the Beaufort
            Sea")
    parser.add_argument('outdir', type=str,
            help='where to save results')
    parser.add_argument('moorings_file', type=str,
            help='path to moorings file')
    parser.add_argument('mesh_file', type=str,
            help='path to nextsim mesh file')
    parser.add_argument('-m', '--make-mp4', action='store_true',
            help='make mp4 movies from figures')
    return parser.parse_args()

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
    os.makedirs(os.path.dirname(figname), exist_ok=True)
    figname = os.path.join(outdir, vname, dto.strftime(figname))
    print(f'Saving {figname}')
    fig, ax = dst_grid.plot(array, clim=clim, cmap=cmap, clabel=clabel)#, land_zorder=1)
    ax.set_title(dto.strftime(f'%Y-%m-%d %H:%M (12h average)'))
    #ax.coastlines()
    gmo.boundary.plot(ax=ax)
    fig.savefig(figname)
    plt.close()

def make_plots(ds, dst_grid, gmfil, outdir, make_mp4=False):
    src_grid = Gridx(
            ds.variables['longitude'][:],
            ds.variables['latitude'][:],
            latlon=True
            )
    gmo = GmshMesh(gmfil, projection=dst_grid.projection)
    igi = src_grid.get_irregular_grid_interpolator(*dst_grid.lonlat, latlon=True)
    print(src_grid.shape, igi.src_shape)
    for i, time in enumerate(ds.variables['time'][:]):
        dto = dt.datetime(1900,1,1) + dt.timedelta(time)
        for vname in _PLOT_INFO:
            arr = ds.variables[vname][i].filled(np.nan)
            print(arr.shape)
            arr = igi.interp_field(arr)
            print(dst_grid.shape, arr.shape)
            plot(dst_grid, gmo, arr, vname, dto, outdir)

    if make_mp4:
        for vname in _PLOT_INFO:
            cmd = ['ffmpeg', '-framerate', '6',
                    '-pattern_type', 'glob', '-i', f'\'{outdir}/{vname}/{vname}*.png\'',
                    '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
                    os.path.join(outdir, vname, f'breakup_2013_wrf10km_ns10km_{vname}.mp4')]
            print(' '.join(cmd))
            subprocess.run(cmd)

args = parse_args()
dst_grid = get_target_grid()
with Dataset(args.moorings_file, 'r') as ds:
    make_plots(ds, dst_grid, args.mesh_file, args.outdir, make_mp4=args.make_mp4)
