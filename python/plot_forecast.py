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

class PlotAtmosphericForcing:
    def __init__(self):
        #should be inputs
        self.data_dir = '/cluster/projects/nn9624k/wrf_init/march2017/'
        self.template = Template(os.path.join(
            self.data_dir, 'od.ans.201302-201303.sfc.${varname}.nc'))
        self.ij_range = [0, 160, 0, 1440]
        self.plot_res = 10
        self.avg_period_hours = 12
        self.mesh_file = os.path.join(os.getenv('NEXTSIM_MESH_DIR'), 'wrf_arctic_10km.msh')
        self.outdir = 'figs/ecmwf_fc'
        #self.varnames = ['t2m', 'sp']
        self.varnames = []
        self.make_wind_plots = True
        self.set_grid_igi()

    @property
    def src_lonlat(self):
        f = self.template.safe_substitute(dict(varname='t2m'))
        return mnu.nc_getinfo(f).get_lonlat(ij_range=self.ij_range)

    def set_grid_igi(self):
        # get target grid
        print('Getting target grid')
        self.gmsh = GmshMesh(self.mesh_file)
        self.grid = self.gmsh.boundary.get_grid(resolution=self.plot_res*1000)

        # get interpolator
        print('Getting interpolator')
        self.igi = self.grid.get_interpolator(self.src_lonlat, interp_from=False, latlon=True)

    @staticmethod
    def floor_time(dto, avg_period_hours=12):
        h = avg_period_hours*(dto.hour//avg_period_hours)
        return dt.datetime(dto.year, dto.month, dto.day, h)

    def average_data(self, nci, varname, dto0, dto1):
        dtimes = np.array(nci.datetimes)
        dtimes = dtimes[(dtimes>=dto0)*(dtimes<=dto1)]
        n = len(dtimes)
        if n == 0:
            return
        indices = [nci.datetimes.index(dto) for dto in dtimes]
        data = 0.
        for i in indices:
            data += nci.get_var(varname,
                    time_index=i, ij_range=self.ij_range).values/n
        return data.filled(np.nan)

    def get_averaged_data(self, nci, varname):
        delt = dt.timedelta(hours=self.avg_period_hours)
        dto0 = self.floor_time(nci.datetimes[0] , avg_period_hours=self.avg_period_hours)
        dto2 = self.floor_time(nci.datetimes[-1], avg_period_hours=self.avg_period_hours)
        while dto0 <= dto2:
            dto1 = dto0 + delt
            yield dto0, dto1, self.average_data(nci, varname, dto0, dto1)
            dto0 = dto1

    def plot_scalar(self, varname, dints=None):
        fig, ax = self.grid.plot(data, add_landmask=False, #clim=[-10,10],
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
        figdir = os.path.join(self.outdir, varname)
        figname = os.path.join(figdir, f'ecmwf_fc_{varname}{datestr2}.png')
        print(f'Saving {figname}')
        os.makedirs(figdir, exist_ok=True)
        fig.savefig(figname)
        plt.close()

    def plot_wind(self, uv, dints, step=10):
        spd = np.hypot(*uv)
        fig, ax = self.grid.plot(spd, add_landmask=False,
                cmap='viridis', clim=[0,10], clabel='Wind speed, m/s')
        ax.coastlines(resolution='50m')

        # wind speed contours
        #cs = ax.contour(*self.grid.xy, spd, levels=[10,15,20], colors='b')
        #ax.clabel(cs, inline=True, fontsize=10)

        # wind vectors
        x = self.grid.xy[0][::step,::step] 
        y = self.grid.xy[1][::step,::step]
        u = uv[0][::step,::step]
        v = uv[1][::step,::step]
        u, v = nsl.rotate_lonlat2xy(self.grid.projection, x, y, u, v)
        ax.quiver(x, y, u, v, units='xy', angles='xy', color='r')

        # tidy up
        datestr = '\n' + ' - '.join([
            dto.strftime('%Y-%m-%d %H:%M') for dto in dints])
        ax.set_title(f'ECMWF FC{datestr}')
        datestr = '_' + '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
        figdir = os.path.join(self.outdir, 'wind')
        figname = os.path.join(figdir, f'ecmwf_fc_wind{datestr}.png')

        #fig.show()
        print(f'Saving {figname}')
        os.makedirs(figdir, exist_ok=True)
        fig.savefig(figname)
        plt.close()

    def test_plot(self):
        ''' test interpolation by plotting lon, lat '''
        #print(lonlat)
        for varname, v in zip(['lon', 'lat'], self.src_lonlat):
            data = self.igi.interp_field(v)
            self.plot_scalar(data, varname)

    def run(self):
        #self.test_plot()

        #for varname in self.varnames:
        #    #scalar fields
        #    print(f'Making plots for {varname}')
        #    f = self.template.safe_substitute(dict(varname=varname))
        #    for dto0, dto1, v in self.get_averaged_data(
        #            mnu.nc_getinfo(f), varname)
        #        if v is None:
        #            continue
        #        data = self.igi.interp_field(v)
        #        if varname in ['t2m', 'd2m']:
        #            data -= 273.15 #kelvin to deg C
        #        self.plot_scalar(data, varname, dints=(dto0, dto1))

        if self.make_wind_plots:
            #wind
            print(f'Making plots for wind speed')
            pairs = []
            for varname in ['u10', 'v10']:
                f = self.template.safe_substitute(dict(varname=varname))
                pairs += [(mnu.nc_getinfo(f), varname)]
            for (dto0, dto1, u10), (_, _, v10) in zip(
                    self.get_averaged_data(*pairs[0]),
                    self.get_averaged_data(*pairs[1]),
                    ):
                if u10 is None or v10 is None:
                    continue
                data = [self.igi.interp_field(v) for v in [u10, v10]]
                self.plot_wind(data, (dto0, dto1))

if __name__ == "__main__":
    obj = PlotAtmosphericForcing()
    obj.run()
