#! /usr/bin/env python
import os
import numpy as np
from string import Template
import datetime as dt
from matplotlib import pyplot as plt
import argparse

from pynextsim.gmshlib import GmshMesh
from pynextsim.gridding import Grid
from pynextsim.nextsim_config import NextsimConfig
from pynextsim.projection_info import ProjectionInfo
import pynextsim.lib as nsl

import mod_netcdf_utils as mnu

PROJECTIONS = dict(
        wrf= ProjectionInfo(proj = 'stere', ecc = 0, a = 6370e3,
            lat_0  = 90., lon_0  = 30., lat_ts = 81.),
        ec2_stere=ProjectionInfo(),
        ec2_arome=ProjectionInfo(),
        )

class PlotAtmForcing:
    """ Base class for plotting the different forcings """

    def __init__(self, cli=None):
        self.args = self.parse_args(cli)
        self.projection = PROJECTIONS.get(self.args.forcing, None)
        self.do_interp = self.projection is None
        self.igi = None
        if self.do_interp:
            # interp from lon, lat to nextsim stereographic projection
            self.projection = ProjectionInfo()
        self.parse_config_file()

    @staticmethod
    def parse_args(cli=None):
        parser = argparse.ArgumentParser("script to plot atmospheric forcing")
        parser.add_argument('config_file', type=str,
                help='path to config file with input settings')
        parser.add_argument('forcing', type=str,
                choices=['wrf', 'era5', 'cfsr', 'ec2', 'ec2_stere', 'ec2_arome'],
                help='type of atmospheric forcing')
        return parser.parse_args(cli)

    def parse_config_file(self):
        opts = NextsimConfig(self.args.config_file)['plot_atm_forcing']
        for att, f_cvt in [
                ('title', str),
                ('data_dir', str),
                ('ij_range', lambda x : [int(i) for i in x.split()]),
                ('plot_res', float),
                ('avg_period_hours', int),
                ('outdir', str),
                ('output_prefix', str),
                ('temp_names', lambda x : x.split()),
                ('wind_names', lambda x : x.split()),
                ('make_wind_plots', nsl.valid_bool),
                ('date0', nsl.valid_date),
                ('date1', nsl.valid_date),
                ]:
            setattr(self, att, None)
            if opts[att] != 'None':
                setattr(self, att, f_cvt(opts[att]))

        self.template = Template(os.path.join(
            self.data_dir, opts['basename_template']))
        self.mesh_file = os.path.join(os.getenv('NEXTSIM_MESH_DIR'),
            opts['mesh_file'])
        self.scalar_vars = []
        svars = opts.get('scalar_vars', 'None')
        if svars == 'None':
            return
        if isinstance(svars, str):
            svars = [svars]
        for s in svars:
            svar = s.split()
            varname, vmin, vmax = svar[:3]
            clabel = ' '.join(svar[3:])
            clim = None
            if vmin != 'None' and vmax != 'None':
                clim = [float(vmin), float(vmax)]
            self.scalar_vars += [(varname, clim, clabel)]

    @property
    def src_lonlat(self):
        f = self.template.safe_substitute(dict(varname=self.wind_names[0]))
        return mnu.nc_getinfo(f).get_lonlat(ij_range=self.ij_range)

    def set_grid_igi(self):
        print('Getting target grid')
        if not self.do_interp:
            # get target grid from eg file
            ncfil = self.template.safe_substitute(
                    dict(varname=self.temp_names[0]))
            kw = dict()
            if self.args.forcing == 'ec2_stere':
                kw = dict(spatial_dim_names=['x', 'y'], latlon=False)
            self.grid = Grid.init_from_netcdf(ncfil, projection=self.projection, **kw)
            return

        # get target grid from mesh file
        gmsh = GmshMesh(self.mesh_file, projection=self.projection)
        self.grid = gmsh.boundary.get_grid(
                resolution=self.plot_res*1000, projection=self.projection)
        print('Getting interpolator')
        self.igi = self.grid.get_interpolator(self.src_lonlat,
                interp_from=False, latlon=True)

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
        if self.date0 is not None and dto0 < self.date0:
            dto0 = self.date0
        dto2 = self.floor_time(nci.datetimes[-1], avg_period_hours=self.avg_period_hours)
        if self.date1 is not None and dto2 > self.date1:
            dto2 = self.date1
        while dto0 <= dto2:
            dto1 = dto0 + delt
            yield dto0, dto1, self.average_data(nci, varname, dto0, dto1)
            dto0 = dto1

    def plot_scalar(self, data, plot_var, dints=None):
        varname, clim, clabel = plot_var
        fig, ax = self.grid.plot(data, add_landmask=False, clim=clim,
                cmap='viridis', clabel=clabel)
        ax.coastlines(resolution='50m')

        # tidy up
        datestr1 = ''
        datestr2 = ''
        if dints is not None:
            datestr1 = '\n' + ' - '.join([dto.strftime('%Y-%m-%d %H:%M') for dto in dints])
            datestr2 = '_' + '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
        ax.set_title(f'{self.title}{datestr1}')

        #fig.show()
        figdir = os.path.join(self.outdir, varname)
        figname = os.path.join(figdir, f'{self.output_prefix}{varname}{datestr2}.png')
        print(f'Saving {figname}')
        os.makedirs(figdir, exist_ok=True)
        fig.savefig(figname)
        plt.close()

    def plot_wind(self, uv, dints):
        '''
        Parameters:
        -----------
        uv : list
            uv = [u, v] with u, v x/y or lon/lat components of wind velocity
        dints: list
            dints = [d0, d1] with d0,d1 datetime.datetime objects
            marking the start and finish of the averaging window
        '''
        spd = np.hypot(*uv)
        fig, ax = self.grid.plot(spd, add_landmask=False,
                cmap='viridis', clim=[0,10], clabel='Wind speed, m/s')
        ax.coastlines(resolution='50m')

        # wind speed contours
        #cs = ax.contour(*self.grid.xy, spd, levels=[10,15,20], colors='b')
        #ax.clabel(cs, inline=True, fontsize=10)

        # wind vectors
        dx_step = 100e3 #interval between wind vectors
        step = int(np.ceil(dx_step/self.grid.dx))
        x = self.grid.xy[0][::step,::step] 
        y = self.grid.xy[1][::step,::step]
        u = uv[0][::step,::step]
        v = uv[1][::step,::step]
        if self.do_interp:
            #u, v = nsl.rotate_lonlat2xy(self.projection, x, y, u, v)
            u, v = nsl.transform_vectors(
                    self.projection.pyproj, x, y, u, v)
        ax.quiver(x, y, u, v, units='xy', angles='xy', color='r')

        # tidy up
        datestr = '\n' + ' - '.join([
            dto.strftime('%Y-%m-%d %H:%M') for dto in dints])
        ax.set_title(f'{self.title}{datestr}')
        datestr = '_' + '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
        figdir = os.path.join(self.outdir, 'wind')
        figname = os.path.join(figdir, f'{self.output_prefix}wind{datestr}.png')

        #fig.show()
        print(f'Saving {figname}')
        os.makedirs(figdir, exist_ok=True)
        fig.savefig(figname)
        plt.close()

    def plot_wind_gradient(self, uv, dints):
        '''
        Parameters:
        -----------
        uv : list
            uv = [u, v] with u, v x/y or lon/lat components of wind velocity
        dints: list
            dints = [d0, d1] with d0,d1 datetime.datetime objects
            marking the start and finish of the averaging window
        '''
        u_x, v_x = [np.gradient(a, axis=1) for a in uv]
        u_y, v_y = [np.gradient(a, axis=0) for a in uv]
        curl = v_x - u_y, ('wcurl', None, 'Wind curl, s$^{-1}$')
        div = u_x + v_y, ('wdiv', None, 'Wind divergence, s$^{-1}$')
        shear = (np.hypot(u_x - v_y, u_y + v_x),
                    ('wshear', None, 'Wind shear, s$^{-1}$'))
        for grad in [curl, div, shear]:
            self.plot_scalar(*grad, dints=dints)

    def get_plot_data(self, arr):
        if not self.do_interp:
            return arr
        return self.igi.interp_field(arr)

    def test_plot(self):
        ''' test interpolation by plotting lon, lat '''
        #print(lonlat)
        plot_vars = [
                ('lon', None, 'Longitude, $^\circ$E'),
                ('lat', None, 'Latitude, $^\circ$N'),
                ]
        for plot_var, arr in zip(plot_vars, self.src_lonlat):
            self.plot_scalar(self.get_plot_data(arr), plot_var)

    def run(self):
        self.set_grid_igi()
        #self.test_plot()

        for plot_var in self.scalar_vars:
            #scalar fields
            varname = plot_var[0]
            print(f'Making plots for {varname}')
            f = self.template.safe_substitute(dict(varname=varname))
            for dto0, dto1, v in self.get_averaged_data(
                    mnu.nc_getinfo(f), varname):
                if v is None:
                    continue
                data = self.get_plot_data(v)
                if varname in self.temp_names:
                    data -= 273.15 #kelvin to deg C
                self.plot_scalar(data, plot_var, dints=(dto0, dto1))

        if self.make_wind_plots:
            #wind
            print(f'Making plots for wind speed')
            pairs = []
            for varname in self.wind_names:
                f = self.template.safe_substitute(dict(varname=varname))
                pairs += [(mnu.nc_getinfo(f), varname)]
            for (dto0, dto1, u10), (_, _, v10) in zip(
                    self.get_averaged_data(*pairs[0]),
                    self.get_averaged_data(*pairs[1]),
                    ):
                if u10 is None or v10 is None:
                    continue
                data = [self.get_plot_data(v) for v in [u10, v10]]
                self.plot_wind(data, (dto0, dto1))
                self.plot_wind_gradient(data, (dto0, dto1))

if __name__ == '__main__':
    obj = PlotAtmForcing()
    obj.run()
