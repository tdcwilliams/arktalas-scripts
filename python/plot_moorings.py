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

from swarp_funs import mod_netcdf_utils as mnu


PROJECTIONS = dict(
        nextsim=ProjectionInfo(),
        arcmfc=ProjectionInfo.topaz_np_stere(),
        )


def parse_args(cli=None):
    """
    parse command line arguments

    Parameters:
    -----------
    cli : list(str)
        list of command line arguments
        default is None, which means they are taken from sys.argv[1:]

    Returns:
    --------
    args : dict
    """
    parser = argparse.ArgumentParser("script to plot moorings file")
    parser.add_argument('filename', type=str,
            help='path to moorings file to plot')
    parser.add_argument('config_file', type=str,
            help='path to config file with input settings')
    parser.add_argument('-o', '--outdir', type=str, default=None,
            help="where to save the figures (default is in the same directory as the moorings file")
    parser.add_argument('-fp', '--figname-prefix', type=str, default="nextsim_",
            help="prefix to figure names")
    parser.add_argument('-d0', '--date0', type=nsl.valid_date, default=None,
            help="plot only dates after this one")
    parser.add_argument('-d1', '--date1', type=nsl.valid_date, default=None,
            help="plot only dates before this one")
    return vars(parser.parse_args(cli))


class PlotMoorings:
    """ Base class for plotting moorings files """

    def __init__(self, filename, config_file, outdir=None, figname_prefix='nextsim_',
            date0=None, date1=None):
        """
        Parameters:
        -----------
        filename : str
            help='path to moorings file to plot
        config_file : str
            path to config file with input settings
        outdir : str
            where to save the figures (default is in the same directory as the moorings file
        figname-prefix : str
            prefix to figure names
        date0 : datetime.datetime
            plot only dates after this one
        date1 : datetime.datetime
            plot only dates before this one
        """
        self.filename = filename
        self.config_file = config_file
        self.outdir = outdir
        self.figname_prefix = figname_prefix
        self.date0 = date0
        self.date1 = date1

        self.parse_config_file()
        # if self.grid_projection is not known, interp from lon, lat to nextsim stereographic projection
        self.projection = PROJECTIONS.get(self.grid_projection, ProjectionInfo())
        self.do_interp = (self.grid_projection not in PROJECTIONS)
        self.igi = None


    def parse_var_list(self, opts, vtype):
        """
        Extract list of variables to plot from config file

        Parameters:
        -----------
        opts : NextsimConfig or dict
            result of reading the config file
        vtype : str
            'scalar_vars' or 'vector_vars' 

        Returns:
        --------
        var_list : list
            each element is a tuple
            (varname, cmap, clim, clabel)
            varname = name of variable in the netcdf file
            clim = (vmin, vmax) = range for the colorbar
            clabel = label for the colorbar
        """
        tmp_vars = opts.get(vtype, 'None')
        if tmp_vars == 'None':
            return []
        if isinstance(tmp_vars, str):
            # only one variable
            tmp_vars = [tmp_vars]

        def split_string(v, vtype):
            sv = v.split()
            if vtype == 'scalar_vars':
                lst = [sv.pop(0)]
            else:
                lst = [sv.pop(0) for _ in range(2)]
            cmap, vmin, vmax = sv[:3]
            clim = None
            if vmin != 'None' and vmax != 'None':
                clim = [float(vmin), float(vmax)]
            clabel = ' '.join(sv[3:])
            return lst + [cmap, clim, clabel]

        return [split_string(v, vtype) for v in tmp_vars]


    def parse_config_file(self):
        """
        parse the config file
        """
        opts = NextsimConfig(self.config_file)['plot_moorings']
        for att, f_cvt in [
                ('title', str),
                ('ij_range', lambda x : [int(i) for i in x.split()]),
                ('grid_projection', str),
                ('plot_resolution', float),
                ('avg_period_hours', int),
                ('plot_deformation', nsl.valid_bool),
                ]:
            setattr(self, att, None)
            if opts[att] != 'None':
                setattr(self, att, f_cvt(opts[att]))
        self.mesh_file = os.path.join(os.getenv('NEXTSIM_MESH_DIR'),
            opts['mesh_file'])
        self.scalar_vars = self.parse_var_list(opts, 'scalar_vars')
        self.vector_vars = self.parse_var_list(opts, 'vector_vars')


    def set_grid_igi(self):
        print('Getting target grid')
        if not self.do_interp:
            # get target grid from eg file
            kw = dict()
            #kw = dict(spatial_dim_names=['x', 'y'], latlon=False)
            self.grid = Grid.init_from_netcdf(self.filename, projection=self.projection, **kw)
            return

        # get target grid from mesh file
        gmsh = GmshMesh(self.mesh_file, projection=self.projection)
        self.grid = gmsh.boundary.get_grid(
                resolution=self.plot_resolution*1000, projection=self.projection)
        print('Getting interpolator')
        self.igi = self.grid.get_interpolator(self.src_lonlat,
                interp_from=False, latlon=True)

    @staticmethod
    def floor_time(dto, avg_period_hours=12):
        """
        Take the average of the data between a range of datetimes

        Parameters:
        -----------
        dto : datetime.datetime
            earlier datetime
        avg_period_hours : int
            hours to average over

        Returns:
        ---------
        data : numpy.ndarray
            missing values are filled with NaN
        """
        h = avg_period_hours*(dto.hour//avg_period_hours)
        return dt.datetime(dto.year, dto.month, dto.day, h)


    def average_data(self, nci, varname, dto0, dto1):
        """
        Take the average of the data between a range of datetimes

        Parameters:
        -----------
        nci : nc_getinfo
        varname : str
        dto0 : datetime.datetime
            earlier datetime
        dto1 : datetime.datetime
            later datetime

        Returns:
        ---------
        data : numpy.ndarray
            missing values are filled with NaN
        """
        dtimes = np.array(nci.datetimes)
        dtimes = dtimes[(dtimes>=dto0)*(dtimes<=dto1)]
        n = len(dtimes)
        if n == 0:
            return
        data = 0.
        for dto in dtimes:
            data += nci.get_var(varname,
                    time_index=nci.datetimes.index(dto), ij_range=self.ij_range).values/n
        return data.filled(np.nan)


    def get_averaged_data(self, nci, varname):
        """
        Loop over the dates in the file, averaging over intervals of length self.avg_period_hours

        Parameters:
        -----------
        nci : nc_getinfo
        varname : str

        Returns:
        ---------
        avg_data : generator
            each element is (dto0, dto1, data), with
                dto0 : datetime.datetime
                    earlier datetime
                dto1 : datetime.datetime
                    later datetime
                data : numpy.ndarray
                    missing values are filled with NaN
        """
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


    def finish_and_save_fig(self, ax, varname, dints=None):
        """
        set title and save figure

        Parameters:
        -----------
        ax : matplotlib.axes._subplots.AxesSubplot
            plot axes
        varname : variable name
        dints : tuple(datetime.datetime)
            (dto0, dto1) time range to use in the title or the figure name
        """

        # date strings
        datestr1 = '' # for title
        datestr2 = '' # for figname
        if dints is not None:
            datestr1 = '\n' + ' - '.join([dto.strftime('%Y-%m-%d %H:%M') for dto in dints])
            datestr2 = '_' + '-'.join([dto.strftime('%Y%m%dT%H%M%SZ') for dto in dints])
        ax.set_title(f'{self.title}{datestr1}')

        # save fig
        figname = os.path.join(self.outdir, f'{self.figname_prefix}{varname}{datestr2}.png')
        print(f'Saving {figname}')
        os.makedirs(self.outdir, exist_ok=True)
        ax.figure.savefig(figname)
        plt.close()


    def plot_scalar(self, data, varname, cmap, clim, clabel, dints=None):
        """
        Plot scalar variable.

        Parameters:
        -----------
        data : numpy.ndarray
            array to plot
        cmap : str
            colormap to use
        clim : tuple(str)
            colorbar range
        clabel : str
            colorbar label
        dints : tuple(datetime.datetime)
            (dto0, dto1) time range to use in the title or the figure name
        """
        fig, ax = self.grid.plot(data, cmap=cmap, clim=clim, clabel=clabel)
        ax.coastlines(resolution='50m')
        self.finish_and_save_fig(ax, varname, dints=dints)


    def plot_vector(self, uv, varname, cmap, clim, clabel, dints):
        '''
        Plot vector variable.

        Parameters:
        -----------
        uv : tuple(numpy.ndarray)
            uv = (u, v) with u, v x/y or lon/lat components of vector
        cmap : str
            colormap to use
        clim : tuple(str)
            colorbar range
        clabel : str
            colorbar label
        dints : tuple(datetime.datetime)
            (dto0, dto1) time range to use in the title or the figure name
        '''
        spd = np.hypot(*uv)
        fig, ax = self.grid.plot(spd, cmap=cmap, clim=clim, clabel=clabel)
        ax.coastlines(resolution='50m')

        # add some direction vectors
        step = int(np.ceil(np.max(self.grid.shape)/50)) # if more than 50 vectors, increase the step between them
        x = self.grid.xy[0][::step,::step] 
        y = self.grid.xy[1][::step,::step]
        u = uv[0][::step,::step]
        v = uv[1][::step,::step]
        if self.do_interp:
            u, v = nsl.transform_vectors(
                    self.projection.pyproj, x, y, u, v)
        ax.quiver(x, y, u, v, units='xy', angles='xy', color='r')

        # add title and save fig
        self.finish_and_save_fig(ax, varname, dints=dints)


    def plot_deformation(self, uv, dints):
        '''
        Plot deformation

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
        clim = (0,10)
        div = (u_x + v_y,
                ('div', 'viridis', clim, r'Sea ice divergence, $\%\textrm{day}^{-1}$'))
        shear = (np.hypot(u_x - v_y, u_y + v_x),
                    ('shear', clim, r'Sea ice shear deformation, $\%\textrm{day}^{-1}$'))
        for grad in [div, shear]:
            self.plot_scalar(*grad, dints=dints)


    def transform_data(self, arr):
        """
        Transform the data to the plot coordinates (do nothing if already in the
        right coordinates ie self.do_interp=False).

        Parameters:
        -----------
        input : numpy.ndarray
            array in original coordinates

        Parameters:
        -----------
        output : numpy.ndarray
            array in plot coordinates
        """
        if not self.do_interp:
            return arr
        return self.igi.interp_field(arr)


    def plot_scalars(self):
        """ make all the scalar plots """

        for plot_var in self.scalar_vars:
            #scalar fields
            varname = plot_var[0]
            print(f'Making plots for {varname}')
            for dto0, dto1, v in self.get_averaged_data(
                    mnu.nc_getinfo(self.filename), varname):
                if v is None:
                    continue
                data = self.transform_data(v)
                self.plot_scalar(data, *plot_var, dints=(dto0, dto1))


    def plot_vectors(self):
        """ make all the vector plots """

        def get_varname(un, vn):
            for i in range(len(un)):
                if un[:i] == vn[:i]:
                    lst = [un[:i], un[i:], vn[i:]]
                else:
                    break
            return ''.join(lst)

        for uname, vname, cmap, clim, clabel in self.vector_vars:
            print(f'Making vector plots for ({uname},{vname})')
            nci = mnu.nc_getinfo(self.filename)
            for (dto0, dto1, u), (_, _, v) in zip(
                    self.get_averaged_data(nci, uname),
                    self.get_averaged_data(nci, vname),
                    ):
                if u is None or v is None:
                    continue
                data = [self.transform_data(a) for a in [u, v]]
                self.plot_vector(data, get_varname(uname, vname), cmap, clim, clabel, (dto0, dto1))
                if uname == "siu" and self.plot_deformation:
                    self.plot_deformation(data, (dto0, dto1))


    def run(self):
        """ main method to make all the plots """
        self.set_grid_igi()
        self.plot_scalars()
        self.plot_vectors()


if __name__ == '__main__':
    PlotMoorings(**parse_args()).run()
