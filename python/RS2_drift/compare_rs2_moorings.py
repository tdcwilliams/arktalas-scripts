#! /usr/bin/env python
import os
from argparse import ArgumentParser
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt

from scipy.interpolate import RegularGridInterpolator

from pynextsim.projection_info import ProjectionInfo
from pynextsim.gridding import Grid
import mod_netcdf_utils as mnu

NS_PROJ = ProjectionInfo()

def get_arg_parser():
    p = ArgumentParser('Script to compare moorings with RS2-derived drift')
    p.add_argument('rs2_dir',
            help='path to directory with npz files from pattern matching')
    p.add_argument('moorings_file', help='path to moorings file')
    p.add_argument('outdir', help='Where to save results')
    p.add_argument('-t', '--test', action='store_true',
            help='Test script on 1 file only')
    p.add_argument('-f', '--force', action='store_true',
            help='Overwrite results')
    return p

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
    print(f'Number of drift vectors : {np.sum(gpi)}')
    return dto1, dto2, pm_results, xy1, gpi

class TrajectoryGenerator:
    def __init__(self, nci, projection, xy1, expansion_factor=2):
        self.nci = nci
        self.projection = projection
        self.expansion_factor = expansion_factor
        x, y = self.projection.pyproj(*nci.get_lonlat())
        xlim = self.get_dim_extent(xy1, 0)
        ylim = self.get_dim_extent(xy1, 1)
        self.x, j0, j1 = self.reduce_range(x[0], *xlim)
        self.y, i0, i1 = self.reduce_range(y[:,0].flatten(), *ylim)
        self.ij_range = [i0, i1+1, j0, j1+1]

    def get_dim_extent(self, xy1, i):
        av = np.mean(xy1[i])
        rng = xy1[i].max() - xy1[i].min()
        return (
                av - .5*self.expansion_factor*rng,
                av + .5*self.expansion_factor*rng,
                )

    @staticmethod
    def reduce_range(a, amin, amax):
        gpi = (a>=amin)*(a<=amax)
        a_ = a[gpi]
        amin_ = a_.min()
        amax_ = a_.max()
        lst = list(a)
        return a_, lst.index(amin_), lst.index(amax_)

    def time_iterator(self, dto1, dto2):

        # default time resolution
        dt_ref = self.nci.datetimes[1] - self.nci.datetimes[0]
        lbnds = np.array(self.nci.datetimes) - .5*dt_ref # lower bounds for averaging window
        ubnds = lbnds + dt_ref                      # upper bounds for averaging window
        dt_ref = dt_ref.total_seconds()

        # 1st time step info
        # - find last interval where dto1 >= the lower bound
        # - integrate until the next lower bound
        time_index1 = [dto1 >=  b for b in lbnds].index(False) - 1 # last True
        u = ubnds[time_index1]
        dt1 = (np.min([u, dto2]) - dto1).total_seconds()

        # Yield the 1st time index and integration time
        # - may not need any more
        yield time_index1, dt1

        if u < dto2:
            # last time step info
            time_index2 = [dto2 >  b for b in lbnds].index(False) - 1 #want last True
            dt2 = (dto2 - lbnds[time_index2]).total_seconds()

            for time_index in range(time_index1+1, time_index2):
                # Yield the intermediary time indices and integration times (if any)
                yield time_index, dt_ref

            # Yield the last time index and integration time
            yield time_index2, dt2

    def load_vars(self, time_index):
        data = dict()
        # load conc and u,v
        for vname in ['sic', 'siu', 'siv']:
            data[vname] = self.nci.get_var(
                    vname, ij_range=self.ij_range, time_index=time_index,
                    ).values.filled(np.nan)
        # mask open water in u,v
        mask = data['sic'] < .15
        for vname in ['siu', 'siv']:
            data[vname][mask] = np.nan
        return data

    def get_interpolated_vars(self, data, dst_p):
        src_p = [self.y, self.x]
        for vname in ['siu', 'siv']:
            rgi = RegularGridInterpolator(src_p, data[vname],
                    bounds_error=False, fill_value=np.nan)
            data[vname] = rgi(dst_p)
        return data

    def integrate_one_time_step(self, x0, y0, dt, data):
        dst_p = np.array([y0, x0]).T
        data = self.get_interpolated_vars(data, dst_p)
        return x0 + data['siu']*dt, y0 + data['siv']*dt

    def integrate_velocities(self, x0, y0, dto1, dto2):
        x = [np.array(x0)]
        y = [np.array(y0)]
        dtimes = [dto1]
        time_indices = []
        sic_av = 0
        for i, dt_i in self.time_iterator(dto1, dto2):
            data = self.load_vars(i)
            sic_av += data['sic']
            xi, yi = self.integrate_one_time_step(x[-1], y[-1], dt_i, data)
            x += [xi]
            y += [yi]
            dtimes += [dtimes[-1] + dt.timedelta(seconds=dt_i)]
            time_indices += [i]
        assert(dtimes[-1] == dto2)
        sic_av /= len(time_indices)
        return np.array(x).T, np.array(y).T, dtimes, time_indices, sic_av

    def get_grid(self):
        return Grid(*np.meshgrid(self.x, self.y),
                projection=self.projection)

def save_fig(fig, filename):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    print(f'Saving {filename}')
    fig.savefig(filename, bbox_inches='tight')
    plt.close()

def save_npz(filename, **kwargs):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    print(f'Saving {filename}')
    np.savez(filename, **kwargs)

def load_npz(args, filename):
    if os.path.exists(filename) and not args.force:
        print(f'Loading {filename}')
        return dict(np.load(filename))

def compare(dx_obs, dy_obs, dx_mod, dy_mod, delta_t, **kwargs):
    unit_fac = 24*3600*1e-3 #m/s to km/day
    errors = dict()
    # bias in speed
    spd_obs = unit_fac*np.hypot(dx_obs, dy_obs)/delta_t
    spd_mod = unit_fac*np.hypot(dx_mod, dy_mod)/delta_t
    errors['bias_speed'] = np.nanmean(spd_mod - spd_obs)
    print(f"Bias in speed = {errors['bias_speed']} km/day")
    # RMSE in speed
    errors['rmse_speed'] = np.sqrt(np.nanmean((spd_mod - spd_obs)**2))
    print(f"RMSE in speed = {errors['rmse_speed']} km/day")
    # Vector RMSE
    vdiff = unit_fac*np.hypot(dx_mod - dx_obs, dy_mod - dy_obs)/delta_t
    errors['vrmse'] = np.sqrt(np.nanmean(vdiff**2))
    print(f"VMRSE = {errors['vrmse']} km/day")
    return errors

def process_1file(args, rs2_file):

    base = os.path.basename(rs2_file).replace('pm', 'comp_moorings')
    npz_file = os.path.join(args.outdir, base)

    # try to load results from file
    # - otherwise create them
    results = load_npz(args, npz_file)
    if results is not None:
        errors = dict()
        for k in ['bias_speed', 'rmse_speed', 'vrmse']:
            errors[k] = results.pop(k)
    else:
        results = dict()
        # read RS2 results
        dto1, dto2, pm_results, xy1, gpi_rs2 = read_rs2_file(rs2_file) 
        results['dx_obs'] = pm_results['upm_clean'][gpi_rs2]
        results['dy_obs'] = pm_results['vpm_clean'][gpi_rs2]
        results['x_obs'], results['y_obs'] = xy1
        results['delta_t'] = (dto2 -dto1).total_seconds()

        # process neXtSIM results
        nci = mnu.nc_getinfo(args.moorings_file)
        tg = TrajectoryGenerator(nci, NS_PROJ, xy1, expansion_factor=1.4)
        (xt, yt, dtimes, time_indices, results['sic_av'],
                ) = tg.integrate_velocities(*xy1, dto1, dto2)
        results['dx_mod'] = xt[:,-1] - xt[:,0]
        results['dy_mod'] = yt[:,-1] - yt[:,0]

        # compare mean differences
        errors = compare(**results)

        # save results
        results['x_grid'], results['y_grid'] = tg.get_grid().xy
        save_npz(npz_file, **results, **errors)

    # plot
    plot(args, rs2_file, **results)

    return errors

def plot(args, rs2_file, x_obs, y_obs, dx_obs, dy_obs,
        dx_mod, dy_mod, x_grid, y_grid, sic_av, **kwargs):
    grid = Grid(x_grid, y_grid)
    fig, ax = grid.plot(sic_av, clabel='neXtSIM concentration')
    s = slice(None, None, 2) #plot every 2nd vector
    ax.quiver(x_obs[s], y_obs[s], dx_obs[s], dy_obs[s],
            color='r', units='xy', scale=.5, label='RS2')
    ax.quiver(x_obs[s], y_obs[s], dx_mod[s], dy_mod[s],
            color='g', units='xy', scale=.5, label='NS')
    ax.legend()
    os.makedirs(args.outdir, exist_ok=True)
    figname = os.path.join(
        args.outdir,
        os.path.basename(rs2_file).replace('npz', 'png')
        )
    save_fig(fig, figname)

def run():
    args = get_arg_parser().parse_args()
    if args.test:
        # if testing just run 1 example
        rs2_file = os.path.join(args.rs2_dir,
                'pm_20130224T023841Z-20130225T020927Z_21.npz')
        process_1file(args, rs2_file)
        return

    for rs2_file in glob.glob(pattern):
        try:
            errors_i = process_1pair(rs2_file)
        except:
            continue


if __name__ == '__main__':
    run()
