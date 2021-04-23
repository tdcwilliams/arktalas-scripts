#! /usr/bin/env python
import os
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import datetime as dt
from subprocess import run
import cmocean

_PLOT_INFO = dict(
        rls = ('Upward longwave flux, W/m$^2$', 'rls_%Y%m%dT%H%M%SZ.png', 'viridis', [0,150]),
        sit = ('Thickness, m', 'sit_%Y%m%dT%H%M%SZ.png', 'viridis', [0,3]),
        sic = ('Concentration', 'sic_%Y%m%dT%H%M%SZ.png', cmocean.cm.ice, [0,1]),
        wspeed = ('Wind speed, m/s', 'wspeed_%Y%m%dT%H%M%SZ.png', 'viridis', [0,10]),
        )
rootdir_5km = '/cluster/work/users/timill/nextsim-stand-alone/wrf_arctic_5km/breakup_2013'
rootdir_10km = '/cluster/work/users/timill/nextsim-stand-alone/wrf_arctic_10km/breakup_2013'
_RUNS = [
        #('WRF 40km', f'{rootdir_10km}/expt_05_wrf_40km/outputs', 'figs/wrf_10km/wrf_40km'),
        ('WRF 20km', f'{rootdir_10km}/expt_06_wrf_20km/outputs', 'figs/wrf_10km/wrf_20km'),
        #('WRF 10km', f'{rootdir_10km}/expt_00_wrf10km/outputs', 'figs/wrf_10km/wrf_10km'),
        #('WRF 10km (Cd=0.004)', f'{rootdir_10km}/expt_01_wrf10km_cd40/outputs', 'figs/wrf_10km/wrf_10km_cd4'),
        #('ERA5', f'{rootdir_10km}/expt_02_era5/outputs', 'figs/wrf_10km/era5'),
        #('CFSR', f'{rootdir_10km}/expt_03_cfsr/outputs', 'figs/wrf_10km/cfsr'),
        #('WRF 10km (5km mesh): Clab=1.5MPa, Pmax=10kPa',
        #    f'{rootdir_5km}/expt_01_wrf_10km-C1.5/outputs', 'figs/wrf_5km/wrf_10km-C1.5'),
        #('WRF 10km (5km mesh): Clab=1.5MPa, Pmax=5kPa',
        #    f'{rootdir_5km}/expt_02_wrf_10km-C1.5-P5/outputs', 'figs/wrf_5km/wrf_10km-C1.5-P5'),
        ]

def plot(array, run_name, vname, dto, outdir):
    clabel, figname, cmap, clim = _PLOT_INFO[vname]
    figname = os.path.join(outdir, dto.strftime(figname))
    print(f'Saving {figname}')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(array, clim=clim, cmap=cmap, origin='lower')
    cbar = fig.colorbar(im, pad=0.01, shrink=0.75)#, format=format, ticks=cticks)
    cbar.set_label(clabel, rotation=270, labelpad=20, fontsize=16)
    ax.set_title(dto.strftime(f'{run_name}\n%Y-%m-%d %H:%M (12h average)'))
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(figname)
    plt.close()

def make_plots(run_name, expt_dir, outdir):
    moorings = os.path.join(expt_dir, 'Moorings.nc')
    with Dataset(moorings, 'r') as ds:
        for i, time in enumerate(ds.variables['time'][:]):
            #break
            os.makedirs(outdir, exist_ok=True)
            dto = dt.datetime(1900,1,1) + dt.timedelta(time)
            for vname in _PLOT_INFO:
                array = ds.variables[vname][i]
                plot(array, run_name, vname, dto, outdir)

    for vname in _PLOT_INFO:
        # make gif file
        tag = os.path.basename(outdir)
        cmd = ['convert', '+map', '-delay', '30', '-loop', '0',
                os.path.join(outdir, f'{vname}*.png'), os.path.join(outdir, f'anim_{vname}.{tag}.gif')]
        print('\n' + ' '.join(cmd))
        run(cmd)

for tup in _RUNS:
    make_plots(*tup)
