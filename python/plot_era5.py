#! /usr/bin/env python
import os
from string import Template
import datetime as dt

from plot_atm_forcing import PlotAtmForcing

class PlotEra5(PlotAtmForcing):
    def __init__(self):
        #should be inputs
        self.title = 'ERA5'
        self.data_dir = '/cluster/projects/nn2993k/sim/data/ERA5'
        self.template = Template(os.path.join(
            self.data_dir, 'ERA5_${varname}_y2013.nc'))
        self.ij_range = None
        self.plot_res = 10
        self.avg_period_hours = 12
        self.mesh_file = os.path.join(os.getenv('NEXTSIM_MESH_DIR'), 'wrf_arctic_10km.msh')
        self.outdir = 'figs/forcing/era5'
        self.output_prefix = 'era5_'
        self.scalar_vars = [
                ('t2m', [-40,0], '2-m air temperature, $^\circ$C'),
                ]
        self.temp_names = ['t2m', 'd2m']
        self.wind_names = ['u10', 'v10']
        self.make_wind_plots = True
        self.date0 = dt.datetime(2013,2,1)
        self.date1 = dt.datetime(2013,3,1)

if __name__ == "__main__":
    obj = PlotEra5()
    obj.run()
