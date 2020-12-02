#! /usr/bin/env python
import os
from string import Template

from plot_atm_forcing import PlotAtmForcing

class PlotEcmwfFC(PlotAtmForcing):
    def __init__(self):
        #should be inputs
        self.title = "ECMWF forecast"
        self.data_dir = '/cluster/projects/nn9624k/wrf_init/march2017/'
        self.template = Template(os.path.join(
            self.data_dir, 'od.ans.201302-201303.sfc.${varname}.nc'))
        self.ij_range = [0, 160, 0, 1440]
        self.plot_res = 10
        self.avg_period_hours = 12
        self.mesh_file = os.path.join(os.getenv('NEXTSIM_MESH_DIR'), 'wrf_arctic_10km.msh')
        self.outdir = 'figs/forcing/ecmwf_fc'
        self.scalar_vars = [
                ('t2m', [-40,0], '2-m air temperature, $^\circ$C'),
                #('sp', None, 'Sea level air pressure, Pa'),
                ]
        self.temp_names = ['t2m', 'd2m']
        self.wind_names = ['u10', 'v10']
        self.make_wind_plots = True
        self.date0 = None
        self.date1 = None

if __name__ == "__main__":
    obj = PlotEcmwfFC()
    obj.run()
