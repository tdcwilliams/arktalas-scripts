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

from plot_atm_forcing import PlotAtmForcing

class PlotEcmwfFC(PlotAtmForcing):
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

if __name__ == "__main__":
    obj = PlotEcmwfFC()
    obj.run()
