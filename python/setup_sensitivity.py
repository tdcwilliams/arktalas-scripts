#!/usr/bin/env python
import os
import glob
import shutil
import argparse
import numpy as np
import itertools
import datetime as dt
from string import Template
import subprocess
import json

from pynextsim.nextsim_config import NextsimConfig

class CustomTemplate(Template):
    delimiter = '%%'
    def __init__(self, filename):
        with open(filename, 'r') as fid:
            super().__init__(fid.read())

class LaunchSensitivity:

    def __init__(self, cli=None):
        args = self.parse_args(cli)
        self.root_dir = args.root_dir
        self.batch_name = args.batch_name
        self.executable = args.executable
        self.parse_config_file(args.config_file)
        self.ns_config_file = args.ns_config_file

    @staticmethod
    def parse_args(cli=None):
        ''' parse input arguments
        
        Parameters:
        -----------
        args : list(str)
            list of command line inputs

        Returns:
        --------
        parsed_args : argparse.ArgumentParser.namespace
            list of command line inputs
        '''
        parser = argparse.ArgumentParser(
                description="Setup nextsim sensitivity experiments")
        parser.add_argument('root_dir', type=str, help='root dir where experiments will be run')
        parser.add_argument('batch_name', type=str, help='name to give batch of experiments')
        parser.add_argument('config_file', type=str,
                help='cfg file with parameters to test')
        parser.add_argument('ns_config_file', type=str,
                help='reference nextsim cfg file')
        parser.add_argument('executable', type=str,
                help='path to nextsim executable')
        return parser.parse_args(cli)

    def parse_config_file(self, config_file):
        nc = NextsimConfig(config_file)
        # init slurm options
        self.slurm_opts = nc['slurm']
        opts = nc['launch_sensitivity']
        self.slurm_template = opts['slurm_template']
        # init sensitivity runs
        sopts = opts['sensitivity']
        if isinstance(sopts, str):
            sopts = [sopts]
        self.configs = []
        params = []
        values = []
        for s in sopts:
            sp = s.split()
            params += [sp[0]]
            values += [sp[1:]]
        for vals in itertools.product(*values):
            self.configs += [dict(zip(params, vals))]

    def write_ns_config(self, opts, fid):
        nc = NextsimConfig(self.ns_config_file)
        for k,v in opts.items():
            sec, opt = k.split('.')
            d = {opt: v}
            if sec not in nc:
                nc[sec] = d
            else:
                nc[sec].update(d)
        nc.write(fid)

    def setup_expt(self, i, i0):
        istr = '%.3i' %(i + i0)
        edir = os.path.join(self.root_dir, f'expt_{istr}')
        for n in ['inputs', 'tmp', 'logs']:
            subdir = os.path.join(edir, n)
            os.makedirs(subdir, exist_ok=True)
        # set the nextsim options
        opts = self.configs[i]
        print(f'\nSetting up {edir}\nwith options:')
        print(json.dumps(opts, indent=4))
        cfg = os.path.join(edir, 'inputs', 'nextsim.cfg')
        with open(cfg, 'w') as fid:
            self.write_ns_config(opts, fid)
        # get slurm template
        t = CustomTemplate(self.slurm_template)
        scr = os.path.join(edir, 'inputs', 'slurm.sh')
        opts = dict(**self.slurm_opts, job_name=f'sens_{self.batch_name}_{istr}')
        with open(scr, 'w') as fid:
            fid.write(t.substitute(opts))
        # get executable
        exe = os.path.join(edir, 'tmp', 'nextsim.exec')
        shutil.copyfile(self.executable, exe)
        return edir

    def run(self):
        os.makedirs(self.root_dir, exist_ok=True)
        runlist = os.path.join(self.root_dir, f'{self.batch_name}.csv')
        i0 = 0
        dlist = sorted(glob.glob(os.path.join(self.root_dir, 'expt_???')))
        if len(dlist) > 0:
            i0 = int(dlist[-1][-3:]) + 1
        with open(runlist, 'w') as fid:
            for i in range(len(self.configs)):
                edir = self.setup_expt(i, i0)
                if i == 0:
                    names = ['\"Experiment Directory\"'] + list(self.configs[i])
                    fid.write(','.join(names) + '\n')
                values = [os.path.basename(edir)] + list(self.configs[i].values())
                fid.write(','.join(values) + '\n')
        print(f"Saved experiment summary to {runlist}")

if __name__ == '__main__':
    obj = LaunchSensitivity()
    obj.run()
