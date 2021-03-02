#! /usr/bin/env python
import os, sys
import glob
from matplotlib import pyplot as plt
import datetime as dt
from pynextsim.time_series import TimeSeries

from fig_base import FigBase

class CompRunsScalars(FigBase):

    def make_plots_one_source(self, src, leg_text_obs, obs_name, pattern, subdir):
        outdir = os.path.join(self.args.outdir, subdir)
        os.makedirs(outdir, exist_ok=True)

        fig1 = plt.figure(figsize=(12,4))
        ax1 = fig1.add_subplot(111)
        plot_info = None
        figname = os.path.join(outdir, 'comp-mean-%s-%s.png' %(src, obs_name))

        fig2 = plt.figure(figsize=(12,4))
        ax2 = fig2.add_subplot(111)
        figname2 = os.path.join(outdir, 'comp-rmse-%s-%s.png' %(src, obs_name))

        fig3 = plt.figure(figsize=(12,4))
        ax3 = fig3.add_subplot(111)
        figname3 = os.path.join(outdir, 'comp-bias-%s-%s.png' %(src, obs_name))

        for fcdir, leg_text in self.run_params:
            pat_path = os.path.join(fcdir, pattern)
            flist = glob.glob(pat_path)
            if len(flist) == 0:
                print(f'\n[WARNING] Empty pattern {pat_path}\n')
                continue
            f_use = flist[0]
            print(f'Using {f_use}')
            TS = TimeSeries.init_from_file(f_use)

            if plot_info is None:
                for ax in [ax1, ax2, ax3]:
                    plot_info = TS.setup_time_series_plot(ax, label_xaxis_dates=True)
                x = plot_info['time_data']
                refdate = plot_info['refdate']

                vname = 'Mean%s_O' %obs_name
                y = TS.data[vname]
                units = TS.units[vname]
                if units != '':
                    units = ', %s' %units
                uname = 'RMS_%s_Uncertainty_O' %obs_name
                if uname in TS.data:
                    yerr = TS.data[uname]
                    ax1.fill_between(x, y-yerr, y+yerr, facecolor=3*[.8], label='Obs. error')

                    # add uncertainty to rmse plot
                    ax2.fill_between(x, 0*x, yerr, facecolor=3*[.8], label='Obs. error')

                    # add uncertainty to bias plot
                    ax3.fill_between(x, -yerr, yerr, facecolor=3*[.8], label='Obs. error')
                    ax3.plot(x, 0*x, ':k')

                # mean of observations
                ax1.plot(x, y, 'k', label=leg_text_obs)

            vname = 'Mean%s_F' %obs_name
            self.plot_line(ax1, TS, refdate, vname, plot_opts=dict(label=leg_text))

            vname = 'RMSE_F'
            self.plot_line(ax2, TS, refdate, vname, plot_opts=dict(label=leg_text))

            vname = 'Bias_F'
            self.plot_line(ax3, TS, refdate, vname, plot_opts=dict(label=leg_text))

        # Save figures
        ax1.set(ylabel='Mean %s%s' %(obs_name, units))
        ax1.legend()
        self.save_fig(fig1, ax1, figname)

        ax2.set(ylabel='RMSE %s%s' %(obs_name, units))
        ax2.legend()
        self.save_fig(fig2, ax2, figname2)

        ax3.set(ylabel='Bias %s%s' %(obs_name, units))
        ax3.legend()
        self.save_fig(fig3, ax3, figname3)

if __name__ == "__main__":
    obj = CompRunsScalars(sys.argv[1:])
    obj.run()
