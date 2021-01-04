#! /usr/bin/env python
import os, sys
import glob
from matplotlib import pyplot as plt
import datetime as dt

from pynextsim.time_series import TimeSeries

from fig_base import FigBase

class CompRunsDrift(FigBase):

    def __init__(self, cli):
        self.args = self.parse_args(cli)
        cfg = self.read_config()
        self.run_params = cfg['run']
        self.src_params = cfg['src']

    def make_plots_one_source(self, src, leg_text_obs, obs_name, pattern, subdir):
        outdir = os.path.join(self.args.outdir, subdir)
        os.makedirs(outdir, exist_ok=True)

        fig1 = plt.figure(figsize=(12,4))
        ax1 = fig1.add_subplot(111)
        plot_info = None
        lines1 = []
        legend_texts1 = []
        figname = os.path.join(outdir, 'comp-bias-speed-%s.png' %src)

        fig2 = plt.figure(figsize=(12,4))
        ax2 = fig2.add_subplot(111)
        figname2 = os.path.join(outdir, 'comp-rmse-speed-%s.png' %src)
        lines2 = []
        legend_texts2 = []

        fig3 = plt.figure(figsize=(12,4))
        ax3 = fig3.add_subplot(111)
        figname3 = os.path.join(outdir, 'comp-vrmse-%s.png' %src)
        lines3 = []
        legend_texts3 = []

        y_factor = .5 #km/day instead of km/2days

        for fcdir, leg_text in self.run_params:
            pat_path = os.path.join(fcdir, pattern)
            flist = glob.glob(pat_path)
            if len(flist) == 0:
                print(f'\nEmpty pattern {pat_path}\n')
                continue
            f_use = flist[0]
            print(f'\nUsing {f_use}\n')
            TS = TimeSeries.init_from_file(f_use)

            if plot_info is None:
                for ax in [ax1, ax2, ax3]:
                    plot_info = TS.setup_time_series_plot(ax, label_xaxis_dates=True)
                x = plot_info['time_data']
                refdate = plot_info['refdate']

                uname = 'RMS_Drift_Uncertainty_O'
                if uname in TS.data:
                    yerr = TS.data[uname]

                    # add uncertainty to bias plot
                    lines1.append(ax1.fill_between(x, -y_factor*yerr, y_factor*yerr, facecolor=3*[.8]))
                    legend_texts1.append('RMS($\sigma_{OSISAF}$)')
                    ax1.plot(x, 0*x, ':k')

                    # add uncertainty to rmse plot
                    lines2.append(ax2.fill_between(x, 0*yerr, y_factor*yerr, facecolor=3*[.8]))
                    legend_texts2.append('RMS($\sigma_{OSISAF}$)')

                    # add uncertainty to vrmse plot
                    lines3.append(ax3.fill_between(x, 0*yerr, y_factor*yerr, facecolor=3*[.8]))
                    legend_texts3.append('RMS($\sigma_{OSISAF}$)')

            vname = 'Bias_Speed'
            lines1.append(self.plot_line(ax1, TS, refdate, vname, y_factor=y_factor))
            legend_texts1.append(leg_text)

            vname = 'RMSE_Speed'
            lines2.append(
                    self.plot_line(ax2, TS, refdate, vname, y_factor=y_factor))
            legend_texts2.append(leg_text)

            vname = 'RMSE'
            lines3.append(
                    self.plot_line(ax3, TS, refdate, vname, y_factor=y_factor))
            legend_texts3.append(leg_text)

        # Save figures
        ax1.set(ylabel='Drift speed bias, km/day')
        ax1.legend(lines1, legend_texts1, loc='center left', bbox_to_anchor=(1, .5))
        self.save_fig(fig1, ax1, figname)

        ax2.set(ylabel='Drift speed RMSE, km/day')
        ax2.legend(lines2, legend_texts2, loc='center left', bbox_to_anchor=(1, .5))
        self.save_fig(fig2, ax2, figname2)

        ax3.set(ylabel='Drift VRMSE, km/day')
        ax3.legend(lines3, legend_texts3, loc='center left', bbox_to_anchor=(1, .5))
        self.save_fig(fig3, ax3, figname3)

if __name__ == "__main__":
    obj = CompRunsDrift(sys.argv[1:])
    obj.run()
