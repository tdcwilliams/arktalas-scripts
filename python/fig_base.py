import os, sys
import glob
from argparse import ArgumentParser
from matplotlib import pyplot as plt
from collections import defaultdict

class FigBase:

    def __init__(self, cli):
        self.args = self.parse_args(cli)
        self.cfg = self.read_config()
        self.run_params = self.cfg.get('run', None)
        self.src_params = self.cfg.get('src', None)

    @staticmethod
    def parse_args(cli):
        parser = ArgumentParser(description='Compare runs using scalar observations.')
        parser.add_argument('cfg_file', type=str, help='config file with info about runs')
        parser.add_argument('-o', '--outdir', type=str, default='output-comp-runs',
                            help='where to save the figs')
        return parser.parse_args(cli)

    def read_config(self):
        '''
        Parameters:
        -----------
        cfg_file : str
            name of config file
        Returns:
        --------
        runs: list
            each element is a list of strings
            [results_dir, linecolor, legend_text]
        sources: list
            each element is a list of strings
            [src_for_evaluate_forecast, legend_text, variable_name, sub_dir_name]
                legend text is for observation line
                sub_dir_name is directory path to evaluation results, relative to each run.results_dir)
        '''

        out = defaultdict(list)
        with open(self.args.cfg_file, 'r') as f:
            for lin in f:   
                if lin[0] == '#':
                    continue
                if '=' not in lin:
                    continue
                k, v = lin.split('=')
                out[k].append(v.split())
        #print(out)
        return out

    def plot_line(self, ax, TS, refdate, vname, plot_opts={}, y_factor=1):
        plot_info = TS.setup_time_series_plot(ax, refdate=refdate)
        x = plot_info['time_data']
        y = y_factor*TS.data[vname]
        return ax.plot(x, y, **plot_opts)[0]

    def save_fig(self, fig, ax, figname):
        print('Saving ' + figname)
        fig.tight_layout()
        #fig.subplots_adjust(wspace=0.5, top=.87, bottom=.2)
        fig.savefig(figname,
                #bbox_extra_artists=artists,
                bbox_inches='tight',
                )
        #plt.show(fig)
        ax.cla()
        plt.close(fig)

    def get_averages(self, dates, x, vec, stat_type):
        y0 = dates[0].year
        y1 = dates[-1].year
        dlims = [
                (dt.datetime(y0, 11, 1), dt.datetime(y0, 12, 31)),
                (dt.datetime(y1, 1, 1), dt.datetime(y1, 2, 28)),
                (dt.datetime(y1, 3, 1), dt.datetime(y1, 4, 30)),
                ]
        #linstyles = ['--', '-.', ':']
        linstyles = 3*[':']
        avgs = []
        i = 0
        for dto0, dto1 in dlims:
            in_rng = (dates>=dto0) * (dates<=dto1)
            assert(stat_type in ['Mean', 'Bias', 'RMSE'])
            if stat_type == 'RMSE':
                yav = np.sqrt(np.mean(vec[in_rng]**2))
            else:
                yav = np.mean(vec[in_rng])
            avgs.append(
                    [x[in_rng], yav + 0*x[in_rng], linstyles[i]])
            i += 1
        return avgs

    def run(self):
        for sp in self.src_params:
            self.make_plots_one_source(*sp)
