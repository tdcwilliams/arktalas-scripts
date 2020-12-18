#! /usr/bin/env python
import os
import glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

from pylatex import (
        #Chapter,
        Command,
        Document,
        Figure,
        Foot,
        Head,
        LargeText,
        LineBreak,
        LongTabu,
        MediumText,
        MiniPage,
        MultiColumn,
        NewPage,
        PageStyle,
        Section,
        StandAloneGraphic,
        Tabu,
        Tabularx,
        TextColor,
        simple_page_number,
        )
from pylatex.math import Math
from pylatex.utils import bold, NoEscape

_GEOMETRY_OPTIONS = {
    "head": "40pt",
    "margin": "0.5in",
    "bottom": "0.6in",
    "includeheadfoot": True
}

_SECTIONS = {
        'eval-cs2smos' : 'Comparison to CS2-SMOS',
        'eval-osisaf-conc': 'Comparison to OSISAF concentration',
        'eval-osisaf-drift': NoEscape(r'Comparison to OSISAF drift ($\varepsilon<1.25$km/day)'),
        'eval-osisaf-drift-mu10kpd': NoEscape(r'Comparison to OSISAF drift ($\varepsilon<10$km/day)'),
        }

def parse_args():
    parser = ArgumentParser("Collect evaluation maps together")
    parser.add_argument('root_dir', type=str,
            help="root folder with evaluation files in subfolders eval-*")
    return vars(parser.parse_args())

class EvalPDF(Document):
    def __init__(self, root_dir, **kwargs):
        super().__init__(**kwargs)
        self.root_dir = os.path.abspath(root_dir)
        #add something to distinguish the name when comparing files in the same window of a PDF
        #reader
        subdir = os.path.basename(root_dir)
        self.pdf_file = os.path.join(root_dir, f'collected_maps_{subdir}')

    def parse_filename(self, figname):
        name1, name2 = os.path.split(figname)[1].split('-')
        date1 = dt.datetime.strptime(name1[-8:], '%Y%m%d')
        date2 = dt.datetime.strptime(name2[:8],  '%Y%m%d')
        return date1, date2

    def add_figs_to_page(self, fignames, new_page=True):
        if new_page:
            self.append(NewPage())
        ncols = len(fignames[0])
        gap = .03
        width = (1 - 2*gap)/ncols
        im_opts = f"width={width}"+r"\textwidth"
        d1 = self.parse_filename(fignames[0][0])[0].strftime('%B %Y')
        d2 = self.parse_filename(fignames[-1][-1])[1].strftime('%B %Y')
        with self.create(Figure(position='H')) as fig:
            with self.create(LongTabu(" ".join(ncols*["X[c]"]))) as fig_table:
                for figrow in fignames:
                    row = []
                    for figname in figrow:
                        row += [ StandAloneGraphic(figname, image_options=im_opts) ]
                    fig_table.add_row(row)
            fig.add_caption(f"{d1} - {d2}")

    def add_summary_figs(self, eval_dir):
        pattern = os.path.join(eval_dir, '*errors.png')
        print(f'Getting figures with pattern:\n\t{pattern}')
        figs_e = sorted(glob.glob(pattern))
        figs_e = [f for f in figs_e if 'mean' not in os.path.split(f)[1]]
        figs_s = [f.replace('errors', 'summary') for f in figs_e]
        new_page = False
        for fig_e, fig_s in zip(figs_e, figs_s):
            self.add_figs_to_page([[fig_e], [fig_s]], new_page=new_page)
            new_page = True

    def add_maps(self, eval_dir):
        # add maps
        d = dict()
        for k in ['bias', 'fcst', 'obs']:
            pattern = os.path.join(eval_dir, f'maps_{k}', '*mean*png')
            print(f'Getting figures with pattern:\n\t{pattern}')
            d[k] = sorted(glob.glob(pattern))
        nrows = 3

        figs = []
        for row in zip(d['bias'], d['fcst'], d['obs']):
            figs += [row]
            if len(figs) == nrows:
                self.add_figs_to_page(figs)
                figs = []

    def append_preamble(self):
        for cmd in [
                Command('usepackage', 'hyperref', options='linktoc=all'),
                Command('usepackage', 'float'),
                Command('usepackage', 'titlesec'),
                Command('usepackage', 'amsmath'),
                Command('author', 'Timothy Williams'),
                Command('title', 'Evaluation')
                ]:
            self.preamble.append(cmd)

    def intro(self):
        for cmd in [
                Command('maketitle'),
                #Command('tableofcontents'),
                ]:
            self.append(cmd)

    def run(self):
        self.append_preamble()
        self.intro()
        pattern = os.path.join(self.root_dir, 'eval-*')
        for eval_dir in sorted(glob.glob(pattern)):
            k = os.path.split(eval_dir)[1]
            label = k.replace('eval-', 'sec:').replace('-', '')
            sname = _SECTIONS.get(k, k)
            with self.create(
                    Section(sname, numbering=False, label=label)
                    ):
                self.add_summary_figs(eval_dir)
                self.add_maps(eval_dir)
            self.append(Command('pagebreak'))
        print(f'Creating {self.pdf_file}.pdf')
        self.generate_pdf(self.pdf_file, clean_tex=True)

if __name__ == "__main__":
    args = parse_args()
    pdf = EvalPDF(**args,
            documentclass='report',
            geometry_options=_GEOMETRY_OPTIONS)
    pdf.run()
