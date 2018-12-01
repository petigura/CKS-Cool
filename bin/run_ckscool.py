#!/usr/bin/env python
import os
import glob
from argparse import ArgumentParser
from collections import OrderedDict

import pandas as pd
from matplotlib import pylab as plt

import ckscool.io     # module for reading and writing datasets
import ckscool.value  # module for computing scalar values for table
import ckscool.table  # module for computing scalar values for table
import ckscool._isoclassify
import ckscool.plot.sample   # submodule for including plots
import ckscool.plot.spectra
import ckscool.plot.hr
import ckscool.plot.planet

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('create-iso-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_batch)

    psr2 = subpsr.add_parser('create-iso-table', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_table)

    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-csv', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_csv)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)


def create_iso_batch(args):
    ckscool._isoclassify.create_iso_batch(args)

def create_iso_table(args):
    ckscool._isoclassify.create_iso_table(args)

def create_spectra_files(args):
    cut = ckscool.plot.spectra.get_spectra()
    fn = 'data/fig_spectra/fig_spectra-files.txt'
    cut.file.to_csv(fn,header=None,index=None) 
    print "created {}".format(fn)

def create_csv(args):
    df = ckscool.io.load_table('ckscool-planets-cuts')
    fn = 'data/ckscool-planets-cuts.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

def create_table(args):
    w = Workflow()
    w.create_file('table', args.name ) 

def create_plot(args):
    w = Workflow()
    w.create_file('plot', args.name ) 

def create_val(args):
    w = Workflow()
    w.create_file('val',args.name) 

def update_paper(args):
    w = Workflow()
    w.update_paper()

class Workflow(object):
    def __init__(self):
        d = OrderedDict()
        # register different plots here
        d['cuts-kepmag-steff'] = ckscool.plot.sample.fig_cuts_kepmag_steff
        d['cuts-period-prad'] = ckscool.plot.sample.fig_cuts_period_prad
        d['cuts-smass-steff'] = ckscool.plot.sample.fig_cuts_smass_steff
        d['cuts-stars-hr'] = ckscool.plot.sample.fig_cuts_stars_hr
        d['cuts-planets-per-prad'] = ckscool.plot.sample.fig_cuts_planets_per_prad
        d['compare-with-cks1'] = ckscool.plot.sample.fig_compare_with_cks1
        d['ferr-hist-star'] = ckscool.plot.hr.fig_ferr_hist_star
        d['ferr-hist-planet'] = ckscool.plot.hr.fig_ferr_hist_planet
        d['planet-per-prad'] = ckscool.plot.planet.fig_per_prad
        d['planet-smass-prad'] = ckscool.plot.planet.fig_smass_prad
        d['planet-smet-prad'] = ckscool.plot.planet.fig_smet_prad

        d['star-steff-srad'] = ckscool.plot.hr.fig_hr
        d['star-smet-smass'] = ckscool.plot.hr.fig_smet_smass

        from ckscool.plot.compare import fig_compare
        d['compare-ckscool-mann13'] = lambda : fig_compare('ckscool-mann13')
        d['compare-ckscool-dressing13'] = lambda : fig_compare('ckscool-dressing13')
        d['compare-ckscool-brewer18'] = lambda : fig_compare('ckscool-brewer18')

        # run_cksgaia create-plot #

        self.plot_dict = d

        d = OrderedDict()
        # register different tables here
        d['star'] = lambda : ckscool.table.tab_star() 
        d['star-stub'] = lambda : ckscool.table.tab_star()[:10]
        d['planet'] = lambda : ckscool.table.tab_planet()
        d['planet-stub'] = lambda : ckscool.table.tab_planet()[:10]
        self.table_dict = d

        d = OrderedDict()
        d['stat'] = ckscool.value.val_stat
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return 'fig_'+key+'.pdf'
        if kind=='table':
            return 'tab_'+key+'.tex'
        if kind=='val':
            return 'val_'+key+'.tex'
            
    def create_file(self, kind, name):
        i = 0
        for key, func in self.all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif key.count(name)==1:
                    func()
                else:
                    continue
                    
                fn = self.key2fn(key, 'plot')
                plt.gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif key.count(name)==1:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'table')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{",
                ]

                lines2 = [
                    "}[XX]",
                    "}"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("%\n".join(lines))

            i+=1

        if i==0:
            assert False, name + " not a valid key"

    def update_paper(self):
        for kind, d in self.all_dict.iteritems():
            for key, val in d.iteritems():
                fn = self.key2fn(key, kind)
                cmd = 'cp {} paper/'.format(fn)
                print cmd
                os.system(cmd)

if __name__=="__main__":
    main()

