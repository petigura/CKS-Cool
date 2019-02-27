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
from ckscool.plot.compare import fig_compare
import ckscool.workflow

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('create-iso-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_batch)

    psr2 = subpsr.add_parser('create-iso-table', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_table)

    psr2 = subpsr.add_parser('create-csv', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_csv)

    psr2 = subpsr.add_parser('build', parents=[psr_parent], )
    psr2.add_argument('type',type=str)
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=build)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.add_argument('paper_dir',type=str)
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

## functions to build plots/tables/val for papers

def build(args):
    w = create_workflow()
    w.create_file(args.type, args.name) 

def update_paper(args):
    w = create_workflow()
    w.build_dir = 'build/'
    w.paper_dir = args.paper_dir
    w.update_paper()

def create_workflow():
    w = ckscool.workflow.Workflow()

    # register different plots here
    w.plot['cuts-kepmag-steff'] = ckscool.plot.sample.fig_cuts_kepmag_steff
    w.plot['cuts-period-prad'] = ckscool.plot.sample.fig_cuts_period_prad
    w.plot['cuts-smass-steff'] = ckscool.plot.sample.fig_cuts_smass_steff
    w.plot['cuts-stars-hr'] = ckscool.plot.sample.fig_cuts_stars_hr
    w.plot['compare-with-cks1'] = ckscool.plot.sample.fig_compare_with_cks1
    w.plot['ferr-hist-star'] = ckscool.plot.hr.fig_ferr_hist_star
    w.plot['ferr-hist-planet'] = ckscool.plot.hr.fig_ferr_hist_planet
    w.plot['star-steff-srad'] = ckscool.plot.hr.fig_hr
    w.plot['star-smet-smass'] = ckscool.plot.hr.fig_smet_smass
    w.plot['compare-ckscool-mann13'] = lambda : fig_compare('ckscool-mann13')
    w.plot['compare-ckscool-dressing13'] = lambda : fig_compare('ckscool-dressing13')
    w.plot['compare-ckscool-brewer18'] = lambda : fig_compare('ckscool-brewer18')

    w.plot['planet-per-prad'] = ckscool.plot.planet.fig_per_prad
    w.plot['planet-per-prad-nopoints'] = lambda : ckscool.plot.planet.fig_per_prad(nopoints=True)
    w.plot['planet-per-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_per_prad(nopoints=True,zoom=True)
    w.plot['planet-per-prad-zoom'] = lambda : ckscool.plot.planet.fig_per_prad(zoom=True)
    w.plot['planet-smass-prad'] = ckscool.plot.planet.fig_smass_prad
    w.plot['planet-smass-prad-nopoints'] = lambda : ckscool.plot.planet.fig_smass_prad(nopoints=True)
    w.plot['planet-smass-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_smass_prad(nopoints=True,zoom=True)
    w.plot['planet-smass-prad-zoom'] = lambda : ckscool.plot.planet.fig_smass_prad(zoom=True)
    w.plot['planet-smet-prad'] = ckscool.plot.planet.fig_smet_prad
    w.plot['planet-smet-prad-nopoints'] = lambda : ckscool.plot.planet.fig_smet_prad(nopoints=True)
    w.plot['planet-smet-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_smet_prad(nopoints=True,zoom=True)
    w.plot['planet-smet-prad-zoom'] = lambda : ckscool.plot.planet.fig_smet_prad(zoom=True)

    # val
    w.val['stat'] = ckscool.value.val_stat

    return w

if __name__=="__main__":
    main()

