#!/usr/bin/env python
import os
import glob
from argparse import ArgumentParser
from collections import OrderedDict

import pandas as pd
from matplotlib import pylab as plt

import ckscool.io     # module for reading and writing datasets
import ckscool.value  # module for computing creating table
import ckscool.table  # module for computing scalar values for table
import ckscool._isoclassify
import ckscool.plot.sample   # submodule for including plots
import ckscool.plot.spectra
import ckscool.plot.hr
import ckscool.plot.planet
from ckscool.plot.compare import fig_compare
import ckscool.workflow

import numpy as np

import sys
#sys.path.append('../Kepler-Radius-Ratio/')
#import keprat.io

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('create-kbc', parents=[psr_parent],)
    psr2.set_defaults(func=create_kbc)

    psr2 = subpsr.add_parser('create-iso-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_batch)

    psr2 = subpsr.add_parser('create-iso-table', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_table)

    psr2 = subpsr.add_parser('create-chain-plots', parents=[psr_parent],)
    psr2.add_argument('id_koicand',type=str)
    psr2.set_defaults(func=create_chain_plots)

    psr2 = subpsr.add_parser('create-chain-plots-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_chain_plots_batch)

    psr2 = subpsr.add_parser('create-chain-summary', parents=[psr_parent],)
    psr2.set_defaults(func=create_chain_summary)

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

def create_kbc(args):
    import cpsutils.kbc
    df = cpsutils.kbc.loadkbc()
    b = ( df.type.str.contains('t') 
          & df.name.str.contains(r'^K\d{5}$|^CK\d{5}$') )
    df = df[b]
    df['id_koi'] = df.name.str.slice(start=-5).astype(int)
    df = df['obs name id_koi'.split()]
    namemap = {'obs':'id_obs','name':'id_name'}
    df = df.rename(columns=namemap)
    fn = ckscool.io.KBCFN
    df.to_csv()
    print "created {}".format(fn)

def create_iso_batch(args):
    sources = ['cks1','smsyn','smemp']
    modes = ['direct','grid-parallax-yes','grid-parallax-no']
    for source in sources:
        stars = ckscool._isoclassify.create_iso_batch_frames(source)
        for i,mode in enumerate(modes):
            fn = 'data/isoclassify-{}-{}.csv'.format(source,mode)
            stars[i].to_csv(fn)
            print "created {}".format(fn)

def create_iso_table(args):
    sources = ['cks1','smsyn','smemp']
    for source in sources:
        inpdir = os.path.join('isoclassify/',source)
        outcsv = os.path.join('data','isoclassify_{}.csv'.format(source))
        ckscool._isoclassify.create_iso_table(inpdir,outcsv)
    
def create_spectra_files(args):
    cut = ckscool.plot.spectra.get_spectra()
    fn = 'data/fig_spectra/fig_spectra-files.txt'
    cut.file.to_csv(fn,header=None,index=None) 
    print "created {}".format(fn)

def create_csv(args):
    df = ckscool.io.load_table(args.name)
    fn = 'data/{}.csv'.format(args.name)
    df.to_csv(fn)
    print "created {}".format(fn)

def create_chain_hdf(args):
    """
    Loop over every candidate in the thompson table and create a chain
    """

    t18 = keprat.io.load_table('t18')
    for i, row in t18.iterrows():
        table = 'chains-dr25-{}'.format(row.id_koicand)
        try:
            cachefn='data/dr25-mcmc-chains.hdf'
            df = keprat.io.load_table(table,cache=1,cachefn=cachefn)
        except IOError:
            with file('data/dr25-mcmc-chains_missing.txt','a') as f:
                s = table+" missing\n"
                print s
                f.write(s)

def create_chain_plots(args):
    import keprat
    from chainconsumer import ChainConsumer
    fmt = {}
    fmt['id_koicand'] = args.id_koicand
    table = "chains-dr25-{id_koicand:}".format(**fmt)
    cachefn='data/dr25-mcmc-chains.hdf'
    dr25 = keprat.io.load_table(table,cachefn=cachefn,cache=1)
    cols = "RHO ZPT EP1 PE1 BB1 RD1".split()
    parameters = [r"$\rho_{\star,\mathrm{circ}}$",r"$F_0$",r"$T_0$", r"$P$",r"$b$",r"$R_P / R_\star$"]
    c = ChainConsumer()
    chain = np.array(dr25[cols].iloc[::10])
    extents = zip(chain.min(axis=0),chain.max(axis=0))
    c.add_chain(chain, parameters=parameters)
    fn = 'mcmc/{id_koicand:}_distributions.png'.format(**fmt)
    c.plotter.plot_walks(parameters=parameters,extents=extents,filename=fn)
    print "created {}".format(fn)

def create_chain_plots_batch(args):
    cachefn='data/dr25-mcmc-chains.hdf'
    t18 = keprat.io.load_table('t18')
    t18.index = t18.id_koicand
    df = pd.read_csv('data/dr25-mcmc-chains_missing.txt',sep=' ',names=['s','x'])
    df['id_koicand'] = df.s.apply(lambda x : x.split('-')[-1]).str.upper()
    t18 = t18.drop(df.id_koicand)
    for i, row in t18.iterrows():
        print "run_ckscool.py create-chain-plots {id_koicand:}".format(**row)

def create_chain_summary(args):
    t18 = keprat.io.load_table('t18')
    df = []
    for i, row in t18.iterrows():
        try:
            d = keprat.io.get_summary(
                row.id_koicand,'dr25',cachefn='data/dr25-mcmc-chains.hdf'
            )
            df.append(d)
        except IOError:
            pass
        except AssertionError:
            pass

    df = pd.DataFrame(df)
    df.to_hdf('data/dr25-mcmc-summary.hdf','dr25')

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
    w.plot['cuts-planets-per-prad'] = ckscool.plot.sample.fig_cuts_planets_per_prad
    w.plot['cuts-planets-per-prad-zoom'] = lambda : ckscool.plot.sample.fig_cuts_planets_per_prad(zoom=True)

    #w.plot['compare-with-cks1'] = ckscool.plot.sample.fig_compare_with_cks1
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

    kw = dict(zoom=True,xerrfac=0.5,yerrfac=1.5)
    w.plot['planet-per-prad-smass-lo-zoom'] = lambda : ckscool.plot.planet.fig_per_prad(query='giso_smass < 0.75',**kw)
    w.plot['planet-per-prad-smass-mi-zoom'] = lambda : ckscool.plot.planet.fig_per_prad(query='0.75 < giso_smass < 0.95',**kw)
    w.plot['planet-per-prad-smass-hi-zoom'] = lambda : ckscool.plot.planet.fig_per_prad(query='0.95 < giso_smass',**kw)

    w.plot['planet-smass-prad'] = ckscool.plot.planet.fig_smass_prad
    w.plot['planet-smass-prad-nopoints'] = lambda : ckscool.plot.planet.fig_smass_prad(nopoints=True)
    w.plot['planet-smass-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_smass_prad(nopoints=True,zoom=True)
    w.plot['planet-smass-prad-zoom'] = lambda : ckscool.plot.planet.fig_smass_prad(zoom=True)

    w.plot['planet-sinc-prad'] = ckscool.plot.planet.fig_sinc_prad
    w.plot['planet-sinc-prad-nopoints'] = lambda : ckscool.plot.planet.fig_sinc_prad(nopoints=True)
    w.plot['planet-sinc-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_sinc_prad(nopoints=True,zoom=True)
    w.plot['planet-sinc-prad-zoom'] = lambda : ckscool.plot.planet.fig_sinc_prad(zoom=True)

    kw = dict(zoom=True,xerrfac=0.5,yerrfac=1.5)
    w.plot['planet-sinc-prad-smass-lo-zoom'] = lambda : ckscool.plot.planet.fig_sinc_prad(query='giso_smass < 0.75',**kw)
    w.plot['planet-sinc-prad-smass-mi-zoom'] = lambda : ckscool.plot.planet.fig_sinc_prad(query='0.75 < giso_smass < 0.95',**kw)
    w.plot['planet-sinc-prad-smass-hi-zoom'] = lambda : ckscool.plot.planet.fig_sinc_prad(query='0.95 < giso_smass',**kw)


    w.plot['planet-smet-prad'] = ckscool.plot.planet.fig_smet_prad
    w.plot['planet-smet-prad-nopoints'] = lambda : ckscool.plot.planet.fig_smet_prad(nopoints=True)
    w.plot['planet-smet-prad-nopoints-zoom'] = lambda : ckscool.plot.planet.fig_smet_prad(nopoints=True,zoom=True)
    w.plot['planet-smet-prad-zoom'] = lambda : ckscool.plot.planet.fig_smet_prad(zoom=True)
    w.plot['occur-contour-three'] = ckscool.plot.occur.fig_contour_three

    # table
    f1 = ckscool.table.tab_star
    w.table['star'] = f1
    w.table['star-stub'] = lambda: f1()[:5] + f1()[-5:]

    f2 = ckscool.table.tab_planet
    w.table['planet'] = f2
    w.table['planet-stub'] = lambda: f2()[:5] + f2()[-5:]

    # val
    w.val['stat'] = ckscool.value.val_stat

    # planet
    w.csv['star'] = ckscool.table.tab_star_csv
    w.csv['planet'] = ckscool.table.tab_planet_csv
    w.csv['planet-cuts2+iso'] = ckscool.table.tab_planet_full_csv
    w.csv['field-cuts'] = ckscool.table.tab_field_full_csv

    return w

if __name__=="__main__":
    main()

