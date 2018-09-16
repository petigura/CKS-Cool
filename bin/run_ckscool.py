#!/usr/bin/env python
from argparse import ArgumentParser
import os
from collections import OrderedDict

import pandas as pd
import glob
import ckscool.io     # module for reading and writing datasets
import ckscool.value  # module for computing scalar values for table
import ckscool.table  # module for computing scalar values for table
import ckscool.plot.sample   # submodule for including plots
from matplotlib import pylab as plt

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

import ckscool.io
def create_iso_batch(args):
    """Create Isoclassify Batch Jobs

    Creates input parameters for two runs
    
       1. The direct method with the following constraints
          - teff, logg, fe, parallax, kmag
       2. The grid method with the following constraints
          - teff, logg, met, kmag [no parallax]

    We default to the direct method. But if the parallax method from
    the grid based method is significantly different than the gaia
    parallax, there is additional flux in the aperture which indicates
    dilution.

    """

    df = ckscool.io.load_table('ckscool+smemp')
    df = df.groupby('id_koi',as_index=False).nth(0)
    df = df.sort_values(by='id_koi')
    # Direct method
    df = df.rename(
        columns={
            'id_name':'id_starname',
            'parallax_error':'parallax_err',
            'ks_m':'kmag',
            'ks_msigcom':'kmag_err',
            'sm_fe':'feh',
            'sm_fe_err':'feh_err',
            'sm_teff':'teff',
            'sm_teff_err':'teff_err',
            'm17_kmag':'kmag',
            'm17_kmag_err':'kmag_err',
            'gaia2_sparallax':'parallax',
            'gaia2_sparallax_err':'parallax_err',
            'gaia2_ra':'ra',
            'gaia2_dec':'dec',
        }
    )

    df['kmag_err'] = df['kmag_err'].fillna(0.02)
    df['band'] = 'kmag'
    df['dust'] = 'green18'
    df['parallax'] /= 1e3
    df['parallax_err'] /= 1e3
    df0 = df.copy() 


    cols = [
        'id_starname','teff','teff_err','logg','logg_err','feh','feh_err',
        'parallax','parallax_err','kmag','kmag_err',
        'ra','dec','band','dust'
    ]

    # Direct method. Don't use spectroscopic logg values so as to not
    # pollute the parallax radii
    df = df0.copy()
    df['logg_err'] = 1 # use large uncertainties
    df['logg'] = 4.7 
    df = df[cols]
    fn = 'data/isoclassify-direct.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

    # Grid method with parallax. This will return model-dependent
    # values of Mstar, Rstar, age, density, luminosity
    df = df0.copy()
    df['logg_err'] = 1 # use large uncertainties
    df['logg'] = 4.7 
    df = df[cols]
    fn = 'data/isoclassify-grid-parallax-yes.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

    # Grid method. Don't set parallax so we can compare later
    df = df0.copy()
    df['logg_err'] = 1 # use large uncertainties
    df['logg'] = 4.7
    df = df[cols]
    df['parallax'] = -99
    df['parallax_err'] = 0
    fn = 'data/isoclassify-grid-parallax-no.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

def create_iso_table(args):
    """
    Read in isochrones csvfiles 
    Args:
        outdir (str): where to look for isochrones.csv files
    """
    import isoclassify.pipeline
    dfd = isoclassify.pipeline.scrape_csv('isoclassify/direct/*/*.csv')
    func = lambda x : x.split('.')[0]
    dfd['id_starname'] = dfd.id_starname.astype(str).apply(func)
    namemap = {
        'id_starname':'id_starname',
        'dir_rad':'gdir_srad',
        'dir_rad_err1':'gdir_srad_err1',
        'dir_rad_err2':'gdir_srad_err2',
    }
    dfd = dfd.rename(columns=namemap)[namemap.values()]
    
    fn = 'isoclassify/grid-parallax-yes/*/*.csv'
    dfg = isoclassify.pipeline.scrape_csv(fn)
    namemap = {
        'id_starname':'id_starname',
        'iso_mass':'giso_smass',
        'iso_mass_err1':'giso_smass_err1',
        'iso_mass_err2':'giso_smass_err2',
        'iso_rad':'giso_srad',
        'iso_rad_err1':'giso_srad_err1',
        'iso_rad_err2':'giso_srad_err2',
        'iso_rho':'giso_srho',
        'iso_rho_err1':'giso_srho_err1',
        'iso_rho_err2':'giso_srho_err2',
        'iso_age':'giso_sage',
        'iso_age_err1':'giso_sage_err1',
        'iso_age_err2':'giso_sage_err2',
    }
    dfg = dfg.rename(columns=namemap)[namemap.values()]
    dfg['id_starname'] = dfg.id_starname.astype(str).apply(func)

    fn = 'isoclassify/grid-parallax-no/*/*.csv'
    dfg2 = isoclassify.pipeline.scrape_csv(fn)
    temp = dfg2['id_starname'].copy()
    dfg2 = dfg2.drop(['id_starname'],axis=1)
    dfg2 = dfg2.convert_objects(convert_numeric=True)
    dfg2['id_starname'] = temp.astype(str)
    
    dfg2['giso2_sparallax'] = 1 / dfg2.iso_dis * 1e3
    dfg2['giso2_sparallax_err1'] = - dfg2['giso2_sparallax'] * dfg2['iso_dis_err2'] / dfg2['iso_dis']
    dfg2['giso2_sparallax_err2'] = - dfg2['giso2_sparallax'] * dfg2['iso_dis_err1'] / dfg2['iso_dis']
    columns = ['id_starname','giso2_sparallax','giso2_sparallax_err1', 'giso2_sparallax_err2',]
    dfg2 = dfg2[columns]
     
    dfm = pd.merge(dfd,dfg,on='id_starname')
    dfm = pd.merge(dfm,dfg2,on='id_starname')
    temp = dfm['id_starname'].copy()
    dfm = dfm.drop(['id_starname'],axis=1)
    dfm = dfm.convert_objects(convert_numeric=True)
    dfm['id_starname'] = temp.astype(str)
    dfm = ckscool.io.order_columns(dfm)
    
    fn = 'data/isoclassify_gaia2.csv'
    dfm.to_csv(fn)
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
        d['compare-with-cks1'] = ckscool.plot.sample.fig_compare_with_cks1
        # run_cksgaia create-plot #

        self.plot_dict = d

        d = OrderedDict()
        # register different tables here
        # d['population-stub'] = lambda : cksmet.tables.population(stub=True)
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
                elif name==key:
                    func()
                else:
                    continue
                    
                fn = self.key2fn(key, 'plot')
                plt.gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif name==key:
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
                    f.writelines("\n".join(lines))

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

