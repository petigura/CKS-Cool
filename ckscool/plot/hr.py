import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy import constants as c
from ckscool.plot.config import *
import ckscool.io
import cksgaia.io
sns.set_style('ticks')
sns.set_color_codes()
texteff = '$\mathregular{T}_{\mathregular{eff}}$'
texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 

def fig_hr():
    fig = figure(figsize=(6,4.5))
    df = cksgaia.io.load_table(
        'cksgaia-planets-filtered',
        cachefn='../CKS-Gaia/load_table_cache.hdf',cache=1
    )
    df = df.groupby('id_koi',as_index=False).nth(0)
    plot(df.cks_steff,df.gdir_srad,'.',label='CKS-I')

    df = ckscool.io.load_table('ckscool-stars-cuts',cache=1)
    df = df.groupby('id_koi',as_index=False).nth(0)
    #plot(df.sm_steff,df.gdir_srad,'.')
    df = df[~df.isany]

    plot(df.sm_steff,df.gdir_srad,'.',label='CKS-Cool')
    xlabel(texteff)
    ylabel('Stellar Radius (Solar-Radii)')
    semilogy()
    xlim(7000,3500)
    legend()

def fig_compare():
    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,
            'xtick.major.size':3.0,
            'xtick.direction': u'in',
            'ytick.direction': u'in',}
    )
    sns.set_context('paper',font_scale=1.1)
    df = cksmet.io.load_table('lamost-cks-calibration-sample')

    dfcal = df.copy()
    calfn  = 'cal_lamo-to-cks.fits'
    namemap = {'teff_new':'teff','logg_new':'logg','fe_new':'fe'}
    dfcal = dfcal.rename(columns=namemap)

    dfcal = cksmet._calibrate.calibrate(dfcal, calfn, mode='uncal')
    dfcal = dfcal.drop(['delta'],axis=1)

    fig = figure(figsize=(7,3.5))
    shape = (10,2)
    ax1 = subplot2grid(shape, (0,0), rowspan=8)
    ax2 = subplot2grid(shape, (8,0), rowspan=2, sharex=ax1)
    ax3 = subplot2grid(shape, (0,1), rowspan=8)
    ax4 = subplot2grid(shape, (8,1), rowspan=2, sharex=ax3)

    kw = dict(label1='CKS',label2='LAMOST',fig0=fig)
    comparison('smet',df.fe_lib,df.fe_new, axL0=[ax1,ax2],**kw)
    comparison('smet',dfcal.fe_lib,dfcal.fe, axL0=[ax3,ax4],**kw)
    fig.set_tight_layout(True)
    fig.set_tight_layout(False)
    fig.subplots_adjust(
        hspace=0.001,left=0.12,top=0.96,right=0.98,wspace=0.4,bottom=0.14
    )
