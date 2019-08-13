import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy import constants as c
from ckscool.plot.config import *
import ckscool.io

sns.set_style('ticks')
sns.set_color_codes()
texteff = '$\mathregular{T}_{\mathregular{eff}}$'
texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 

cks1kw = dict(
    marker='.', lw=0, ms=7, color=sns.xkcd_rgb["dark sky blue"], 
    zorder=1, label='CKS-I', mew=1,mfc='none'
)

rstarticks = [0.5,0.7,1.0,1.4,2.0]

def fig_hr():
    sns.set_context('paper',font_scale=1.3)
    fig,axL = subplots(figsize=(7,5))

    df = ckscool.io.load_table('planets-cuts2+iso',cache=1)
    df = df[~df.isany]
    df = df.groupby('id_koi',as_index=False).nth(0)
    plot(df.cks_steff,df.gdir_srad,'.')
    xlabel(texteff)
    ylabel('Stellar Radius (Solar-Radii)')
    semilogy()
    xlim(7000,3500)
    minorticks_off()
    yticks(rstarticks,rstarticks)
    fig.set_tight_layout(True)

def fig_smet_smass():
    sns.set_context('paper',font_scale=1.3)
    fig,axL = subplots(figsize=(7,5))
    df = ckscool.io.load_table('planets-cuts2+iso',cache=1)
    df = df[~df.isany]
    df = df.groupby('id_koi',as_index=False).nth(0)
    semilogy()
    plot(df.cks_smet,df.giso_smass,'.')
    ylabel('Stellar Radius (Solar-Radii)')
    tight_layout(True)
    xlabel('[Fe/H] (dex)')
    yticks(rstarticks,rstarticks)
    minorticks_off()
    fig.set_tight_layout(True)


def fig_ferr_hist_star():
    fig = subplots(figsize=(4,3))
    df = ckscool.io.load_table('planets-cuts2+iso')
    df = df[~df.isany]
    df = df.groupby('id_koi',as_index=False).nth(0)

    bins = arange(0,0.1,0.0025)
    kw = dict(lw=2, bins=bins, histtype='step')

    rat = df.m17_kmag_err.dropna()
    hist(rat/2,label='$\sigma(m_K$) /2',**kw)

    rat = (df.gaia2_sparallax_err /  df.gaia2_sparallax).dropna()
    hist(rat,label='$\sigma(\pi) / \pi$',color='m',**kw)

    rat = 60 / df.cks_steff
    hist([2 * rat],label='2 $\sigma$(teff) / teff',**kw)

    rat = df.gdir_srad_err1/df.gdir_srad
    rat = rat.dropna().tolist()
    hist([rat],label='$\sigma(R_\star) / R_\star$',**kw)

    legend()
    xlabel('Fractional Precision')
    ylabel('Number of Stars')
    tight_layout()

def fig_ferr_hist_planet():
    fig = subplots(figsize=(4,3))
    df = ckscool.io.load_table('planets-cuts2+iso')

    bins = arange(0,0.2,0.005)
    kw = dict(lw=2, bins=bins, histtype='step')
    rat = df.gdir_srad_err1/df.gdir_srad
    rat = rat.dropna().tolist()
    hist([rat],label='$\sigma(R_\star) / R_\star$', **kw)

    rat = (df.koi_ror_err1 /df.koi_ror).dropna()
    hist(rat,label='$\sigma(R_p/R_\star) / (R_p/R_\star)$',color='m',**kw)

    rat = (df.gdir_prad_err1/df.gdir_prad).dropna()
    hist(rat,label='$\sigma(R_p) / (R_p)$',**kw)
    legend()
    xlabel('Fractional Precision')
    ylabel('Number of Stars')
    tight_layout()
