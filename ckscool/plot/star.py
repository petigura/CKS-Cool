import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from ckscool.plot.config import *
import ckscool.io
from matplotlib.patches import Rectangle

sns.set_style('ticks')
sns.set_color_codes()

rstarticks = [0.5,0.7,1.0,1.4,2.0]


class Plotter():
    
    labels = {
        'steff' :'$T_{\mathregular{eff}}$ (K)',
        'srad' : 'Stellar Radius (Solar-Radii)',
        'smass' : 'Stellar Mass (Solar-Mass)',
        'smet' : '[Fe/H] (dex)',
        'sage' : 'Age (Gyr)'
    }

    limits = {
        'steff':(6750,3750),
        'srad': (0.4,2.4),
        'smet': (-0.6,0.6),
        'smass': (0.45,2.4),
        'sage': (0,15)
    }

    from matplotlib.colors import ListedColormap
    #cmap = ListedColormap(sns.color_palette("mako").as_hex())

    cmap = ListedColormap(sns.color_palette("viridis_r").as_hex())

    kw = dict(s=10,cmap=cmap,edgecolor='k',linewidths=0.5)
    
    def __init__(self):
        df = ckscool.io.load_table('planets-cuts2',cache=1)
        df = df[~df.isany]
        df = df.groupby('id_koi',as_index=False).nth(0)
        self.steff = df.cks_steff
        self.srad = df.gdir_srad
        self.smet = df.cks_smet
        self.sage = df.giso_sage
        self.slogage = df.giso_slogage
        self.slogage_err = df.eval('(giso_slogage_err1 - giso_slogage_err2)/2').tolist()
        self.smass = df.giso_smass

    def hr(self):
        scatter(self.steff,self.srad, c=self.slogage_err,vmin=0,vmax=0.5,**self.kw)
        semilogy()
        xlim()
        minorticks_off()
        yticks(rstarticks,rstarticks)
        self.setp('steff','srad')
        
    def smet_smass(self):
        semilogy()
        scatter(self.smet,self.smass, c=self.slogage,vmin=9,vmax=10.2,**self.kw)
        yticks(rstarticks,rstarticks)
        minorticks_off()
        self.setp('smet','smass')
        
    def smet_sage(self):
        scatter(self.smet,self.sage,c=self.smass, vmin=0.5,vmax=1.4, **self.kw)
        self.setp('smet','sage')
    
    def smass_sage(self):
        semilogx()
        scatter(self.smass,self.sage, c=self.slogage_err,vmin=0,vmax=0.5,**self.kw)
        self.setp('smass','sage')
        xticks(rstarticks,rstarticks)
        minorticks_off()
        
    def setp(self,xk,yk):
        setp(gca(),xlabel=self.labels[xk], ylabel=self.labels[yk],
             xlim=self.limits[xk],ylim=self.limits[yk])

def fig_sample():
    sns.set_context('paper',font_scale=1.0)
    fig,axL = subplots(ncols=2,nrows=2,figsize=(7,6))
    p = Plotter()

    # Metallicity and Mass
    sca(axL[0,0])
    p.smet_smass()
    fig_label('a')
    ax = axes([0.4, 0.85, 0.01, 0.1])
    cbar = colorbar(cax=ax,shrink=20)
    cbar.set_label('log(age)',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # HR diagram
    sca(axL[0,1])
    p.hr()
    fig_label('b')
    ax = axes([0.9, 0.85, 0.01, 0.1])
    cbar = colorbar(cax=ax,shrink=20)
    cbar.set_label('log(age) error',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Metallicity and age
    sca(axL[1,0])
    p.smet_sage()
    fig_label('c')
    ax = axes([0.4, 0.35, 0.01, 0.1])
    cbar = colorbar(cax=ax,shrink=20)
    cbar.set_label('mass',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Mass and age
    sca(axL[1,1])
    p.smass_sage()
    fig_label('d')
    ax = axes([0.9, 0.35, 0.01, 0.1])
    cbar = colorbar(cax=ax,shrink=20)
    cbar.set_label('log(age) error',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    tight_layout(True)



def fig_ferr_hist_star():
    fig = subplots(figsize=(4,3))
    df = ckscool.io.load_table('planets-cuts2')
    df = df[~df.isany]
    df = df.groupby('id_koi',as_index=False).nth(0)

    bins = arange(0,0.1,0.0025)
    kw = dict(lw=2, bins=bins, histtype='step')

    rat = df.m17_kmag_err.dropna() # magnitudes so don't divde
    hist(rat,label='$\sigma(m_K$)',**kw)

    rat = (df.gaia2_sparallax_err /  df.gaia2_sparallax).dropna()
    hist(rat,label='$\sigma(\pi) / \pi$',color='m',**kw)

    rat = df.cks_steff_err / df.cks_steff
    hist([rat],label='$\sigma$(teff) / teff',**kw)

    rat = df.gdir_srad_err1/df.gdir_srad
    rat = rat.dropna().tolist()
    hist([rat],label='$\sigma(R_\star) / R_\star$',**kw)


    legend()
    xlabel('Fractional Precision')
    ylabel('Number of Stars')
    tight_layout()

def fig_ferr_hist_planet():
    fig = subplots(figsize=(4,3))
    df = ckscool.io.load_table('planets-cuts2')

    bins = arange(0,20,0.5)
    kw = dict(lw=2, bins=bins, histtype='step')
    rat = df.gdir_srad_err1/df.gdir_srad *100
    rat = rat.dropna().tolist()
    hist([rat],label=r'$R_{\star}$', **kw)

    rat = (df.koi_ror_err1 /df.koi_ror).dropna() * 100
    hist(rat,label=r'$R_p/R_\star$',color='m',**kw)

    rat = (df.gdir_prad_err1/df.gdir_prad).dropna() *100
    hist(rat,label=r'$R_p$',**kw)
    legend()
    xlabel('Fractional Precision (%)')
    ylabel('Number of Planets')
    tight_layout()
