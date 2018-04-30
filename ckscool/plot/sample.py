import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy import constants as c
from ckscool.plot.config import *

errorbar_kw = dict(fmt='.',markersize=5,color='b')

sns.set_style('whitegrid')
sns.set_color_codes()

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
texteff = '$\mathregular{T}_{\mathregular{eff}}$'
texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
perticks = [0.3,1,3,10,30,100,300]

letter_bbox_props = dict(
    boxstyle="round,pad=0.,rounding_size=0.2",fc='w',alpha=0.7,
    ec='none'
)
letter_text_props = dict(size='large', weight='bold')

import ckscool.cuts

def fig_cuts_kepmag_steff():
    cuts('kic_kepmag','m17_steff',nrows=2,ncols=4,plot_func=None,stars=True)
    fig = gcf() 
    axL = fig.get_axes()
    setp(axL,ylim=(3000,7000))
    setp(axL[0],xlabel='Kp (mag)', ylabel='Teff (K)')
    
def fig_cuts_period_prad():
    cuts('koi_period','koi_prad',nrows=2,ncols=4,plot_func=loglog)
    fig = gcf() 
    axL = fig.get_axes()
    setp(
        axL[0],xlabel='Orbital Period (days)',
        ylabel='Planet Size (Earth-radii)'
    )

def fig_cuts_smass_steff():
    cuts('m17_smass','m17_steff',nrows=2,ncols=4,plot_func=semilogx,stars=True)
    fig = gcf() 
    axL = fig.get_axes()
    setp(axL,ylim=(3000,7000))
    setp(
        axL[0],xlabel='Stellar Mass (Solar-masses)',
        ylabel='Teff (K)'
    )

def cuts(xk,yk,nrows=None,ncols=None,plot_func=None, stars=False):
    sns.set_style('ticks')
    sns.set_context('paper')

    if type(plot_func)!=type(None):
        plot = plot_func
    else:
        plot = plt.plot

        
    table = 'koi-mullally15'
    df = ckscool.io.load_table('ckscool-cuts')
    cuttypes = df.cuttypes

    width = ncols * 2
    height = nrows * 2

    fig, axL = subplots(
        nrows=nrows,ncols=ncols,figsize=(width,height),
        sharex=True,sharey=True
    )
    axL = axL.flatten()

    bpass = np.zeros(len(df))
    iplot = 0
    
    for cuttype in df.cuttypes:
        ax = axL[iplot] 
        sca(ax)
        
        key = 'is'+cuttype
        obj = ckscool.cuts.get_cut(cuttype)
        cut = obj(df,table)
        df[key] = cut.cut()
        bpass += df[key].astype(int)
            
        plotkw = dict(ms=4,rasterized=True)

        x = df[xk]
        y = df[yk]

        plot(x, y,'.',color='LightGray',**plotkw)
        dfcut = df[bpass==0]
        
        plot(dfcut[xk], dfcut[yk],'.', **plotkw)
        _text = cut.plotstr + ' ({})'.format(len(dfcut))

        if stars:
            _text = cut.plotstr + ' ({})'.format(len(dfcut.id_kic.drop_duplicates()))


        textkw = dict(fontsize='small',transform=ax.transAxes, ha='right')
        text(0.95,0.05, _text, **textkw)

        letter = string.ascii_lowercase[iplot]
        at = AnchoredText(letter,loc=2, frameon=True, prop=letter_text_props)
        ax.add_artist(at)
        setp(at.patch,**letter_bbox_props)
        
        #minorticks_off()
        iplot+=1

    axL = axL.reshape(nrows,ncols)
    #setp(axL[-1,:],xlabel='$R_p\, (R_{\oplus})$')
    #setp(axL[:,0],ylabel='[Fe/H]')
    fig.set_tight_layout(True)
    #fig.subplots_adjust(hspace=0.001,wspace=0.001,left=0.07,right=0.99,top=0.95,bottom=0.12)


class ComparisonPlotter(object):
    def __init__(self):

        kw = dict(lw = 2,histtype='step')
        self.cks1kw = dict(hatch='//',color='g',label='CKS-I',**kw)
        self.cksckw = dict(hatch='\\',color='b',label='CKS-Cool',**kw)

        df = ckscool.io.load_table('ckscool-cuts')
        df = df[~df.isany]
        df = df.groupby('id_kic',as_index=False).first()

        self.cksc = df

        df = ckscool.io.load_table('j17+m17')
        df = df[df.m17_kepmag < 14.2]
        df = df.groupby('id_kic',as_index=False).first()
        self.cks1 = df
    
    def _provision_figure(self):
        figsize = (2.5,2.5)
        fig = figure(figsize=figsize)
        return fig 

    def hist_kepmag(self):
        #fig = self._provision_figure()
        bins = arange(10,16.001,0.2)
        hist(self.cks1.m17_kepmag,bins=bins,**self.cks1kw)
        hist(self.cksc.m17_kepmag,bins=bins,**self.cksckw)

        xlabel('Kepmag')
        ylabel('Stars per bin')
        legend(loc='upper left',fontsize='small')
        #fig.set_tight_layout(True)
        #fig.savefig('proposal/fig_kepmag-hist.pdf',transparent=True)

    def hist_smass(self):
        #fig = self._provision_figure()
        bw = 0.05
        bins = arange(0.,1.5,bw)
        bins = np.logspace(log10(.2),log10(2.0),30)

        semilogx()
        hist(self.cks1.m17_smass.dropna(),bins=bins,**self.cks1kw)
        hist(self.cksc.m17_smass.dropna(),bins=bins,**self.cksckw)
        xt = [0.3,0.5,0.7,1.0,1.5,2.0]
        xticks(xt,xt)
        minorticks_off()
        legend(loc='upper left')
        xlabel('Mass (Solar Masses)')
        ylabel('Number of stars per bin')
        xlim(0.25,2)
        #fig.set_tight_layout(True)

    def hist_steff(self):
        #fig = self._provision_figure()
        bw = 200
        bins = arange(3000,7000,bw)

        hist(self.cks1.m17_steff.dropna(),bins=bins,**self.cks1kw)
        hist(self.cksc.m17_steff.dropna(),bins=bins,**self.cksckw)
        legend(loc='upper left',fontsize='medium')
        xlabel('Teff')
        ylabel('Number of Stars per bin')
        #fig.set_tight_layout(True)
        #fig.savefig('proposal/fig_teff-hist.pdf',transparent=True)

def fig_compare_with_cks1():
    """Compare CKS-Cool with CKS-I
    """
    fig, axL = subplots(nrows=1,ncols=3,figsize=(9.5,3))
    cp = ckscool.plot.sample.ComparisonPlotter()
    sca(axL[0])
    cp.hist_kepmag()
    sca(axL[1])
    cp.hist_smass()
    sca(axL[2])
    cp.hist_steff()
    fig.set_tight_layout(True)

