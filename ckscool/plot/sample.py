import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy import constants as c
from ckscool.plot.config import *

errorbar_kw = dict(fmt='.',markersize=5,color='b')

#sns.set_style('whitegrid')
sns.set_color_codes()

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
texteff = '$\mathregular{T}_{\mathregular{eff}}$'
texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 
texrs = '$R_\star (R_\odot)$' 

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
perticks = [0.3,1,3,10,30,100,300]

letter_bbox_props = dict(
    boxstyle="round,pad=0.,rounding_size=0.2",fc='w',alpha=0.7,
    ec='none'
)
letter_text_props = dict(size='large', weight='bold')

import ckscool.cuts

#def fig_cut_field_steff_srad():
evalMG = 'gaia2_gmag - 5 * log10(bj_dist/10)'
evalBR ='gaia2_bpmag - gaia2_rpmag' 

def fig_cuts_all_multi():
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.0)
    sns.set_color_codes()

    fig,axL = subplots(nrows=2,ncols=2,figsize=(7,5))
    sca(axL[0,0])
    fig_cuts_all('field-cuts','m17_steff','m17_kepmag')
    sca(axL[0,1])
    fig_cuts_all('planets-cuts1','m17_steff','m17_kepmag')

    sca(axL[1,0])
    fig_cuts_all('field-cuts','m17_steff','gaia2_srad')
    semilogy()
    sca(axL[1,1])
    fig_cuts_all('planets-cuts1','m17_steff','gaia2_srad')
    semilogy()
    tight_layout()


def fig_cuts_all_multi2():
    """
    Same as multi1, but plotting 
    """
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.0)
    sns.set_color_codes()


    fig,axL = subplots(nrows=2,ncols=2,figsize=(7,5))
    sca(axL[0,0])
    fig_cuts_all('field-cuts',evalBR,'m17_kepmag')
    sca(axL[0,1])
    fig_cuts_all('planets-cuts1',evalBR,'m17_kepmag')

    sca(axL[1,0])
    fig_cuts_all('field-cuts',evalBR,evalMG)
    sca(axL[1,1])
    fig_cuts_all('planets-cuts1',evalBR,evalMG)
    tight_layout()


def fig_cuts_all_multi3():
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.0)
    #sns.set_color_codes()
    sns.color_palette("bright")

    fig,axL = subplots(nrows=2,ncols=2,figsize=(7,5))
    starkw = dict(marker=',',lw=0,rasterized=True,color='DarkGray')
    starselkw = dict(marker=',',lw=0,rasterized=True,color='r')

    planetkw = dict(marker='.',lw=0,ms=4,mew=0,color='DarkGray')
    planetselkw = dict(marker='.',lw=0,ms=4,mew=0,color='r')
    
    sca(axL[0,0])
    xk = 'm17_steff'
    yk = 'm17_kepmag'
    df = ckscool.io.load_table('field-cuts') 
    plot(df.eval(xk),df.eval(yk),**starkw)
    df = df[~df.isany]
    plot(df.eval(xk),df.eval(yk),**starselkw)

    sca(axL[0,1])
    df = ckscool.io.load_table('planets-cuts1')
    plot(df.eval(xk),df.eval(yk),**planetkw)
    df = df[~df.isany]
    plot(df.eval(xk),df.eval(yk),**planetselkw)

    sca(axL[1,0])
    cut = ckscool.cuts.occur.CutGiantCMD(df,'field')
    xk = 'gaia2_bpmag - gaia2_rpmag' 
    yk = 'gaia2_gmag - 5 * log10(bj_dist/10)'
    df = ckscool.io.load_table('field-cuts') 
    plot(df.eval(xk),df.eval(yk),**starkw)
    df = df[~df.isany]
    plot(df.eval(xk),df.eval(yk),**starselkw)
    x = linspace(0.5,2.5,100)
    
    sca(axL[1,1])
    df = ckscool.io.load_table('planets-cuts1')
    plot(df.eval(xk),df.eval(yk),**planetkw)
    df = df[~df.isany]
    plot(df.eval(xk),df.eval(yk),**planetselkw)
    x = linspace(0.5,2.5,100)

    for i in range(2):
        sca(axL[1,i])
        plot(x,cut.fabove(x),color='r',ls='--',lw=0.5)
        plot(x,cut.fbelow(x),color='r',ls='--',lw=0.5)

    setp(axL[0,:],ylim=(17,8),xlim=(7000,3000),xlabel=texteff,ylabel='$Kp$ (mag)')
    setp(axL[1,:],ylim=(13,-2),xlim=(0,3),xlabel='$B_p - R_p$ (mag)',ylabel='$M_G$ (mag)')
    tight_layout()


def fig_cuts_all(sample,xk,yk):
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.1)
    sns.set_color_codes()

    df = ckscool.io.load_table(sample,cache=2)
    x = df.eval(xk)
    y = df.eval(yk)
    plotkw = dict(ms=3,rasterized=True,alpha=0.4)
    plot(x, y,'.',color='LightGray', **plotkw)
    dfcut = df[~df.isany]
    xcut = dfcut.eval(xk)
    ycut = dfcut.eval(yk)
    plot(xcut, ycut,'.', color='b', **plotkw)
    ax = gca()

    xlim=None
    ylim=None
    xlabel=None
    ylabel=None
    if xk.count('steff'):
        xlim=(7000,3000)
        xlabel=texteff
    if xk==evalBR:
        xlabel='$B_p - R_p$ (mag)'
        xlim=(0,3)
    
    if yk.count('srad'):
        ylabel=texrs
        ylim=(0.1,10)        

    if yk=='m17_kepmag':
        ylabel='$Kp$ (mag)'
        ylim=(17,8)        

    if yk==evalMG:
        ylabel='$M_G$ (mag)'
        ylim=(13,-2)


    setp(ax,ylim=ylim,xlim=xlim,xlabel=xlabel,ylabel=ylabel)
    
    

def fig_cuts_all_field_steff_srad():
    df = ckscool.io.load_table('field-cuts',cache=2)
    xk = 'ber19_steff'
    yk = 'ber19_srad'


def fig_cuts_all_plnt_steff_srad():
    df = ckscool.io.load_table('planet-cuts1',cache=2)
    xk = 'ber19_steff'
    yk = 'ber19_srad'
    x = df.eval(xk)
    y = df.eval(yk)
    plotkw = dict(ms=2,rasterized=True,alpha=0.1)
    plot(x, y,'.',color='LightGray',**plotkw)
    dfcut = df[~df.isany]
    xcut = dfcut.eval(xk)
    ycut = dfcut.eval(yk)
    semilogy()
    plot(xcut, ycut,'.', color='blue',**plotkw)
    ax = gca()
    setp(ax,ylim=(0.1,10),xlim=(8000,3000),xlabel=texteff,ylabel=texrs)



def fig_cuts_field_bmv_gmag():
    nrows=2
    ncols=4
    df = ckscool.io.load_table('field-cuts',cache=2)

    df['dmod'] = 5 * log10(df.bj_dist/10)
    cuts(df,'gaia2_bpmag - gaia2_rpmag','gaia2_gmag - dmod',nrows=nrows,ncols=ncols,stars=False)
    axL = gcf().get_axes()
    axL = np.array(axL).reshape(nrows,ncols) 
    setp(axL,ylim=(15,-5),xlim=(-1,5))
    #setp(axL,xlim=(10,0))
    setp(axL[1,:],xlabel='$Bp - Rp$ (mag)')
    setp(axL[:,0],ylabel='$M_G$ (mag)')

def fig_cuts_field_steff_srad():
    nrows=2
    ncols=4
    df = ckscool.io.load_table('field-cuts',cache=2)
    cuts(df,'ber19_steff','ber19_srad',nrows=nrows,ncols=ncols,stars=False)
    axL = gcf().get_axes()
    axL = np.array(axL).reshape(nrows,ncols) 
    setp(axL,ylim=(0.1,10),xlim=(8000,3000))
    #setp(axL,xlim=(10,0))
    semilogy()
    setp(axL[1,:],xlabel=texteff)
    setp(axL[:,0],ylabel=texrs)

def fig_cuts_kepmag_steff():
    nrows=2
    ncols=4
    df = ckscool.io.load_table('planets-cuts1',cache=2)
    cuts(df,'m17_kepmag','m17_steff',nrows=nrows,ncols=ncols,plot_func=None,stars=False)
    axL = gcf().get_axes()
    axL = np.array(axL).reshape(nrows,ncols) 
    setp(axL,ylim=(3000,7000))
    setp(axL[1,:],xlabel='Kp (mag)')
    setp(axL[:,0],ylabel='Teff (K)')
    
def fig_cuts_period_prad():
    nrows=2
    ncols=4
    df = ckscool.io.load_table('planets-cuts1')
    cuts(df, 'koi_period','koi_prad',nrows=nrows,ncols=ncols,plot_func=loglog)
    axL = gcf().get_axes()
    axL = np.array(axL).reshape(nrows,ncols) 
    setp(axL[1,:],xlabel='Orbital Period (days)',xlim=(0.1,1000),ylim=(0.1,100))
    setp(axL[:,0],ylabel='Planet Size (Earth-radii)')

def fig_cuts_smass_steff():
    df = ckscool.io.load_table('planets-cuts1')
    cuts(df, 'm17_smass','m17_steff',nrows=2,ncols=4,plot_func=semilogx)
    axL = gcf().get_axes()
    setp(axL,ylim=(3000,7000))
    setp(
        axL[0],xlabel='Stellar Mass (Solar-masses)',
        ylabel='Teff (K)'
    )

def fig_cuts_stars_hr():
    df = ckscool.io.load_table('planets-cuts2',cache=2)
    nplots = len(df.cuttypes)
    cuts(
        df, 'cks_steff', 'gdir_srad', nrows=1, ncols=nplots, plot_func=semilogy
    )
    fig = gcf()
    axL = fig.get_axes()
    setp(axL,xlabel="{} (K)".format(texteff),xlim=(7000,3500),ylim=(0.4,2))
    setp(axL[0],ylabel='Stellar radius (Solar-radii)')
    for ax in axL:
        sca(ax)
        yt = [0.5,0.7,1.0,1.4,2.0]
        yticks(yt,yt)
        minorticks_off()
        ylim(0.4,3)


def fig_cuts_planets_per_prad(zoom=False):
    xk = 'koi_period'
    yk = 'gdir_prad'
    nrows = 2
    ncols = 4
    df = ckscool.io.load_table('planets-cuts2',cache=2)
    cuts(df, xk, yk , nrows=nrows, ncols=ncols, stars=False, plot_func=loglog)
    axL = gcf().get_axes()
    axL = np.array(axL).reshape(nrows,ncols) 
    yt = [0.5,1,2,4,8,16]
    xt = [0.3,1,3,10,30,100,300]
    for ax in axL.flatten():
        minorticks_off()
        sca(ax)
        xticks(xt,xt)
        yticks(yt,yt)
        grid()
    

    if zoom:
        ylim=(1,4)
    else:
        ylim=(0.5,20)
    setp(axL[1,:],xlabel='Orbital Period (days)',xlim=(0.3,300),ylim=ylim)
    setp(axL[:,0],ylabel='Planet Size (Earth-radii)')

def cuts(df, xk,yk,nrows=None,ncols=None,plot_func=None, stars=False):
    sns.set_style('ticks')
    sns.set_context('paper')

    if type(plot_func)!=type(None):
        plot = plot_func
    else:
        plot = plt.plot

    table = 'koi-mullally15'
    cuttypes = df.cuttypes

    width = ncols * 2
    height = nrows * 2.25

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
        obj = ckscool.cuts.occur.get_cut(cuttype)
        cut = obj(df,table)
        #df[key] = cut.cut()
        bpass += df[key].astype(int)
            
        plotkw = dict(ms=2,rasterized=True,alpha=0.1)

        x = df.eval(xk)
        y = df.eval(yk)

        plot(x, y,'.',color='LightGray',**plotkw)
        dfcut = df[bpass==0]

        xcut = dfcut.eval(xk)
        ycut = dfcut.eval(yk)

        plot(xcut, ycut,'.', color='blue',**plotkw)
        _text = cut.plotstr + ' ({})'.format(len(dfcut))

        if stars:
            _text = cut.plotstr + ' ({})'.format(len(dfcut.id_kic.drop_duplicates()))

        textkw = dict(fontsize='small',transform=ax.transAxes, ha='right')
        text(0.95,0.05, _text, **textkw)

        letter = string.ascii_lowercase[iplot]
        at = AnchoredText(letter,loc=2, frameon=True, prop={'weight':'bold','size':'large'})
        ax.add_artist(at)
        setp(at.patch,**letter_bbox_props)
        
        #minorticks_off()
        iplot+=1

    axL = axL.reshape(nrows,ncols)
    fig.set_tight_layout(True)



class ComparisonPlotter(object):
    def __init__(self):

        kw = dict(lw = 2,histtype='step')
        self.cks1kw = dict(hatch='//',color='g',label='CKS-I',**kw)
        self.cksckw = dict(hatch='\\',color='b',label='CKS-Cool',**kw)

        #df = ckscool.io.load_table('ckscool-cuts')
        df = ckscool.io.load_table('ckscool-stars')
        df = df[~df.isany]
        df = df.groupby('id_kic',as_index=False).first()

        self.cksc = df

        df = ckscool.io.load_table('j17+m17')
        df = df[df.kic_kepmag < 14.2]
        df = df.groupby('id_kic',as_index=False).first()
        self.cks1 = df
    
    def _provision_figure(self):
        figsize = (2.5,2.5)
        fig = figure(figsize=figsize)
        return fig 

    def hist_kepmag(self):
        #fig = self._provision_figure()
        bins = arange(10,16.001,0.2)
        hist(self.cks1.kic_kepmag,bins=bins,**self.cks1kw)
        hist(self.cksc.kic_kepmag,bins=bins,**self.cksckw)

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
        hist(self.cksc.giso_smass.dropna(),bins=bins,**self.cksckw)
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
        hist(self.cksc.sm_steff.dropna(),bins=bins,**self.cksckw)
        legend(loc='upper left',fontsize='medium')
        xlabel('Teff')
        ylabel('Number of Stars per bin')
        #fig.set_tight_layout(True)
        #fig.savefig('proposal/fig_teff-hist.pdf',transparent=True)

def fig_compare_with_cks1():
    """Compare CKS-Cool with CKS-I
    """
    fig, axL = subplots(nrows=1,ncols=3,figsize=(9.5,3))
    cp = ComparisonPlotter()
    sca(axL[0])
    cp.hist_kepmag()
    sca(axL[1])
    cp.hist_smass()
    sca(axL[2])
    cp.hist_steff()
    fig.set_tight_layout(True)


'''
def samples():
    ncols = 3
    nrows = 3
    height = 2.5 * nrows
    width = 3.0 * ncols

    fig, axL = subplots(nrows=nrows, ncols=ncols, figsize=(width, height))

    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamoc = lamo[~lamo.isany]

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cksc = cks.query('isany==False')
    
    cks = cks.groupby('id_kic').first()
    cksc = cksc.groupby('id_kic').first()
    
    field = cksmet.io.load_table('field-cuts',cache=1)
    fieldc = field.query('isany==False')

    kwpts = dict(color='LightGray',rasterized=True,marker='o',ms=1.5,lw=0)
    kwptsc = dict(color='b',rasterized=True,marker='o',ms=1.5,lw=0)

    bins = np.arange(9,16.0001,0.2)
    histkw = dict(bins=bins,histtype='stepfilled',color='LightGray')
    histkwc = dict(bins=bins,histtype='stepfilled',color='b')
    
    bins = np.arange(-1.0,0.5001,0.1)
    smethistkw = dict(bins=bins,histtype='stepfilled',color='LightGray')
    smethistkwc = dict(bins=bins,histtype='stepfilled',color='b')

    jcks = 0
    jfield = 1 
    jlamo = 2

    ihr = 0 
    ikepmag = 1 
    ismet = 2

    # HR diagrams
    sca(axL[0,jcks])
    plot(cks.cks_steff, cks.cks_slogg, **kwpts)
    plot(cksc.cks_steff, cksc.cks_slogg, **kwptsc)

    sca(axL[0,jfield])
    print list(field.columns)
    plot(field.m17_steff, field.m17_slogg, **kwpts)
    plot(fieldc.m17_steff, fieldc.m17_slogg, **kwptsc)

    sca(axL[0,2])
    plot(lamo.lamo_steff, lamo.lamo_slogg,**kwpts)
    plot(lamoc.lamo_steff, lamoc.lamo_slogg,**kwptsc)

    # KepMag
    sca(axL[1,jcks])
    hist(cks.m17_kepmag,**histkw)
    hist(cksc.m17_kepmag,**histkwc)

    sca(axL[1,jfield])
    hist(field.kepmag,**histkw)
    hist(fieldc.kepmag,**histkwc)

    sca(axL[1,2])
    hist(lamo.m17_kepmag,**histkw)
    hist(lamoc.m17_kepmag,**histkwc)
    
    # Metal
    sca(axL[ismet,jcks])
    hist(cks.cks_smet.dropna(),**smethistkw)
    hist(cksc.cks_smet.dropna(),**smethistkwc)

    sca(axL[ismet,jfield])
    hist(field.m17_smet.dropna(),**smethistkw)
    hist(fieldc.m17_smet.dropna(),**smethistkwc)

    sca(axL[ismet,jlamo])
    hist(lamo.lamo_smet,**smethistkw)
    hist(lamoc.lamo_smet,**smethistkwc)

    setp(
        axL[0,:], xlim=(7000,4000), ylim=(5,3.0), xlabel='Teff (K)', 
        ylabel='logg (dex)'
    )
    fig.set_tight_layout(True)

    setp(axL[1,:],ylabel='Number of Stars',xlabel='Kepmag')
    setp(axL[2,:],ylabel='Number of Stars',xlabel='[Fe/H]')

    for ax in axL[:,2]:
        sca(ax)
        add_anchored('LAMOST',loc=1)
        
    for ax in axL[:,jcks]:
        sca(ax)
        add_anchored('CKS',loc=1)

    for ax in axL[:,jfield]:
        sca(ax)
        add_anchored('Field',loc=1)

    for ax, letter in zip(axL.flatten(),string.ascii_lowercase):
        sca(ax)
        add_anchored(
            letter,loc=2, frameon=False, 
            prop=dict(size='large', weight='bold')
        )

'''
