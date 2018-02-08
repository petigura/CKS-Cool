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

def cuts(xk,yk,nrows=None,ncols=None,plot_func=None):
    sns.set_style('ticks')
    #df =  cksmet.io.load_table('cks-cuts',cache=1)

    if type(plot_func)!=type(None):
        plot = plot_func
    else:
        plot = plt.plot

    table = 'koi-thompson18'
    cuttypes = ['none','notreliable','badteffphot','faint',]
    df = ckscool.io.load_table(table)
    df = ckscool.cuts.add_cuts(df,cuttypes,table)

    #cuttypes = cksmet.cuts.plnt_cuttypes
    width = ncols * 2
    height = nrows * 2

    fig, axL = subplots(
        nrows=nrows,ncols=ncols,figsize=(width,height),
        sharex=True,sharey=True
    )
    axL = axL.flatten()

    bpass = np.zeros(len(df))
    iplot = 0
    
    for cuttype in cuttypes:
        ax = axL[iplot] 
        sca(ax)
        
        key = 'is'+cuttype
        obj = ckscool.cuts.get_cut(cuttype)
        cut = obj(df,table)
        df[key] = cut.cut()
        bpass += df[key].astype(int)
        plotkw = dict(ms=4)
        plot(df[xk], df[yk],'.',color='LightGray',**plotkw)
        dfcut = df[bpass==0]
        
        plot(dfcut[xk], dfcut[yk],'.', **plotkw)
        _text = cut.plotstr + ' ({})'.format(len(dfcut))
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
    #tight_layout()

    fig.subplots_adjust(hspace=0.001,wspace=0.001,left=0.07,right=0.99,top=0.95,bottom=0.12)
    


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

def lamo_detectability():
    querys = ['lamo_smet < 0','0 < lamo_smet']
    colors = ['b','r']
    labels = ['[Fe/H] < 0','[Fe/H] > 0']

    lamo = cksmet.io.load_table('lamost-cal-cuts+cdpp',cache=1)
    lamo = lamo[~lamo.isany]
    fig1 = figure(figsize=(4,4))
    semilogy()

    bins = arange(11,14.1,0.5)
    #querys = ['-0.5 < lamo_smet < -0.1','0.1 < lamo_smet < 0.5']
    for i in range(len(querys)):
        query = querys[i]
        color=colors[i]
        label=labels[i]
        cut = lamo.query(query) 
        plot(
            cut.kepmag,cut.cdpp3,'.',color=color,zorder=1,label=label,alpha=0.8
        )
        g = cut.groupby(pd.cut(cut.kepmag,bins))
        i = 0

        for interval, cdpp in g['cdpp3'].median().iteritems():
            if i==0:
                label=query
            else:
                label=None
            plot([interval.left,interval.right], [cdpp,cdpp],color='w',lw=6)
            plot([interval.left,interval.right], [cdpp,cdpp],color=color,lw=3)
            i+=1
            print query, interval.left,interval.right, cdpp

    xlim(11,14)
    ylim(15,150)
    legend()

    xlabel('kepmag')
    ylabel('CDPP3 (ppm)')
    minorticks_off()
    yticks([20,30,40,50,60,70,80,90,100],[20,30,40,50,60,70,80,90,100])
    fig1.set_tight_layout(True)

    fig2 = figure(figsize=(4,4))
    kepmag = (12,14)
    cut2 = lamo[lamo.kepmag.between(*kepmag)]
    fsamp =  100 * len(cut2) / len(lamo)
    print "for kepmag = {}, {:.2f}\% of stellar sample median CDPP3 is ".format(kepmag, fsamp)

    fig2 = figure(figsize=(4,4))
    for i in range(2):
        query = querys[i]
        color=colors[i]
        label=labels[i]
        cut = lamo.query(query) 
        plot(cut.lamo_steff,cut.lamo_slogg,'.',color=color)

    xlim(6500,4700)
    ylim(5.0,4.0)
    xlabel('Effective Temp. (K)')
    ylabel('logg (dex)')
    fig2.set_tight_layout(True)
    return fig1,fig2

def smet_snr():
    sns_set_style('ticks')
    fig,axL = subplots(figsize=(4,4))
    df = pd.read_csv('isoclassify-lamost-dr2.csv')
    huber14 = cksmet.io.load_table('huber14+cdpp',cache=1)
    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    lamo = pd.merge(lamo,huber14['id_kic kepmag cdpp3'.split()],on='id_kic')
    df = pd.merge(lamo,df,left_on='id_kic',right_on='id_starname')

    semilogy()
    srad = df.iso_srad * c.R_sun
    prad = 1 * c.R_earth
    df['snr'] = (prad / srad)**2 / (1e-6 * df.cdpp3)
    plot(df.lamo_smet,df.snr,',')

    bins = arange(-0.3,0.31,0.15,)
    g = df.groupby(pd.cut(df.lamo_smet,bins))
    i = 0
    for interval, snr in g['snr'].median().iteritems():
        x = [interval.left,interval.right]
        y = [snr,snr]
        plot(x,y,color='w',lw=6)
        plot(x,y,color='b',lw=3)
        i+=1
        print interval.left,interval.right, snr

    xlim(-0.4,0.4)
    ylim(0.1,10)
    yt = [0.1,0.3,1,3,10]
    yticks(yt,yt)
    xlabel('[Fe/H] (dex)')
    ylabel('Single Transit SNR')
    fig.set_tight_layout(True)
    return fig
