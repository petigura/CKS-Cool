import string
from collections import Iterable

import pandas as pd
import seaborn as sns
from matplotlib.pylab import *
from matplotlib.ticker import FuncFormatter, MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import cksphys.io
import cksmet.io
import cksmet.cuts
from cksmet.plotting.samples import add_anchored
from cksmet.plotting import *
from cksmet.plotting.config import *

texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
perticks = [0.3,1,3,10,30,100,300]

letter_bbox_props = dict(
    boxstyle="round,pad=0.,rounding_size=0.2",fc='w',alpha=0.7,
    ec='none'
)
letter_text_props = dict(size='large', weight='bold')

def prad_fe_label():
    # median planet radius error
    xlabel('Planet size (Earth-radii)')
    ylabel('[Fe/H]')
    xticks(rpticks,rpticks)
    xlim(0.5,25)
    ylim(-0.6,0.5)

def per_fe_label():
    # median planet radius error
    xlabel('Orbital Period (days)')
    ylabel('[Fe/H]')
    xticks(perticks,perticks)
    xlim(0.3,300)
    ylim(-0.6,0.5)

def prad_fe_errbar():
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    ebarx = 8.
    ebary = -0.51
    xferr = (cks.iso_prad_err1 / cks.iso_prad).median()    
    xerr = [[ ebarx - ebarx/(1+xferr) ], [ebarx*(1+xferr) - ebarx]]
    yerr = [[0.04],[0.04]]
    s='Median   \nUncert.   '
    errorbar([ebarx], [ebary], yerr=yerr, xerr=xerr,fmt='o',zorder=10,color='g')
    errorbar([ebarx], [ebary], yerr=yerr, xerr=xerr,fmt='o',lw=0,ms=0,label=s,color='g')
    legend(markerscale=0,frameon=True,scatterpoints=0,loc='lower right')

def per_fe_errbar():
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    ebarx = 0.50
    ebary = -0.51
    xferr = 1e-5
    xerr = [[ ebarx - ebarx/(1+xferr) ], [ebarx*(1+xferr) - ebarx]]
    yerr = [[0.04],[0.04]]
    s='Median   \nUncert.   '
    errorbar([ebarx], [ebary], yerr=yerr, xerr=xerr,fmt='o',zorder=10,color='g')
    errorbar([ebarx], [ebary], yerr=yerr, xerr=xerr,fmt='o',lw=0,ms=0,label=s,color='g')
    legend(markerscale=0,frameon=True,scatterpoints=0,loc='lower left')

def bins_to_xerr(bin0,binc,bin1):
    xerr = [[np.array(binc - bin0)], [np.array(bin1 - binc)] ]
    xerr = np.vstack(xerr)
    return xerr

def _plot_prad_fe(cks, **kw):
    x = np.array(cks.iso_prad)
    y = np.array(cks.cks_smet)
    plot(x,y,'.', **kw)

def _plot_per_fe(cks, **kw):
    x = np.array(cks.koi_period)
    y = np.array(cks.cks_smet)
    plot(x,y,'.', **kw)
    
def prad_fe():
    """Plot of planet radius vs stellar metallicty"""
    cks = cksmet.io.load_table('cks-cuts', cache=1)
    yk = 'cks_smet'
    cks = cks[~cks.isany]
    bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 8.0, 16]
    cksbin = table_bin(cks, bins, yk)

    #figure(figsize=(6,4))
    semilogx()
    _plot_prad_fe(cks)
    prad_fe_label()
    prad_fe_errbar()
    yk = 'cks_smet_mean'
    yerrk = 'cks_smet_mean_err'
    gcf().set_tight_layout(True)

def prad_smet_percentiles():
    """Plot of planet radius vs stellar metallicty"""
    cks = cksmet.io.load_table('cks-cuts', cache=1)
    cks = cks[~cks.isany]
    bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 8.0, 16]
    xk = 'binc'
    yk = 'cks_smet'
    yerrk = 'cks_smet_mean_err'
    binkey = 'iso_prad'
    cksbin = table_bin(cks, bins, yk, binkey)
    #figure(figsize=(6,4))
    semilogx()
    _plot_prad_fe(cks,color='LightGray')
    prad_fe_label()
    prad_fe_errbar()

    i = 0 
    for _i, row in cksbin.iterrows():
        x = row[xk]
        xerr = [[x - row.bin0], [row.bin1 - x] ]
        yerr = row[yerrk]
        yk = 'cks_smet_mean'
        y = row[yk]
        errorbar(x, row[yk], xerr=xerr, yerr=yerr, color='r',zorder=10)
        yk = 'cks_smet'
        for p in [25,75]:
            pk = yk+'_{:02d}'.format(p)
            errorbar(x,row[pk],xerr=xerr,color='m',zorder=10)
        i+=1
    gcf().set_tight_layout(True)

def per_smet_percentiles():
    """Plot of planet radius vs stellar metallicty"""
    cks = cksmet.io.load_table('cks-cuts', cache=1)
    cks = cks[~cks.isany]
    bins = 10**np.arange(np.log10(1),np.log10(150)+1e-5,0.25)
    xk = 'binc'
    yk = 'cks_smet'
    yerrk = 'cks_smet_mean_err'
    bink= 'koi_period'
    cksbin = table_bin(cks, bins, yk, bink)
    semilogx()
    _plot_per_fe(cks,color='LightGray')
    per_fe_label()
    per_fe_errbar()

    i = 0 
    for _i, row in cksbin.iterrows():
        x = row[xk]
        xerr = [[x - row.bin0], [row.bin1 - x] ]
        yerr = row[yerrk]
        yk = 'cks_smet_mean'
        y = row[yk]
        errorbar(x, row[yk], xerr=xerr, yerr=yerr, color='r',zorder=10)
        yk = 'cks_smet'
        for p in [25,75]:
            pk = yk+'_{:02d}'.format(p)
            errorbar(x,row[pk],xerr=xerr,color='m',zorder=10)
        i+=1
    gcf().set_tight_layout(True)


def fig_percentiles():
    sns_set_style('ticks')
    fig,axL = subplots(nrows=1, ncols=2, sharey=True,figsize=(7.25,3.5))
    sca(axL[0])
    prad_smet_percentiles()
    fig_label("a")
    sca(axL[1])
    per_smet_percentiles()
    setp(axL[1],ylabel='')
    fig_label("b")
    tight_layout(True)

def table_bin(df, bins, key, binkey):
    g = df.groupby(pd.cut(df[binkey],bins=bins))
    g = g[key]
    dfbin = pd.DataFrame(index=g.first().index)
    dfbin['bin0'] = bins[0:-1]
    dfbin['bin1'] = bins[1:]
    dfbin['binc'] = dfbin.eval('sqrt(bin0*bin1)')
    dfbin['count'] = g.count()
    dfbin[key + '_mean'] = g.mean()
    dfbin[key + '_std'] = g.std()
    dfbin[key +'_mean_err'] = dfbin[key + '_std']/ np.sqrt(dfbin['count'])

    for p in [1,5,25,50,75,95,99]:
        k = key+'_{:02d}'.format(p)
        dfbin[k] = g.quantile(q=p * 0.01)
    return dfbin

def prad_fe_fit():
    cksbin = cksmet.io.load_table('cksbin-fe')
    x = log10(cksbin.binc)
    y = cksbin.fe_mean
    yerr = cksbin.fe_mean_err
    import lmfit

    p0 = lmfit.Parameters()
    p0.add('p1',0.1)
    p0.add('p0',0.1)

    def model_linear(params,x):
        return params['p1'] * x + params['p0']

    def model_step(params, x):
        if isinstance(x,Iterable):
            out = map(lambda x : model_step(params,x), x)
            out = np.array(out)
            return out

        if x < log10(4):
            return params['p0']
        elif x >=  log10(4):
            return params['p1']

    def resid(params, model, x, y, yerr):
        return (y - model(params,x))/yerr

    res_linear = lmfit.minimize(resid,p0,args=(model_linear, x,y,yerr,))
    p1_linear = res_linear.params
    print "linear model"
    lmfit.report_fit(res_linear)

    res_step = lmfit.minimize(resid,p0,args=(model_step, x,y,yerr,))
    p1_step = res_linear.params
    print "step model"
    lmfit.report_fit(res_step)


    print "BIC(linear) - BIC(step) = {}".format(res_linear.bic - res_step.bic)

    semilogx()
    xerr = bins_to_xerr(cksbin.bin0,cksbin.binc,cksbin.bin1)
    print xerr
    errorbar(10**x,y,xerr=xerr, yerr=yerr,fmt='o')
    plot(10**x,model_linear(p1_linear,x))

    prad_fe_label()
    xlim(0.5,16)


def cuts():
    sns_set_style('ticks')
    df =  cksmet.io.load_table('cks-cuts',cache=1)

    cuttypes = cksmet.cuts.plnt_cuttypes
    nrows = len(cuttypes)

    nrows = 2
    ncols = 4

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
        obj = cksmet.cuts.get_cut(cuttype)
        cut = obj(df,'cks')
        df[key] = cut.cut()
        bpass += df[key].astype(int)
        plotkw = dict(ms=4)
        plot(df.iso_prad, df.cks_smet,'.',color='LightGray',**plotkw)
        dfcut = df[bpass==0]
        plot(dfcut.iso_prad, dfcut.cks_smet,'.',**plotkw)

        
        _text = cut.plotstr + ' ({})'.format(len(dfcut))
        textkw = dict(fontsize='small',transform=ax.transAxes, ha='right')
        text(0.95,0.05, _text, **textkw)

        semilogx()
        xlim(0.3333,48) # Spaceing between plots
        xt = [0.5,1,2,4,8, 16, 32]
        xticks(xt,xt)

        letter = string.ascii_lowercase[iplot]
        at = AnchoredText(letter,loc=2, frameon=True, prop=letter_text_props)
        ax.add_artist(at)
        setp(at.patch,**letter_bbox_props)
        
        minorticks_off()
        iplot+=1

    axL = axL.reshape(nrows,ncols)
    setp(axL[-1,:],xlabel='$R_p\, (R_{\oplus})$')
    setp(axL[:,0],ylabel='[Fe/H]')
    #tight_layout()

    fig.subplots_adjust(hspace=0.001,wspace=0.001,left=0.07,right=0.99,top=0.95,bottom=0.12)

import ckscool.io
def period_prad_slices(mode='tall'):
    sns_set_style('ticks')
    sns.set_color_codes()

    yt = [0.5, 1, 2, 4, 8, 16, 32]
    xt = [0.3, 1, 3, 10, 30, 100, 300]

    # Provision figure
    height = 3.375
    width = 3.6
    if mode=='four-equal-smet':
        smet_bins = [-0.75, -0.45, -0.15, 0.15, 0.45]
        
        labels = []
        for i in range(4):
            labels+='[Fe/H] = ({:+.2f},{:+.2f})'.format(
                smet_bins[i],smet_bins[i+1]
            )
        
        nrows = 2
        ncols = 2
        fig,axL = subplots(
            nrows=nrows,ncols=ncols,figsize=(ncols*width,nrows*height),
            sharex=True,sharey=True
        )

    if mode=='four-equal-stars':

        lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1,cachefn='../CKS-Metallicity/load_table_cache.hdf')
        lamo = lamo[~lamo.isany]
        quantiles = [0.25, 0.5, 0.75]
        lamo = lamo.lamo_smet
        lamoq = lamo.quantile(quantiles)
        smet_bins = [-10] + list(lamoq) + [10]
        labels = [
            '[Fe/H] < ${:.3f}$'.format(smet_bins[1]),
            '[Fe/H] = (${:+.3f}$,${:+.3f}$)'.format(*smet_bins[1:3]),
            '[Fe/H] = (${:+.3f}$,${:+.3f}$)'.format(*smet_bins[2:4]),
            '[Fe/H] > ${:+.3f}$'.format(smet_bins[3]),
        ]
        nrows = 2
        ncols = 2
        fig,axL = subplots(
            nrows=nrows,ncols=ncols,figsize=(ncols*width,nrows*height),
            sharex=True,sharey=True
        )

        setp(axL[1,:], xlabel='Orbital Period (days)')
        setp(axL[:,0], ylabel='Planet Size (Earth-radii)')

    # Load up LAMOST metallicities
    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]

    plnt = ckscool.io.load_table('ckscool-planets-cuts')
    plnt = plnt[~plnt.isany]
    i = 0
    nbins = len(smet_bins) - 1
    axL = axL.flatten()

    while True:
        if i==nbins:
            break
        ax = axL[i]
        sca(ax)
        loglog()
        smet1 = smet_bins[i]
        smet2 = smet_bins[i+1]
        plnt_cut = plnt[plnt.sm_smet.between(smet1,smet2)]
        lamo_cut = lamo[lamo.lamo_smet.between(smet1,smet2)]
        
        f_planet = 1.0 * len(plnt_cut) / len(plnt)        
        f_stars = 1.0 * len(lamo_cut) / len(lamo)

        label = labels[i]
        label += "\n$f_p = {:.0f}\%$".format(100*f_planet)

        # Whole sample for comparison
        kwpts = dict(marker='o',ms=4, lw=0)
        kw = dict(label=label,color='LightGray', **kwpts)
        plot(plnt.koi_period,plnt.gdir_prad, **kw)
        
        # Whole sample for comparison
        kw = dict(label='',color='b', **kwpts)
        plot(plnt_cut.koi_period,plnt_cut.gdir_prad, **kw)

        #legend(frameon=True,markerscale=0,framealpha=0.5,fontsize='small',handletextpad=-2,loc='upper right')

        bbox_props = dict(
            boxstyle="round,pad=0.,rounding_size=0.2",fc='w',alpha=0.7,
            ec='none'
        )

        at = AnchoredText(label,loc=1, frameon=True,prop=dict(size='small'))
        ax.add_artist(at)
        setp(at.patch,**bbox_props)

        yticks(yt,yt)
        xticks(xt,xt)
        i+=1

    # Errorbar 
    for ax, letter in zip(axL.flatten(),string.ascii_lowercase):
        sca(ax)
        at = AnchoredText(letter,loc=2, frameon=True, prop=letter_text_props)
        ax.add_artist(at)
        setp(at.patch,**letter_bbox_props)
        grid()


    # Errorbar 
    sca(axL[0])
    ebarx = sqrt(3) * 100.
    ebary = 0.7
    yferr = 0.1
    yerr = [[ ebary - ebary/(1+yferr) ], [ebary*(1+yferr) - ebary]]
    xerr = [[0.00],[0.00]]
    errorbar([ebarx], [ebary], yerr=yerr, zorder=10, fmt='o',c='g',ms=4)
    xlim(0.3,1000)
    ylim(0.5,32)
    fig.set_tight_layout(True)
    text(ebarx,ebary,'  Median\n  Uncert.',size='small',va='center')
    minorticks_off()
