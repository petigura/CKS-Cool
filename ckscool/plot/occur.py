import lmfit
from matplotlib.pylab import *
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
import pandas as pd
import seaborn as sns
from scipy import ndimage as nd

import ckscool.io
import ckscool.plot.planet

sns.set_style('ticks')
sns.set_color_codes()

def contour(cp, plot_interval=False,
            draw_colorbar=True,cax=None,plot_completeness=True,
            normalize=False, ntrials_min=50, levels=None):
    """
    Args:
       cp : contour plotter object

    """

    ax = gca()
    tax = gca().transAxes
    min_ntrials = 50 # minimum number of N trials to show
    cmap = 'YlGn' #,None #'hot_r'

    # convert into an x-array for plotting
    ds = cp.rate.groupby(['per1','prad1']).first().to_xarray()
    norm = cp.rate.query('10 < perc < 100 and 2 < pradc < 4').rate.sum()
    rate = ds.rate
    rate = rate.fillna(1e-10)
    rate = nd.gaussian_filter(rate,(4,2))

    # Smooth out the contours a bit
    eps = 1e-10
    X = np.array(np.log10(ds.perc))
    Y = np.array(np.log10(ds.pradc))
    Z = np.array(rate)
    if levels==None:
        #levels = arange(0,5e-2+eps,0.0025) 
        b = (
            (ds.ntrial > min_ntrials) 
            & (ds.perc > 1)
            & (ds.perc < 300)
            & (ds.pradc > 1.0) 
            & (ds.pradc < 4.0)
        )
        b = array(b)
        maxz = np.max(rate[b])
        maxz = np.round(maxz*1.1,3)
        levels = linspace(0,maxz+eps,20)

    cbarticks = levels[::2]
    cbarticklabels = ["{:.1f}".format(1e2*_yt) for _yt in cbarticks]
    kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)
    cbarlabel = r"""Planets per 100 Stars per $P-R_P$ interval"""
    qcs = contourf(X,Y,Z, **kw)

    # Completeness
    if plot_completeness:
        Z = np.array(ds.ntrial)
        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(X,Y,Z,[0,min_ntrials],zorder=2.5,cmap=cmap,vmax=1)
        '''
        text(
            0.95,0.15,'Low Completeness',rotation=12,zorder=5,size='small',
            transform=ax.transAxes,ha='right'
        )
        '''


    # plot straight contours
    #kw.pop('cmap')
    #kw.pop('extend')
    #qcs = contour(X,Y,Z,)
    #plt.clabel(qcs, inline=1, fmt='%.3f', colors='w', fontsize=1)
    if draw_colorbar:
        cbar = colorbar(qcs,cax=cax,ticks=cbarticks,)
        t = cbar.ax.set_yticklabels(cbarticklabels)
        setp(t,size='x-small')
        cbar.set_label(cbarlabel,size='small')


    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16,32]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    xlim(log10(1),log10(300))
    ylim(log10(1),log10(4))
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')

    if plot_interval:
        #axis_to_data = ax.transAxes + ax.transData.inverted()
        #xy = axis_to_data.transform((0.1,0.9))
        xy = 0.1, 0.52
        w = cp.perwid
        h = cp.pradwid
        rect = Rectangle(xy, w, h, lw=1, ec='r', fc='none', zorder=10)
        ax.add_patch(rect)
        s = """\
$P-R_P$ Interval
$\Delta \log P$ = {:.2f} dex
$\Delta \log R_P$ = {:.2f} dex
""".format(w,h)
        kw = dict(
            size='x-small',zorder=5,va='top',ha='left',transform=ax.transAxes,
        )
#        text(xyaxes[0]+0.07,xyaxes[1],s,**kw)
# ---------------------------------------------------------------------------- #

def contour_sinc(cp, plot_interval=False,
                 draw_colorbar=True,cax=None,plot_completeness=True,label=False,
                 normalize=False, ntrials_min=50):
    """
    Args:
       cp : contour plotter object

    """

    ax = gca()
    tax = gca().transAxes

    # convert into an x-array for plotting
    ds = cp.rate.groupby(['sinc1','prad1']).first().to_xarray()
    norm = cp.rate.query('10 < sincc < 1000 and 2 < pradc < 4').rate.sum()
    rate = ds.rate
    rate = rate.fillna(4e-6)

    # Smooth out the contours a bit
    rate = nd.gaussian_filter(rate,(4,2))

    eps = 1e-10
    X, Y = np.log10(ds.sincc), np.log10(ds.pradc)
    cmap = 'YlGn' #,None #'hot_r'
    levels = None
    cbarlabel=''

    if normalize:
        Z = rate / norm
        levels = linspace(0,Z.max(),20)
        kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)
        cbarticks = levels[::2]
        cbarticklabels = ["{:.0f}".format(1e2*_yt) for _yt in cbarticks]

    else:
        levels = arange(0,5e-2+eps,0.0025)
        Z = rate
        cbarticks = levels[::2]
        cbarticklabels = ["{:.0f}".format(1e2*_yt) for _yt in cbarticks]
        kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)

    cbarlabel = r"""Planets per 100 Stars per $Sinc-R_P$ interval"""

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    qcs = contourf(X,Y,Z, **kw)

    # plot straight contours
    if draw_colorbar:
        cbar = colorbar(qcs,cax=cax,ticks=cbarticks,)
        t = cbar.ax.set_yticklabels(cbarticklabels)
        setp(t,size='x-small')
        cbar.set_label(cbarlabel,size='small')

    # Completeness
    if plot_completeness:
        Z = np.array(ds.ntrial)
        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(X,Y,Z,[0,ntrials_min],zorder=2.5,cmap=cmap,vmax=1)

    if plot_interval:
        #inv = ax.transAxes.inverted()
        xyaxes = (0.1,0.9)
        xy = ax.transLimits.inverted().transform(xyaxes)
        #xy = ax.transAxes.transform((0.9,0.9))
        w = cp.sincwid
        h = cp.pradwid
        rect = Rectangle(xy, w, h,lw=1,ec='r',fc='none',zorder=4)
        ax.add_patch(rect)
        s = """\
$P-R_P$ Interval
$\Delta \log Sinc$ = {:.2f} dex
$\Delta \log R_P$ = {:.2f} dex
""".format(w,h)
        kw = dict(
            size='x-small',zorder=5,va='top',ha='left',transform=ax.transAxes,
        )
#        text(xyaxes[0]+0.07,xyaxes[1],s,**kw)

    xt = [1,3,10,30,100,1000,10000]
    yt = [0.5,1,2,4,8,16,32]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    xlim(log10(10000),log10(1))
    ylim(log10(1),log10(4))
    xlabel('Stellar Incident Flux (Earth Units)')
    ylabel('Planet Size (Earth-radii)')

def fig_contour_three():
    cp0 = ckscool.io.load_object('cp_smass=0.5-0.7',cache=1)
    cp1 = ckscool.io.load_object('cp_smass=0.7-1.0',cache=1)
    cp2 = ckscool.io.load_object('cp_smass=1.0-1.4',cache=1)
    fig, axL = subplots(ncols=3,figsize=(8,3))
    sca(axL[0])
    contour(cp0,plot_planetpoints=False,plot_interval=True,draw_colorbar=False)
    title('$M_\star = 0.5-0.7\, M_\odot$ ')

    sca(axL[1])
    contour(cp1,plot_planetpoints=False,plot_interval=True,draw_colorbar=False)
    title('$M_\star = 0.7-1.0\, M_\odot$ ')

    sca(axL[2])
    contour(cp2,plot_planetpoints=False,plot_interval=True,draw_colorbar=False)
    title('$M_\star = 1.0-1.4\, M_\odot$ ')

    xlim=log10(1),log10(300)
    ylim=log10(1),log10(4)
    setp(axL,xlim=xlim,ylim=ylim)
    fig.subplots_adjust(hspace=0.2)

def gradient_arrays(cp):

    """
    Args:
       cp : contour plotter object

    """
    # convert into an x-array for plotting
    ds = cp.rate.groupby(['per1','prad1']).first().to_xarray()
    rate = ds.rate
    rate = rate.fillna(4e-6)

    # Smooth out the contours a bit
    #rate = nd.gaussian_filter(rate,(4,2))
    rate = nd.gaussian_filter(rate,(2,2))
    X, Y = ds.perc, ds.pradc
    ntrial = np.array(ds.ntrial)

    return X, Y, rate, ntrial

def fig_contour_six():
    sns.set_context('paper')
    mass1 = [0.5,0.8,1.1]
    mass2 = [0.8,1.1,1.4]
    fig, axL = subplots(nrows=3,ncols=2,figsize=(8.5,9))

    i = 0
    for _mass1, _mass2 in zip(mass1,mass2):
        key = 'cp_smass={}-{}'.format(_mass1,_mass2)
        cp = ckscool.io.load_object(key,cache=1)

        sca(axL[i,0])
        df = cp.occ.plnt.copy().rename(columns={'prad':'gdir_prad','per':'koi_period'})
        ckscool.plot.planet._per_prad(df,nopoints=False,zoom=False,query=None,yerrfac=1,xerrfac=1)

        sca(axL[i,1])
        contour(cp,plot_interval=True,draw_colorbar=True,normalize=False)
        title('$M_\star = {}-{}\, M_\odot$ '.format(_mass1,_mass2))
        i+=1

    for ax in axL.flatten():
        sca(ax)
        xt = [1,3,10,30,100,300]
        yt = [1.0,1.4,2.0,2.8,4.0]
        xticks([log10(_xt) for _xt in xt],xt)
        yticks([log10(_yt) for _yt in yt],yt)
        grid()

    xlim=log10(1),log10(300)
    ylim=log10(1),log10(4)
    setp(axL,xlim=xlim,ylim=ylim)
    setp(axL[:,1:],ylabel='')
    setp(axL[:-1,:],xlabel='')
    tight_layout(True)

def fig_contour_six_sinc():
    sns.set_context('paper')
    mass1 = [0.5,0.8,1.1]
    mass2 = [0.8,1.1,1.4]
    fig, axL = subplots(nrows=3,ncols=2,figsize=(8.5,9))
    i = 0
    for _mass1, _mass2 in zip(mass1,mass2):
        key = 'cp-sinc_smass={}-{}'.format(_mass1,_mass2)
        cp = ckscool.io.load_object(key,cache=1)

        sca(axL[i,0])
        df = cp.occ.plnt.copy().rename(columns={'prad':'gdir_prad','sinc':'giso_sinc'})
        ckscool.plot.planet._sinc_prad(df,nopoints=False,zoom=False,query=None,yerrfac=1,xerrfac=1)

        sca(axL[i,1])
        contour_sinc(cp,plot_interval=True,draw_colorbar=True,normalize=False, ntrials_min=100)
        title('$M_\star = {}-{}\, M_\odot$ '.format(_mass1,_mass2))
        i+=1

    for ax in axL.flatten():
        sca(ax)
        xt = [3000,1000,300,100,30,10,3]
        yt = [1.0,1.4,2.0,2.8,4.0]
        xticks([log10(_xt) for _xt in xt],xt)
        yticks([log10(_yt) for _yt in yt],yt)
        grid()

    xlim=log10(3000),log10(3)
    ylim=log10(1),log10(4)
    setp(axL,xlim=xlim,ylim=ylim)
    setp(axL[:,1:],ylabel='')
    setp(axL[:-1,:],xlabel='')
    tight_layout(True)

class obj(object):
    def __init__(self):
        pass

def load_contour_plotter(occ):
    logperc = linspace(log10(0.1),log10(300),80)
    logpradc = linspace(log10(0.5),log10(4),80)
    perwid = 0.25
    pradwid = 0.05
    df = []
    for i in range(len(logperc)):
        for j in range(len(logpradc)):
            d = {}
            d['logperc'] = logperc[i]
            d['logpradc'] = logpradc[j]
            d['logper1'] = d['logperc'] - perwid / 2
            d['logper2'] = d['logperc'] + perwid / 2
            d['logprad1'] = d['logpradc'] - pradwid / 2
            d['logprad2'] = d['logpradc'] + pradwid / 2
            d['per1'] = 10**d['logper1']
            d['per2'] = 10**d['logper2']
            d['perc'] = 10**d['logperc']
            d['prad1'] = 10**d['logprad1']
            d['prad2'] = 10**d['logprad2']
            d['pradc'] = 10**d['logpradc']
            limits = dict([(k,d[k]) for k in 'per1 per2 prad1 prad2'.split()])
            d = dict(d,**occ.occurence_box(limits))
            df.append(d)

    df = pd.DataFrame(df)
    cp = obj()
    cp.rate = df
    cp.perwid = perwid
    cp.pradwid = pradwid
    cp.occ = occ
    return cp

def load_contour_plotter_sinc(occ):
    logsincc = linspace(log10(0.1),log10(30000),80)
    logpradc = linspace(log10(0.5),log10(4),80)
    sincwid = 0.25
    pradwid = 0.05
    df = []
    for i in range(len(logsincc)):
        for j in range(len(logpradc)):
            d = {}
            d['logsincc'] = logsincc[i]
            d['logpradc'] = logpradc[j]
            d['logsinc1'] = d['logsincc'] - sincwid / 2
            d['logsinc2'] = d['logsincc'] + sincwid / 2
            d['logprad1'] = d['logpradc'] - pradwid / 2
            d['logprad2'] = d['logpradc'] + pradwid / 2
            d['sinc1'] = 10**d['logsinc1']
            d['sinc2'] = 10**d['logsinc2']
            d['sincc'] = 10**d['logsincc']
            d['prad1'] = 10**d['logprad1']
            d['prad2'] = 10**d['logprad2']
            d['pradc'] = 10**d['logpradc']
            limits = dict([(k,d[k]) for k in 'sinc1 sinc2 prad1 prad2'.split()])
            d = dict(d,**occ.occurence_box_sinc(limits))
            df.append(d)

    df = pd.DataFrame(df)
    cp = obj()
    cp.rate = df
    cp.sincwid = sincwid
    cp.pradwid = pradwid
    cp.occ = occ
    return cp

def load_surface_smass(per1, per2):
    logsmassc = linspace(log10(0.5),log10(1.4),10)
    logsmassc = logsmassc[1:]
    logpradc = linspace(log10(0.5),log10(4),80)
    smasswid = 0.1
    pradwid = 0.05
    df = []
    for i in range(len(logsmassc)):
        d = {}
        d['logsmassc'] = logsmassc[i]
        d['logsmass1'] = d['logsmassc'] - smasswid / 2
        d['logsmass2'] = d['logsmassc'] + smasswid / 2
        d['smass1'] = 10**d['logsmass1']
        d['smass2'] = 10**d['logsmass2']
        d['smassc'] = 10**d['logsmassc']
        key = 'occur_smass={smass1:.3f}-{smass2:.3f}'.format(**d)
        occ = ckscool.io.load_object(key,cache=1)

        for j in range(len(logpradc)):
            d['logsmassc'] = logsmassc[i]
            d['logpradc'] = logpradc[j]
            d['logprad1'] = d['logpradc'] - pradwid / 2
            d['logprad2'] = d['logpradc'] + pradwid / 2
            d['prad1'] = 10**d['logprad1']
            d['prad2'] = 10**d['logprad2']
            d['pradc'] = 10**d['logpradc']
            d['nstars'] = occ.nstars
            limits = dict(per1=per1,per2=per2,prad1=d['prad1'],prad2=d['prad2'])
            d = dict(d,**occ.occurence_box(limits))
            df.append(d)

    df = pd.DataFrame(df)
    return df

def fig_contour_smass(normalize=False):
    sns.set_context('talk')
    figure(figsize=(8,6,))
    df = load_surface_smass(10,30)
    if normalize:
        df = df.query('1 < pradc < 4')
        ds = df.groupby(['smass1','prad1']).nth(0).to_xarray()
        X, Y = np.array(np.log10(ds.smassc)), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(1.5,1))
        Z = Z / Z.sum(axis=1)[:,newaxis]
        levels = arange(0.0,0.05,0.001)
        _title = 'Normalized Occurrence \n $\Delta M_\star = 0.1$ dex = 0.2 mag $\Delta R_P$ = 0.05 dex, $P$ = 10-30 day'
    else:
        ds = df.groupby(['smass1','prad1']).nth(0).to_xarray()
        X, Y = np.array(np.log10(ds.smassc)), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(1.5,1))
        levels = arange(0.0,0.08,0.003)
        _title = 'Occurrence \n $\Delta M_\star = 0.1$ dex $\Delta R_P$ = 0.05 dex, $P$ = 10-30 day'

    cmap = 'YlGn' #,None #'hot_r'
    qcs = contourf(X,Y,Z,levels=levels,cmap=cmap,zorder=0)
    colorbar()
    Z = np.array(ds.ntrial)
    cmap = sns.light_palette("gray",as_cmap=True)
    contourf(X,Y,Z,[0,25],zorder=2.5,cmap=cmap,vmax=1)
    xt = [0.5,0.7,1.0,1.4]
    yt = [1,1.4,2,2.8,4]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    ylim(log10(1,),log10(4))
    title(_title)
    xlabel('Stellar Mass (Solar Masses)')
    ylabel('Planet Size (Earth Radii)')

def load_surface_bmr(per1,per2):
    bmrc = arange(0.5,2.0,0.1)
    logpradc = linspace(log10(0.5),log10(4),80)
    bmrwid = 0.2
    pradwid = 0.05
    df = []
    for i in range(len(bmrc)):
        d = {}
        d['bmrc'] = bmrc[i]
        d['bmr1'] = d['bmrc'] - bmrwid / 2
        d['bmr2'] = d['bmrc'] + bmrwid / 2
        key = 'occur_bmr={bmr1:.3f}-{bmr2:.3f}'.format(**d)
        occ = ckscool.io.load_object(key,cache=1,verbose=0)
        for j in range(len(logpradc)):
            d['bmrc'] = bmrc[i]
            d['logpradc'] = logpradc[j]
            d['logprad1'] = d['logpradc'] - pradwid / 2
            d['logprad2'] = d['logpradc'] + pradwid / 2
            d['prad1'] = 10**d['logprad1']
            d['prad2'] = 10**d['logprad2']
            d['pradc'] = 10**d['logpradc']
            d['nstars'] = occ.nstars
            limits = dict(per1=per1,per2=per2,prad1=d['prad1'],prad2=d['prad2'])
            d = dict(d,**occ.occurence_box(limits))
            df.append(d)
    df = pd.DataFrame(df)
    return df

def fig_contour_bmr(normalize=False):
    sns.set_context('talk')
    figure(figsize=(8,6,))
    df = load_surface_bmr(10,30)

    if normalize:
        df = df.query('1 < pradc < 4')
        ds = df.groupby(['bmr1','prad1']).nth(0).to_xarray()
        X, Y = np.array(ds.bmrc), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(2,2))
        Z = Z / Z.sum(axis=1)[:,newaxis]
        levels = arange(0.0,0.06,0.003)
        _title = 'Normalized Occurrence \n $\Delta Bp-Rp$ = 0.2 mag $\Delta R_P$ = 0.05 dex, $P$ = 10-30 day'
    else:
        ds = df.groupby(['bmr1','prad1']).nth(0).to_xarray()
        X, Y = np.array(ds.bmrc), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(2,2))
        levels = arange(0.0,0.1,0.003)
        _title = 'Occurrence \n $\Delta Bp - Rp$ = 0.2 mag $\Delta R_P$ = 0.05 dex, $P$ = 10-30 day'

    qcs = contourf(X,Y,Z,levels=levels,cmap='YlGn',zorder=0)
    colorbar()
    Z = np.array(ds.ntrial)
    cmap = sns.light_palette("gray",as_cmap=True)
    contourf(X,Y,Z,[0,25],zorder=2.5,cmap=cmap,vmax=1)
    yt = [1,2,4]
    yticks([log10(_yt) for _yt in yt],yt)
    ylim(log10(1),log10(4))
    title(_title) 
    xlabel('B-V (mag)')
    ylabel('Planet Size (Earth-radii)')
