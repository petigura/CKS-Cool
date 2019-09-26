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

def contour(cp, plot_interval=False, draw_colorbar=True, cax=None, 
            plot_completeness=True, ntrials_min=50, levels=None):
    """
    Args:
       cp : contour plotter object

    """

    ax = gca()
    tax = gca().transAxes
    cmap = 'YlGn' #,None #'hot_r'

    # convert into an x-array for plotting
    ds = cp.rate.groupby(['per1','prad1']).first().to_xarray()
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
            (ds.ntrial > ntrials_min) 
            & (ds.pradc > 1.0) 
            & (ds.pradc < 4.0)
        )
        b = array(b)
        maxz = np.max(rate[b])
        maxz = np.round(maxz*1.1,3)
        levels = linspace(0,maxz+eps,14)

    cbarticks = levels[::2]
    cbarticklabels = ["{:.1f}".format(1e2*_yt) for _yt in cbarticks]
    kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)
    cbarlabel = r"""Planets per 100 Stars per $P-R_P$ interval"""
    qcs = contourf(X, Y, Z, **kw)

    # Completeness
    if plot_completeness:
        Zt = np.array(ds.ntrial)
        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(X,Y,Zt,[0,ntrials_min],zorder=2.5,cmap=cmap,vmax=1)
        '''
        text(
            0.95,0.15,'Low Completeness',rotation=12,zorder=5,size='small',
            transform=ax.transAxes,ha='right'
        )
        '''


    # plot straight contours
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


    return X, Y, Z

#        text(xyaxes[0]+0.07,xyaxes[1],s,**kw)
# ---------------------------------------------------------------------------- #

def contour_sinc(cp, plot_interval=False, draw_colorbar=True,cax=None,
                 plot_completeness=True,label=False, normalize=False, 
                 ntrials_min=50):
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
    Z = np.array(rate)
    cmap = 'YlGn' #,None #'hot_r'
    levels = None
    cbarlabel=''
    if levels==None:
        #levels = arange(0,5e-2+eps,0.0025) 
        b = (
            (ds.ntrial > ntrials_min) 
            & (ds.pradc > 1.0) 
            & (ds.pradc < 4.0)
        )
        b = array(b)
        maxz = np.max(rate[b])
        maxz = np.round(maxz*1.1,3)
        levels = linspace(0,maxz+eps,14)

    cbarticks = levels[::2]
    cbarticklabels = ["{:.1f}".format(1e2*_yt) for _yt in cbarticks]
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

    xt = [1,3,10,30,100,1000,10000]
    yt = [0.5,1,2,4,8,16,32]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    xlim(log10(10000),log10(1))
    ylim(log10(1),log10(4))
    xlabel('Stellar Incident Flux (Earth Units)')
    ylabel('Planet Size (Earth-radii)')

    if plot_interval:
        #inv = ax.transAxes.inverted()
        xy = 3.1, 0.52
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
        df = cp.occ.plnt.copy()
        df = df.rename(columns={'prad':'gdir_prad','per':'koi_period'})
        ckscool.plot.planet._per_prad(
            df, nopoints=False, zoom=False, query=None, yerrfac=1, xerrfac=1
        )
        sca(axL[i,1])
        contour(cp,plot_interval=True,draw_colorbar=True)
        title = '$M_\star = {}-{}\, M_\odot$ '.format(_mass1,_mass2)
        setp(axL[i,:],title=title)
        i+=1

    for ax in axL.flatten():
        sca(ax)
        xt = [1,3,10,30,100,300]
        yt = [1.0,1.4,2.0,2.8,4.0]
        xticks([log10(_xt) for _xt in xt],xt)
        yticks([log10(_yt) for _yt in yt],yt)
        grid()

    xlim = log10(1),log10(300)
    ylim = log10(1),log10(4)
    setp(axL, xlim=xlim, ylim=ylim)
    setp(axL[:,1:], ylabel='')
    setp(axL[:-1,:], xlabel='')
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
        df = cp.occ.plnt.copy()
        df = df.rename(columns={'prad':'gdir_prad','sinc':'giso_sinc'})
        ckscool.plot.planet._sinc_prad(df,nopoints=False,zoom=False,query=None,yerrfac=1,xerrfac=1)

        sca(axL[i,1])
        contour_sinc(cp,plot_interval=True,draw_colorbar=True, ntrials_min=100)
        title = '$M_\star = {}-{}\, M_\odot$ '.format(_mass1,_mass2)
        setp(axL[i,:],title=title)
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

def fig_contour_smass(per1, per2, normalize=False,):
    sns.set_context('talk')
    figure(figsize=(8,6,))
    df = load_surface_smass(per1,per2)
    if normalize:
        df = df.query('0.25 < pradc < 4')
        ds = df.groupby(['smass1','prad1']).nth(0).to_xarray()
        X, Y = np.array(np.log10(ds.smassc)), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(2,4))
        Z = Z / Z.sum(axis=1)[:,newaxis]
        levels = np.round(linspace(0.0,1.2*Z.max(),20),3)
        _title = 'Normalized Occurrence \n $\Delta M_\star = 0.1$ dex $\Delta R_P$ = 0.05 dex, $P$ = {}-{} day'.format(per1,per2)
    else:
        ds = df.groupby(['smass1','prad1']).nth(0).to_xarray()
        X, Y = np.array(np.log10(ds.smassc)), np.array(np.log10(ds.pradc))
        Z = ds.rate.fillna(0)
        Z = nd.gaussian_filter(Z,(1,4))
        
        # Choose levels based on the maximum flu
        Zcut = Z[array((1 < ds.pradc) & (ds.pradc < 4))]
        levels = linspace(0,1.2*Zcut.max(),20)
        _title = 'Occurrence \n $\Delta M_\star = 0.1$ dex $\Delta R_P$ = 0.05 dex, $P$ = {}-{} day'.format(per1,per2)

    cmap = 'YlGn' #,None #'hot_r'
    qcs = contourf(X,Y,Z,levels=levels,cmap=cmap,zorder=0)
    colorbar()
    Z = np.array(ds.ntrial)
    cmap = sns.light_palette("gray",as_cmap=True)
    contourf(X,Y,Z,[0,100],zorder=2.5,cmap=cmap,vmax=1)
    xt = [0.5,0.7,1.0,1.4]
    yt = [1,1.4,2,2.8,4]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)

    df = ckscool.io.load_table('planets-cuts2+iso')
    df = df[~df.isany]
    df = df[df.koi_period.between(per1,per2)]
    plot(log10(df.giso_smass), log10(df.gdir_prad),'.')

    ylim(log10(0.5),log10(4))
    xlim(log10(0.5),log10(1.4))
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


def per(per1, per2, prad1, prad2, smass1, smass2, fmtkey):
    dlogper_bin = 0.25
    dlogper_fit = 0.05 # Size of the bins used in the fitting

    key = 'occur_smass={}-{}'.format(smass1,smass2)
    occ = ckscool.io.load_object(key,cache=1)

    dx = [dlogper_fit]
    cut = occ.occurrence_grid(
        per1=per1, per2=per2, dlogper=dlogper_fit, 
        prad1=prad1, prad2=prad2
    )

    occ = ckscool.io.load_object(key,cache=1)
    key = 'fitper_per={}-{}_prad={}-{}_smass={}-{}'.format(
        per1,per2,prad1,prad2,smass1,smass2
    )
    fit = ckscool.io.load_object(key,cache=1)
    sampler = ckscool.plot.occur.SamplerPer(fit, fmtkey, dlogper_bin)
    sampler.plot_all()

    # Plot binned occurrence to guide the eye
    df = occ.occurrence_grid(
        per1=per1, per2=per2, dlogper=dlogper_bin, prad1=prad1, prad2=prad2
    )
    fac = 1
    plot_rates('perc', df, fmtkey, fac=fac)

def sinc(sinc1, sinc2, prad1, prad2, smass1, smass2, fmtkey):
    dlogsinc_bin = 0.5
    dlogsinc_fit = 0.05 # Size of the bins used in the fitting

    key = 'occur-sinc_smass={}-{}'.format(smass1,smass2)
    occ = ckscool.io.load_object(key,cache=1)

    dx = [dlogsinc_fit]
    cut = occ.occurrence_grid(
        sinc1=sinc1, sinc2=sinc2, dlogsinc=dlogsinc_fit, 
        prad1=prad1, prad2=prad2
    )

    occ = ckscool.io.load_object(key,cache=1)
    key = 'fitsinc_sinc={}-{}_prad={}-{}_smass={}-{}'.format(
        sinc1,sinc2,prad1,prad2,smass1,smass2
    )
    fit = ckscool.io.load_object(key,cache=1)
    sampler = ckscool.plot.occur.SamplerSinc(fit, fmtkey, dlogsinc_bin)
    sampler.plot_all()
    # Plot binned occurrence to guide the eye
    df = occ.occurrence_grid(
        sinc1=sinc1, sinc2=sinc2, dlogsinc=dlogsinc_bin, prad1=prad1, prad2=prad2
    )
    fac = 1
    plot_rates('sincc', df, fmtkey, fac=fac)

def plot_rates(xk, occur, fmtkey, fac=1.0, **kw):
    """
    Args
        xk (str): x value
        occur (pd.DataFrame): must contain rate, rate_err1, rate_err2
    """
    _ptcolor = ptcolor[fmtkey]

    efac = 0.2
    cfac = 1.05
    ebkw = {}

    kw['ms'] = 5
    kw['mew'] = efac * kw['ms']
    kw['mfc'] = 'w'
    kw['capsize'] = cfac * kw['ms']
    kw['capthick'] = kw['mew']
    kw['zorder'] = 4
    kw['color'] = _ptcolor
    kw['fmt'] = 'o'

    ebkw1 = kw
    ebkw2 = dict(**ebkw1)
    ebkw2['mfc'] = 'none'
    ebkw2['zorder'] = 5
    ebkw2['lw'] = 0
    ebkw2['zorder'] = 5

    # Points as upperlimit remove errorbar specific kw
    ulkw = dict(**kw)
    ulkw.pop('fmt')
    ulkw.pop('capthick')
    ulkw.pop('capsize')
    ulkw['color'] = _ptcolor
    ulkw['marker'] = 'v'
    ulkw['lw'] = 0
    ulkw['zorder'] = 6
    
    yerr = np.array(occur['rate_err2 rate_err1'.split()]).T
    yerr[0] *= -1 

    #semilogy()
    
    x = occur[xk]
    y = occur.rate
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw1)
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw2)

    occurul = occur.dropna(subset=['rate_ul'])
    if len(occurul) >0:
        plot(x,occur.rate_ul*fac,**ulkw)

sns.set_style('ticks')
sns.set_color_codes()

ptcolor = {
    'se':'g',
    'sn':'b',
    'ss':'y',
    'jup':'r',
    'smass1':'r',
    'smass2':'g',
    'smass3':'b',
}

bdcolor = {
    'se':'light green',
    'sn':'light blue',
    'ss':'light mustard',
    'jup':'light pink'
}

namedict = {
    'se':'Super-Earths',
    'sn':'Sub-Neptunes',
    'ss':'Sub-Saturns',
    'jup':'Jupiters',
    'sub':'[Fe/H] < 0',
    'sup':'[Fe/H] > 0',
}

sizedict = {
    'se':'$R_P$ = 1.0$-$1.7 $R_E$',
    'sn':'$R_P$ = 1.7$-$4.0 $R_E$',
    'ss':'$R_P$ = 4.0$-$8.0 $R_E$',
    'jup':'$R_P$ = 8.0$-$24.0 $R_E$'
}

class Sampler(object):
    nsamples = 1000
    def __init__(self, fit, fmtkey, dx):
        """
        dx : size of phase-space volume to integrate occurrence 
        """
        self.fit = fit
        self.fmtkey = fmtkey
        self.dx = dx

    def plot_band(self):
        p16, p50, p84 = np.percentile(self.fit_samples,[16,50,84],axis=0)
        _bdcolor = ptcolor[self.fmtkey]
        fill_between(self.x, p16, p84, color=_bdcolor,alpha=0.3)
        
    def plot_best(self):
        plot(self.x, self.fit_best, color=ptcolor[self.fmtkey])

    def plot_all(self):
        self.compute_samples()
        self.compute_best()
        self.plot_band()
        self.plot_best()

    def dict_to_params(self, params):
        _params = lmfit.Parameters()
        for k in params.keys():
            _params.add(k, value=params[k] )
        return _params 

    def compute_samples(self):
        psamples = self.fit.sample_chain(self.nsamples)
        fit_samples = []
        for i, params in psamples.iterrows():
            params = self.dict_to_params(params)
            fit_samples.append(self.model(params))
        self.fit_samples = np.vstack(fit_samples)

    def model(self,params):
        """Compute planet occurrence integrated over volumes of size dx""" 
        return self.fit.model(params,self.x) * self.dx 

    def compute_best(self):
        params = self.fit.pfit
        self.fit_best = self.model(params)

class SamplerPer(Sampler):
    dx = 0.25 / 10
    x = arange(np.log10(0.1) + 0.5*dx,np.log10(1000),dx)
    x = 10**x

class SamplerSinc(Sampler):
    dx = 0.25 / 10
    x = arange(np.log10(0.1) + 0.5*dx,np.log10(1e4),dx)
    x = 10**x
