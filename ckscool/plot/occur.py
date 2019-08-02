import seaborn as sns
import pandas as pd 
from scipy import ndimage as nd
import lmfit 

from matplotlib.pylab import *
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects

sns.set_style('ticks')
sns.set_color_codes()

ptcolor = {
    'se':'g',
    'sn':'b',
    'ss':'y',
    'jup':'r',
    'sub':'b',
    'sup':'r',
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


def fig_contour_all(scale='log'):
    sns_set_style('ticks')
    fig, axL = subplots(nrows=2,ncols=2,figsize=(8,6.5))

    sca(axL[0,0])
    fig_label("a")
    contour(scale='linear',draw_colorbar=False,plot_interval=True)
    
    sca(axL[0,1])
    fig_label("b")
    caxheight = 0.25
    caxcenter = 0.8
    caxleft = 0.92
    cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
    sca(axL[0,1])    
    contour(scale='linear',plot_planetpoints=False, draw_colorbar=True,cax=cax,label=True)

    sca(axL[1,0])
    fig_label("c")
    contour(scale='log',draw_colorbar=False)
    
    sca(axL[1,1])
    fig_label("d")
    caxcenter = 0.3
    cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
    sca(axL[1,1])    
    contour(scale='log',plot_planetpoints=False, draw_colorbar=True,cax=cax,label=True)

    tight_layout(rect=[0.01,0.01,caxleft-0.01,0.99],pad=0.00)
    setp(axL[0,:],xlabel='')
    setp(axL[:,1],ylabel='')
    
def contour(cp,scale='linear', plot_planetpoints=True, plot_interval=False, 
            draw_colorbar=True,cax=None,plot_completeness=True,label=False):
    """
    Args:
       cp : contour plotter object
    
    """

    ax = gca()
    tax = gca().transAxes

    # convert into an x-array for plotting
    ds = cp.rate.groupby(['per1','prad1']).first().to_xarray()
    rate = ds.rate
    rate = rate.fillna(4e-6)

    # Smooth out the contours a bit
    #rate = nd.gaussian_filter(rate,(4,2))
    rate = nd.gaussian_filter(rate,(2,2))

    eps = 1e-10
    X, Y = np.log10(ds.perc), np.log10(ds.pradc)
    cmap = None
    cmap = 'YlGn' #,None #'hot_r'
    #cmap = 'hot_r'
    levels = None
    cbarlabel=''
    if scale=='linear':
        levels = arange(0,4.0e-2+eps,0.0025) 
        Z = rate
        kw = dict(extend='neither',cmap=cmap)
        cbarticks = levels[::2]
        cbarticklabels = ["{:.0f}".format(1e2*_yt) for _yt in cbarticks]
        kw = dict(levels=levels,extend='neither',cmap=cmap)

    # log scale to show Hot-J
    if scale=='log':
        Z = np.log10(rate)
        #levels = np.arange(-3.75,-1.5+eps,0.125) 
        #cbarticks = levels[::2]
        #cbarticklabels = ["{:.2f}".format(1e2*10**_yt) for _yt in cbarticks]
        levels = np.arange(-4,-1+eps,0.25) 
        cbarticklabels = [0.01, 0.03, 0.1, 0.3, 1, 3,10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        kw = dict(extend='min',cmap=cmap,levels=levels,vmin=-3.99)

    cbarlabel = r"""Planets per 100 Stars per $P-R_P$ interval"""

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    qcs = contourf(X,Y,Z, **kw)

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

    if plot_planetpoints:
        cks = cksmet.io.load_table('cks-cuts',cache=1)
        cks = cks[~cks.isany]
        x,y = log10(cks.koi_period),log10(cks.iso_prad)
        plot(x,y,'.',mfc='none',mec='red',mew=0.3,ms=5)

    # Completeness
    if plot_completeness:
        Z = np.array(ds.ntrial)
        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(X,Y,Z,[0,50],zorder=2.5,cmap=cmap,vmax=1)
        '''
        text(
            0.95,0.15,'Low Completeness',rotation=12,zorder=5,size='small',
            transform=ax.transAxes,ha='right'
        )
        '''
    if plot_interval:
        #inv = ax.transAxes.inverted()
        xyaxes = (0.1,0.9)
        xy = ax.transLimits.inverted().transform(xyaxes)
        #xy = ax.transAxes.transform((0.9,0.9))
        w = cp.perwid
        h = cp.pradwid
        rect = Rectangle(xy, w, h,lw=1,ec='r',fc='none',zorder=4)
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

    if label:
        if scale=='linear':
            kw = dict(
                size='x-small',zorder=5,va='center',ha='center',color='red'
            )
            text(log10(30),log10(3),'Warm Sub-Neptunes', **kw)
            text(log10(20),log10(1),'Warm Super-Earths',**kw)
            text(log10(30),log10(1.7),'Radius Gap',rotation=-10,**kw)
            text(log10(150),log10(16),'Cool Jupiters',**kw)

        if scale=='log':
            kw = dict(
                size='x-small',zorder=5,va='center',ha='center',color='red'
            )
            text(log10(3),log10(20),'Hot Jupiters',**kw)
            x = [1,3.0,15,1]
            y = [1.7,1.7,10.0,10.0]
            x = np.log10(np.array(x))
            y = np.log10(np.array(y))
            plot(x,y,linestyle='--',color='red',lw=1)
            kw['va'] = 'top'
            text(log10(3),log10(10)-0.03,'Hot Planet Desert',**kw)

    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16,32]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    xlim(log10(1),log10(300))
    ylim(log10(0.5),log10(32))
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    
def fig_checkerboard(plot_planetpoints=True):
    key = 'occur-per=0.25-prad=twoperoctave'
    occ = cksmet.io.load_object(key,cache=1)
    sns.set_context('paper',font_scale=1.0)
    sns.set_style('ticks')
    fig,ax = subplots(figsize=(7,5.5))
    loglog()

    logper1 = np.log10(1)
    logper2 = np.log10(350)
    logprad1 = np.log10(0.5)
    logprad2 = np.log10(32)
    eps = 1e-10
    bin_per = 0.25 # period bin width in dex
    bin_prad = 0.5 * log10(2) # prad bin width in dex
    per_bins = 10**(np.arange(logper1,logper2+eps,bin_per))
    prad_bins = 10**(np.arange(logprad1,logprad2+eps,bin_prad))
    eps = 1e-3
    smet_bins = [-1,0.5]

    occ = cksmet.analysis.compute_occur_bins(
        per_bins, prad_bins, smet_bins
    )

    prob_det_min = 0.1
    cmap='YlGn'
    # Completeness contour
    # levels = [prob_det_min,prob_det_min + epsilon]
    # contour(dsf.perc, dsf.pradc, dsf.prob_det, levels=levels,color='k')
    ds = occ.df.groupby(['perc','pradc']).first().to_xarray()
    ds = ds.transpose('pradc','perc')

    per = np.hstack([np.array(ds.per1[0]),np.array(ds.per2[0,-1])])
    prad = np.hstack([np.array(ds.prad1[:,0]),np.array(ds.prad2[-1,0])])
    X,Y = meshgrid(per,prad)

    ds = ds.where((ds.per2 < 350) & (ds.prob_det_mean > 0.25))
    ar = ma.masked_invalid(np.array(ds.rate))
    ar = log10(ar)
    qcs = pcolormesh(X,Y,ar,cmap=cmap,vmin=-4,vmax=-1)

    if plot_planetpoints:
        cks = cksmet.io.load_table('cks-cuts',cache=1)
        cks = cks[~cks.isany]
        x,y = cks.koi_period,cks.iso_prad
        plot(x,y,'.',mfc='none',mec='red',mew=0.5,ms=5)

    if 1:
        caxheight = 0.5
        caxcenter = 0.5
        caxleft = 0.88
        cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
        sca(ax)    

        cbarticklabels = [0.01,0.03, 0.1, 0.3, 1, 3, 10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        cbar = colorbar(qcs,cax=cax,ticks=cbarticks)
        #cbar = colorbar(qcs,cax=cax)
        t = cbar.ax.set_yticklabels(cbarticklabels) 
        setp(t,size='x-small')
        cbarlabel = r"""Planets per 100 Stars per Bin"""
        cbar.set_label(cbarlabel,size='small')

        #plot(occ.plnt.per,occ.plnt.prad,'o',ms=5,color='Tomato')
    df = ds.to_dataframe()
    annotate_checkerboard(df)
    yt = [0.5, 1, 2, 4, 8, 16, 32]
    xt = [0.3,1, 3, 10, 30, 100, 300]
    xticks(xt,xt)
    yticks(yt,yt)
    xlim(1,300)
    ylim(0.5, 32)
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    minorticks_off()
    tight_layout(rect=[0.01,0.01,0.85,0.99],pad=0)
    return occ

def annotate_checkerboard(df,print_ul=False):
    ax = gca()
    for i, row in df.iterrows():
        xy = np.array([row.per1, row.prad1])
        width = row.per2 - row.per1
        height = row.prad2 - row.prad1
        if ~np.isnan(row.rate_ul):
            if print_ul:
                txt = "<{:.2f}".format(row.rate_ul * 100)
            else:
                txt=""
            prate = 0
        else:
            prate = row.rate * 100
            txt = "{:.2f}".format(prate)

        at =  matplotlib.patches.Rectangle(
            xy, width, height,ec='LightGray',lw=0.5, fc='none'
        )

        ax.add_artist(at)

        text = ax.annotate(txt, xy=xy, xytext=xy*1.05,size='x-small')
        text.set_path_effects(
            [
                path_effects.Stroke(linewidth=1.5, foreground='white'),
                path_effects.Normal()
            ]
        )

def fig_per_smet():
    sns.set_context('paper',font_scale=1.0)
    sns.set_style('ticks')
    fig,axL = subplots(ncols=2,figsize=(7.0,3.5))
    sca(axL[0])

    xk = 'perc'
    fac = 1.0

    # Super Earths Metal rich
    sizes = 'se se sn sn'.split()
    smets = 'sub sup sub sup'.split()
    querys = ['1 <  perc < 100']*2  + ['1 < perc < 300']*2
    iax = [0,0,1,1]
    
    dper = 0.25
    for i in range(4):
        sca(axL[iax[i]])
        size = sizes[i]
        smet = smets[i]
        query = querys[i]
        occkey = 'occur-per={}-prad=physical-smet={}'.format(dper,smet)
        fitkey = 'fit_per-{}-{}'.format(smet,size)
        occ = cksmet.io.load_object(occkey)
        cut = occ.df.query(query)
        sampler = cksmet.plotting.occur.SamplerPer(fitkey, smet, dper)
        sampler.plot_all()
        plot_rates(xk,cut.ix[size],smet,fac=fac)

    letters = ['a','b']
    titles = ['Super-Earths','Sub-Neptunes']
    for i in range(2):
        ax = axL[i] 
        sca(ax)
        loglog()
        fig_label(letters[i])
        yticks_planets_per_100_stars()
        xt = [1,3,10,30,100,300]
        xticks(xt,xt)
        title(titles[i])

        i = 0
        for smet in 'sub sup'.split():
            s =  "\n"*i + namedict[smet] 
            kw=dict(color=ptcolor[smet],va='top',ha='left',transform=ax.transAxes,size='small')
            text(0.75, 0.12, s, **kw)
            i+=1

    setp(axL, xlabel='Period (days)', ylim=(3e-4,3e-1), xlim = (1,300))
    setp(axL, ylabel='Planets per 100 Stars per 0.25 dex $P$ Interval')

    fig.set_tight_layout(True)


def fig_summary():
    sns_set_style('ticks')
    fig,ax = subplots(figsize=(5,4))
    perbins = [1, 10, 100]
    pers = ['hot','warm']
    sizes = ['se','sn','ss','jup']
    pradbins = [1.0, 1.7, 4.0, 8.0, 24.0]
    loglog()
    d = cksmet.values.val_samp(return_dict=True)
    for i in range(len(perbins)-1):
        for j in range(len(pradbins)-1):
            per = pers[i]
            size = sizes[j]
            xy = perbins[i], pradbins[j]
            h = pradbins[j+1] - pradbins[j]
            w = perbins[i+1] - perbins[i]
            xytext = xy[0],xy[1]+h
            
            print xy
            fitkey = 'fit_persmet-{}-{}'.format(per,size)
            fit = cksmet.io.load_object(fitkey,cache=1)

            rect = Rectangle(xy, w, h, lw=1,ec='r',fc='none',zorder=10)
            ptype = "%s %s" % (per.capitalize(),namedict[size])
            s = ""
            s+= "%s \n" % ptype
            s+=r"$\left<\mathrm{[Fe/H]}\right>$ = $%s \pm %s$" % (d[ptype+" mean"],d[ptype+" sem"])
            s+="\n"
            pstr = fit.to_string()
            for _pstr in pstr:
                if _pstr.count('beta')==1:
                    _pstr = _pstr.replace(r'{beta}','')
                    _pstr = _pstr.replace(r"{$","").replace("$}","")
                    _pstr = r"$\beta = %s$" % _pstr

            if fitkey=='fit_persmet-warm-jup':
                _pstr = r"$\beta$ = Unconstrained" 

            s+=_pstr 


            ax.annotate(s, xy=xytext, xytext=(3, -2.5),textcoords='offset points',size='x-small',va='top')

            #text(xy[0],xy[1],s,zorder=10)
            ax.add_patch(rect)

    yt = [1, 1.7, 4.0, 8.0, 24.0]
    xt = [1, 10, 100]
    xticks(xt,xt)
    yticks(yt,yt)
    xlabel('Orbital Period (days)')
    ylabel('Planet size (Earth-radii)')
    setp(ax,xlim=(1,100),ylim=(1,24))
    minorticks_off()
    fig.set_tight_layout(True)

def fig_per():
    sns.set_context('paper',font_scale=1.1)
    fig,axL = subplots(ncols=1,sharex=True,sharey=True,figsize=(5,4))
    loglog()
    xk = 'perc'
    dper = 0.25
    # Super Earths
    key = 'occur-per=0.25-prad=physical-smet=all'.format(dper)
    occ = cksmet.io.load_object(key,cache=1)
    cut = occ.df.query('1 <  perc < 100')
    size='se'
    sampler = cksmet.plotting.occur.SamplerPer('fit_per-all-se', size, dper)
    sampler.plot_all()
    fac = 1
    plot_rates(xk,cut.ix[size],size,fac=fac)

    # Sub Neptunes
    key = 'occur-per=0.25-prad=physical-smet=all'.format(dper)
    occ = cksmet.io.load_object(key,cache=1)
    cut = occ.df.query('1 <  perc < 350')
    size='sn'
    sampler = cksmet.plotting.occur.SamplerPer('fit_per-all-sn', size, dper)
    sampler.plot_all()
    plot_rates(xk, cut.ix[size], size, fac=fac)

    # Just plot the rates for Sub-Sat and hot Jup
    size='ss'
    plot_rates(xk, cut.ix[size], size, fac=fac)

    # Sub Satruns 
    size='jup'
    plot_rates(xk, cut.ix[size], size, fac=fac)

    ylabel('Planets per 100 Stars per 0.25 dex $P$ Interval')
    yticks_planets_per_100_stars()


    xt = [1,3,10,30,100,300]
    xticks(xt,xt)
    yticks_planets_per_100_stars()
    setp(axL, xlabel='Period (days)', ylim=(1e-4,0.3), xlim = (1,300))

    i = 0
    for size in 'se sn ss jup'.split():
        s =  "\n"*i + namedict[size] 
        kw=dict(color=ptcolor[size],va='top',ha='left',transform=axL.transAxes,size='small')
        text(0.05, 0.95, s, **kw)
        i+=1

    fig.set_tight_layout(True)

def plot_rates(xk, occur, fmtkey,fac=1.0, **kw):
    """
    Args
        xk (str): x value
        occur (pd.DataFrame): must contain rate, rate_err1, rate_err2
    """
    if xk=='smetc':
        _ptcolor = ptcolor[fmtkey]
    if xk=='perc':
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

    semilogy()

    x = occur[xk]
    y = occur.rate
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw1)
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw2)

    occurul = occur.dropna(subset=['rate_ul'])
    if len(occurul) >0:
        plot(x,occur.rate_ul*fac,**ulkw)

class Sampler(object):
    nsamples = 1000
    def __init__(self, key, fmtkey, dx):
        """
        dx : size of phase-space volume to integrate occurrence 
        """
        self.fit = cksmet.io.load_object(key)
        self.fmtkey = fmtkey
        self.dx = dx

    def plot_band(self):
        p16, p50, p84 = np.percentile(self.fit_samples,[16,50,84],axis=0)
        _bdcolor = ptcolor[self.fmtkey]
        fill_between(self.x, p16, p84, color=_bdcolor,alpha=0.3)
        
    def plot_best(self):
        plot(self.x,self.fit_best,color=ptcolor[self.fmtkey])

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

class SamplerSmet(Sampler):
    x = np.arange(-0.5,0.6,0.02) 
    
class SamplerPerSmet(Sampler):
    smet1 = -0.4
    smet2 = 0.4
    x = np.arange(-0.5,0.6,0.2) 
    dsmet = 0.2

    def __init__(self, key, fmtkey, perlims, smet_field):
        self.fit = cksmet.io.load_object(key)
        self.fmtkey = fmtkey
        self.smet_field = smet_field
        self.perlims = perlims

    def compute_best(self):
        p = self.fit.pfit.valuesdict()
        self.fit_best = self.model(p)

    # Overload the model method
    def model(self, params):
        fit_sample = []
        for _smeti in self.x:
            _smet1 = _smeti - 0.5*self.dsmet
            _smet2 = _smeti + 0.5*self.dsmet
            smet = self.smet_field[self.smet_field.between(_smet1,_smet2)]
            val = cksmet.fit.per_powerlaw_smet_exp_integral(params, self.perlims, smet)
            fit_sample.append(val)
        fit_sample = np.array(fit_sample)
        return fit_sample

class SamplerPer(Sampler):
    dx = 0.25 / 10
    x = arange(np.log10(0.1) + 0.5*dx,np.log10(1000),dx)
    x = 10**x

def sum_cells_per(df0):
    df2 = []
    for smetc in df0.smetc.drop_duplicates():
        df = df0[df0.smetc==smetc] 
        rate = cksmet.stats.sum_cells(df.ntrial,df.nplnt)
        rate['smetc'] = smetc
        df2.append(rate)
    
    df2 = pd.DataFrame(df2)
    return df2
    
def fig_smet_warm():
    key = 'occur-per=hotwarm-prad=physical-smet=0.2'
    occ = cksmet.io.load_object(key,cache=1)
    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]
    xk = 'smetc'
    dist = 'warm'

    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    smet_field = lamo.lamo_smet

    sizes = 'se sn ss jup'.split()
    plot_fits = [True, True, True, False]
    for i in range(4):
        size = sizes[i]
        plot_fit = plot_fits[i]
        cut = df.ix[dist,size]
        plot_rates(xk, cut, size)
        fitkey = 'fit_persmet-{}-{}'.format(dist, size)
        sampler = SamplerPerSmet(fitkey, size, smet_field)
        if plot_fit:
            sampler.plot_all()

    ylabel('Planets per 100 Stars')
    ylim(1e-3,1)

def fig_smet_hot():
    """
    Display the integrated planet occurrence from P = 1-10 days 
    """

    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    smet_field = lamo.lamo_smet

    xk = 'smetc'
    key = 'occur-per=0.25-prad=physical-smet=0.2'
    occ = cksmet.io.load_object(key)

    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]
    sizes = 'se sn ss jup'.split()
    for i in range(4):
        size = sizes[i]
        cut = df.ix[size].query('1 < perc < 10')
        cut = sum_cells_per(cut) 
        plot_rates(xk, cut, size)
        fitkey = 'fit_persmet-hot-{}'.format(size)
        sampler = SamplerPerSmet(fitkey, size, smet_field)
        sampler.plot_all()

    ylabel('Planets per 100 Stars')
    ylim(1e-5,1)

def yticks_planets_per_100_stars():
    yt = np.array([1e-5,3e-5,1e-4,3,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1])
    syt = map(lambda x: x if x < 1 else "%.0f" % x , yt*100)
    yticks(yt,syt)

def fig_smet():
    """Figure showing total occurrence in hot and cold regions for
    different planet classes. Also show the best fit curves.  Because
    the fits are surfaces in smet and perc we cannot display on a 2D
    graph. We integrate over period to compute the occurrence just as
    a function of metallicity. The bins are added up.
    """
    sns_set_style('ticks')
    fig, axL = subplots(ncols=2,nrows=1,figsize=(7,5),sharex=True)
    sca(axL[0])

    xk = 'smetc'
    key = 'occur-per=0.25-prad=physical-smet=0.2'
    occ = cksmet.io.load_object(key)
    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]

    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    smet_field = lamo.lamo_smet

    sizes = 'se sn ss jup se sn ss jup'.split()
    dists = 'hot hot hot hot warm warm warm warm'.split()
    plot_fits = [True, True, True, True, True, True, True, False]
    for i in range(8):
        size = sizes[i]
        plot_fit = plot_fits[i]
        dist = dists[i]
        if dist=='hot':
            perlim = [1,10]
            sca(axL[0])
        elif dist=='warm':
            perlim = [10,100]
            sca(axL[1])
            
        cut = df.ix[size]
        cut = cut[cut.perc.between(*perlim)]
        cut = sum_cells_per(cut) 
        plot_rates(xk, cut, size)
        fitkey = 'fit_persmet-{}-{}'.format(dist, size)
        sampler = SamplerPerSmet(fitkey, size, perlim, smet_field)
        if plot_fit:
            sampler.plot_all()

    labels = ["$P$ = 1$-$10 days","$P$ = 10$-$100 days"]
    for i in range(2):
        ax = axL[i] 
        sca(ax)
        letters = ["a","b"]
        fig_label(letters[i])
        add_anchored(labels[i],1,frameon=False)

        j = 0
        for size in 'se sn ss jup'.split():
            s =  "\n"*j + namedict[size] 
            text(0.65, 0.15, s, color=ptcolor[size],va='top',ha='left', transform=ax.transAxes,size='small')
            j+=1

    for ax in axL:
        sca(ax)
        yticks_planets_per_100_stars()
        ylabel('Planets per 100 Stars')
        ylim(1e-4,3)
        xlim(-0.4,0.4)
        xlabel('[Fe/H]')

    fig.set_tight_layout(True)


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
