import lmfit
from matplotlib.pylab import *
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
import pandas as pd
import seaborn as sns
from scipy import ndimage as nd

import ckscool.io
import ckscool.plot.planet
import ckscool.gradient 
from .planet import NDPlotter
from .config import *

sns.set_style('ticks')
sns.set_color_codes()

# Functions that make pdf plots

def fig_contour_six_per(plot_gradient=False):
    sns.set_context('paper', font_scale=1.0)
    pl = SixPlotterPerPrad()
    pl.occur_prefix = 'occur-per-prad'
    pl.plot(plot_gradient)

def fig_contour_six_sinc(plot_gradient=False):
    sns.set_context('paper')
    pl = SixPlotterSincPrad()
    pl.plot(plot_gradient)

def fig_mean_planet_size(annotate=True):
    sns.set_context('paper')
    fig, axL = subplots(nrows=1,figsize=(3.5,3.5),sharex=True)

    loglog()
    mps = ckscool.io.load_object('mps_size-se',cache=1)
    errorbar(mps.smassc,mps.mn,yerr=mps.std,fmt='o',color='g')

    xi = 10**mps.logsmassci
    fill_between(xi, mps.q16,mps.q84 ,alpha=0.5,color='g')


    mps = ckscool.io.load_object('mps_size-sn',cache=1)
    errorbar(mps.smassc,mps.mn,yerr=mps.std,fmt='o',color='b')
    fill_between(xi, mps.q16,mps.q84 ,alpha=0.5,color='b')

    yt = np.round(arange(1.0,3.201,0.2),1)
    xt = [0.5,0.7,1.0,1.4] 

    minorticks_off()
    yticks(yt,yt)
    xticks(xt,xt)

    setp(axL,ylabel='Mean Planet Size (Earth-radii)')
    text(0.5,2.8,'Sub-Neptunes (P < 100)')
    text(0.5,1.4,'Super-Earths (P < 30)')

    #setp(axL,title='Super-Earths (P < 30)')
    setp(axL,xlabel='Stellar Mass (Solar-Masses)')
    ylim(1.0,3.2)

    if annotate:
        plot(xi, 2.5 * xi**0.25,'k--')
        text(1.35,1.55,r'$\alpha=0.25$',ha='right',size='small',rotation=16)
        plot(xi, 1.4 * xi**0.25,'k--')
        text(1.35,2.75,r'$\alpha=0.25$',ha='right',size='small',rotation=16)

        
    tight_layout()

#def fig_mean_planet_size_ann()

def fig_occur_per():
    sns.set_context('paper',font_scale=1.2)
    fig, axL = subplots(ncols=2,figsize=(8,3.75),sharey=True)
    sizes = ['sn','se']
    for i in range(2):
        loglog()
        size = sizes[i]
        ax = axL[i] 
        sca(ax)
        smass = [0.5, 0.7, 1.0, 1.4]
        if size=='se':
            prad1, prad2 = 1.0, 1.7
            _title = 'Super-Earths'
            xlim_points = [1,100]
        elif size=='sn':
            prad1, prad2 = 1.7, 4.0
            _title = 'Sub-Neptunes'
            xlim_points = [1,300]

        _title += ' ($R_p$ = {}$-${} $R_\oplus$)'.format(prad1,prad2) 
            
        dlogper = 0.25
        
        for i in range(3):
            smass1, smass2 = smass[i],smass[i+1]
            key = 'fitper_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            fmtkey = 'smass{}'.format(i+1)
            plot_per_rates(
                fit, log10(xlim_points[0]), log10(xlim_points[1])+0.1,
                dlogper,fmtkey
            )
            if size=='sn':
                s = '$M_\star$ = {}$-${} $M_\odot$'.format(smass1, smass2)
                text(0.6, 0.15 - i*0.05, s, color=ptcolor[fmtkey], transform=ax.transAxes, size='small')                
        title(_title)
            
    sca(axL[0])
    setp(axL, xlabel='Orbital Period (days)', ylim=(1e-4,1e0))
    setp(
        axL[0],
        ylabel = 'Planets per Star'.format(dlogper),
        xlim =(1,300),
    )
    setp(axL[1], xlim =(1,300))

    for i in range(2):
        sca(axL[i])
        fig_label('ab'[i])
    
    tight_layout(True)

def fig_occur_sinc():
    sns.set_context('paper',font_scale=1.2)
    fig, axL = subplots(ncols=2,figsize=(8,3.75),sharey=True)
    sizes = ['sn','se']
    for i in range(2):
        size = sizes[i]
        sca(axL[i])
        smass = [0.5,0.7,1.0,1.4]
        if size=='se':
            prad1, prad2 = 1.0, 1.7
            _title = 'Super-Earths'
            xlim_points = [1,1e4]
        elif size=='sn':
            prad1, prad2 = 1.7, 4.0
            _title = 'Sub-Neptunes'
            xlim_points = [1,1e4]

        _title += ' ($R_p$ = {}$-${} $R_\oplus$)'.format(prad1,prad2) 
        dlogsinc = 0.5
        for i in range(3):
            smass1, smass2 = smass[i],smass[i+1]
            key = 'fitsinc_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            fmtkey = 'smass{}'.format(i+1)
            plot_sinc_rates(
                fit, log10(xlim_points[0]), log10(xlim_points[1])+0.1,
                dlogsinc, fmtkey
            )
            loglog()
        title(_title)
        
    sca(axL[0])
    setp(axL, xlabel='Incident Stellar Flux (Earth-units)', ylim=(1e-4,1e0))
    setp(
        axL[0],
        ylabel = 'Planets per Star'.format(dlogsinc),
        xlim =(1e4,1)
    )
    setp(axL[1], xlim =(1e4,1))
    for i in range(2):
        sca(axL[i])
        fig_label('cd'[i])

    tight_layout(True)

def fig_occur_per3():
    sns.set_context('paper',font_scale=1.1)
    fig, axL = subplots(ncols=3,figsize=(8,2.75),sharey=True)
    sizes = ['se','sn']
    for i in range(2):
        size = sizes[i]
        smass = [0.5, 0.7, 1.0, 1.4]

        if size=='se':
            prad1, prad2 = 1.0, 1.7
            _s = 'Super-Earths'
            xlim_points = [1,100]
        elif size=='sn':
            prad1, prad2 = 1.7, 4.0
            _s = 'Sub-Neptunes'
            xlim_points = [1,300]
            
        dlogper = 0.25
       
        for j in range(3):
            ax = axL[j] 
            sca(ax)
            smass1, smass2 = smass[j],smass[j+1]
            _title = '$M_\star$ = {}$-${} $M_\odot$'.format(smass1, smass2)
            title(_title)
            key = 'fitper_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            fmtkey = 'prad{}'.format(i+1)
            plot_per_rates(
                fit, log10(xlim_points[0]), log10(xlim_points[1]) + 0.1,
                dlogper,fmtkey
            )
            fit.res.flatchain['logx0']
            if j==0:
                text(0.6, 0.15 - i*0.07, _s, color=ptcolor[fmtkey], transform=ax.transAxes, size='small')

            lo,mi,hi = 10**fit.res.flatchain['logx0'].quantile([0.16,0.5,0.84])
            axvline(mi,color=ptcolor[fmtkey],alpha=1,zorder=0)
            axvspan(lo,hi,color=ptcolor[fmtkey],alpha=0.5,zorder=0)

    sca(axL[0])
    setp(axL, xlabel='Orbital Period (days)', ylim=(1e-4,1e0))
    setp(
        axL[0],
        ylabel = 'Planets per Star'.format(dlogper),
        xlim =(1,300),
    )
    setp(axL[1], xlim =(1,300))
    for i in range(3):
        sca(axL[i])
        loglog()
        fig_label('abc'[i])

    tight_layout(True)
    
def fig_occur_sinc3():
    sns.set_context('paper',font_scale=1.0)
    fig, axL = subplots(ncols=3,figsize=(8,2.75),sharey=True)
    sizes = ['se','sn']
    for i in range(2):
        size = sizes[i]
        sca(axL[i])
        smass = [0.5,0.7,1.0,1.4]
        if size=='se':
            prad1, prad2 = 1.0, 1.7
            _title = 'Super-Earths'
            xlim_points = [1,1e4]
        elif size=='sn':
            prad1, prad2 = 1.7, 4.0
            _title = 'Sub-Neptunes'
            xlim_points = [1,1e4]
            
        dlogsinc = 0.5
       
        for j in range(3):
            ax = axL[j] 
            sca(ax)
            smass1, smass2 = smass[j],smass[j+1]
            key = 'fitsinc_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            fmtkey = 'prad{}'.format(i+1)
            plot_sinc_rates(
                fit, log10(xlim_points[0]), log10(xlim_points[1]) + 0.1,
                dlogsinc,fmtkey
            )
            fit.res.flatchain['logx0']

            lo,mi,hi = 10**fit.res.flatchain['logx0'].quantile([0.16,0.5,0.84])
            axvline(mi,color=ptcolor[fmtkey],alpha=1,zorder=0)
            axvspan(lo,hi,color=ptcolor[fmtkey],alpha=0.5,zorder=0)
            
            
                
    sca(axL[0])
    setp(axL, xlabel='Incident Stellar Flux \n(Earth-units)', ylim=(1e-4,1e0))
    setp(
        axL[0],
        ylabel = 'Planets per Star'.format(dlogsinc),
        xlim =(1e4,1)
    )
    setp(axL,xlim =(1e4,1))
    for i in range(3):
        sca(axL[i])
        loglog()
        fig_label('def'[i])

    tight_layout(True)
    
from chainconsumer import ChainConsumer

def fig_occur_violin():
    sns.set_context('paper',font_scale=1.1)
    fig, axL = subplots(ncols=2,nrows=3,figsize=(6,5.5),sharex=True)

    smass = [0.5,0.7,1.0,1.4]
    prad = [1.0,1.7,4.0]
    for j in range(2):
        prad1, prad2 = prad[j],prad[j+1]
        data_f = []
        data_per0 = []
        data_sinc0 = []
        
        for i in range(3):
            smass1, smass2 = smass[i],smass[i+1]

            key = 'fitper_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            chain = fit.res.flatchain
            chain['x0'] = 10**chain['logx0']
            chain['f'] = 10**chain['logf']
            chain = chain[chain.x0.between(*chain.x0.quantile([0.01,0.99]))]
            #data.append(10**chain['logx0'])
            data_f.append(chain['logf'])
            chain['x0'] = chain['logx0']
            data_per0.append(chain['x0'])

            key = 'fitsinc_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            chain = fit.res.flatchain
            #chain['x0'] = 10**chain['logx0']
            chain['x0'] = chain['logx0']
            chain = chain[chain.x0.between(*chain.x0.quantile([0.01,0.99]))]
            data_sinc0.append(chain['x0'])

        col = 1-j    
        sca(axL[0,col])
        kw = dict(color='LightGray',scale='width',cut=0,inner='quartile')
        sns.violinplot(data=data_f, **kw)
        sca(axL[1,col])
        sns.violinplot(data=data_per0, **kw)
        sca(axL[2,col])
        sns.violinplot(data=data_sinc0, **kw)

    setp(axL[-1,:],xlabel=r'Stellar Mass (Solar-Masses)')
    xt=['0.5-0.7','0.7-1.0','1.0-1.4']
    setp(axL[0,:],xticklabels=xt)

    yt = [0.1,0.3,1,3]
    _yt = np.log10(yt)
    setp(axL[0,:],yticks=_yt)
    setp(axL[0,:],yticklabels=yt)
    setp(axL[0,1],title='Super-Earths ($R_p$ = 1.0-1.7 $R_\oplus$)')
    setp(axL[0,0],title='Sub-Neptunes ($R_p$ = 1.7-4.0 $R_\oplus$)')

    setp(axL[0,0],ylabel=r'$f_\star$ (planets per star)')
    yt = [1,3,10,30]
    _yt = np.log10(yt)
    setp(axL[1,0],ylabel=r'$P_0$ (days)')
    setp(axL[1,:],yticks=_yt)
    setp(axL[1,:],yticklabels=yt)

    setp(axL[2,0],ylabel=r'$S_\mathrm{inc,0}$ (Earth-units)')
    yt = [10,30,100,300,1000]
    _yt = np.log10(yt)
    setp(axL[2,:],yticks=_yt)
    setp(axL[2,:],yticklabels=yt)

    for i in range(6):
        sca(axL.T.flatten()[i])
        fig_label('abcdef'[i])
    tight_layout(True)

# Core plotting code


class SixPlotter(object):
    def __init__(self):
        self.mass1 = [0.5,0.7,1.0]
        self.mass2 = [0.7,1.0,1.4]
        
    def plot(self, plot_gradient):
        fig, axL = subplots(nrows=3,ncols=2,figsize=(6.5,6.5))
        i = 0

        #maxlevel = [6, 6 *(0.79 / 1.37), 6*(0.41 / 1.37)]
        #maxlevel = [6, 6, 6]
        maxlevel = [15, 3, 1]
        
        for _mass1, _mass2 in zip(self.mass1,self.mass2):
            key = '{}_smass={}-{}'.format(self.occur_prefix,_mass1,_mass2)
            occ = ckscool.io.load_object(key,cache=1)
            axd = axL[i,0]
            axo = axL[i,1]

            # Detected planets
            sca(axd)
            df = occ.plnt.copy()
            df = df.rename(
                columns={
                    'prad':'gdir_prad','per':'koi_period','sinc':'giso_sinc'
                }
            )
            pl = NDPlotter(df,self.xk,smass_lims=[_mass1,_mass2],zoom=True)
            pl.plot()

            gradkey = key.replace('occur','grad').replace('prad','prad-det')
            grad = ckscool.gradient.Gradient(gradkey)
            grad.load_csv()
            grad.plot_gradients('band')
            cbar = colorbar(pl.qc,shrink=0.8,format='%.1f')
            cbar.set_label(self.clabel,size='small')
            cbar.ax.tick_params(labelsize='x-small')
           
            
            # Occurrence rate 
            sca(axo)
            cp = pl.cp
            pl = self.ORDPlotter(occ, cp)
            #levels=linspace(0,maxlevel[i],14)
            #levels=logspace(np.log10(0.01),np.log10(maxlevel[i]),14)
            levels=None
            pl.plot_ord(levels=levels)
            

            pl.plot_completeness()
            gradkey = key.replace('occur','grad').replace('prad','prad-occ')
            grad = ckscool.gradient.Gradient(gradkey)
            grad.load_csv()
            grad.plot_gradients('band')
            title = '$M_\star = {}-{}\, M_\odot$ '.format(_mass1,_mass2)
            setp(axL[i,:],title=title)
            i+=1

        for ax in axL.flatten():
            sca(ax)
            xticks([log10(_xt) for _xt in self.xt],self.xt)
            yticks([log10(_yt) for _yt in self.yt],self.yt)
            grid()

        setp(axL, xlim=self.xlim, ylim=self.ylim)
        setp(axL, ylabel='')
        setp(axL, xlabel='')
        setp(axL[:,0], ylabel=self.ylabel)
        setp(axL[-1,:], xlabel=self.xlabel)
        tight_layout(True)

class SixPlotterPerPrad(SixPlotter):
    def __init__(self):
        super(SixPlotterPerPrad,self).__init__()
        self.xk = 'koi_period'
        self.xt = [1,3,10,30,100,300]
        self.yt = [1.0,1.4,2.0,2.8,4.0]
        self.xlim = log10(1),log10(300)
        self.ylim = log10(1),log10(4)
        self.ORDPlotter = ORDPlotterPerPrad
        self.xlabel='Orbital Period (days)'
        self.ylabel='Planet Size (Earth-radii)'
        self.clabel='$\mathrm{d}N/\mathrm{d}\log P/\mathrm{d}\log R_p$'
        
class SixPlotterSincPrad(SixPlotter):
    def __init__(self):
        super(SixPlotterSincPrad,self).__init__()
        self.xk = 'giso_sinc'
        self.xt = [3000,1000,300,100,30,10,3]
        self.yt = [1.0,1.4,2.0,2.8,4.0]
        self.xlim=log10(3000),log10(3)
        self.ylim=log10(1),log10(4)
        self.occur_prefix = 'occur-sinc-prad'
        self.ORDPlotter = ORDPlotterSincPrad
        self.xlabel='Incident Flux (Earth-units)'
        self.ylabel='Planet Size (Earth-radii)'
        self.clabel='$\mathrm{d} N / \mathrm{d} \log S_\mathrm{inc} / \mathrm{d} \log R_p$'

class ORDPlotter(object):
    """
    Occurrence rate density plotter
    """
    def __init__(self, occ, cp):
        """
        occ : occurrence object. can either be 
        """
        ds = cp.meshgrid()
        Z = occ.occurrence_rate_density_idem(array(ds.kxc), array(ds.kyc))
        Z = Z.reshape(ds.kxc.shape)
        ds['occrd'] = (['kx','ky'],Z)
        Z = occ.comp.prob_trdet_interp(ds.xc, ds.yc)
        ds['prob_trdet'] = (['kx','ky'],Z)
        ds['ntrial'] = occ.nstars * ds.prob_trdet
        self.ds = ds
        self.ntrials_min = 50

    def plot_ord(self,levels=None, gradient_array=False):
        """
        plot occurrence rate density
        """
        ds = self.ds
        if levels is None:
            eps = 1e-4
            maxz = ds.occrd.where(ds.ntrial > self.ntrials_min).max()
            maxz = np.round(maxz*1.1,3)
            levels = linspace(0,maxz+eps,14)
        if gradient_array:
            return ds.kxc.values, ds.kyc.values, ds.occrd.values

        cbarticks = levels[::2]
        cbarticklabels = ["{:.1f}".format(1e2*_yt) for _yt in cbarticks]
        cmap = 'YlGn' 
        kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)
        qc = contourf(ds.kxc,ds.kyc,ds.occrd,**kw)
        cbar = colorbar(qc,shrink=0.8,format='%.2f')
        cbar.set_label(self.cbarlabel,size='small')
        cbar.ax.tick_params(labelsize='x-small')

        
    def plot_completeness(self, gradient_array=False):
        """
        Gray out region of low completeness
        """
        ds = self.ds

        if gradient_array:
            return ds.ntrial.values, self.ntrials_min

        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(ds.kxc, ds.kyc, ds.ntrial, [0,self.ntrials_min], zorder=2.5,
                 cmap=cmap,vmax=1)

    def label(self):
        xlabel(self.xlabel)
        ylabel(self.ylabel)
        

class ORDPlotterPerPrad(ORDPlotter):
    """
    """
    def __init__(self, occ, cp):
        super(ORDPlotterPerPrad,self).__init__(occ,cp)
        self.xlabel = 'Orbital Period (days)'
        self.ylabel = 'Planet size (Earth-radii)'
        self.cbarlabel='$\mathrm{d}f/\mathrm{d}\log P/\mathrm{d}\log R_p$'

class ORDPlotterSincPrad(ORDPlotter):
    """
    """
    def __init__(self, occ, cp):
        super(ORDPlotterSincPrad,self).__init__(occ,cp)
        self.xlabel = 'Incident Flux (Earth-units)'
        self.ylabel = 'Planet size (Earth-radii)'
        self.cbarlabel ='$\mathrm{d} f / \mathrm{d} \log S_\mathrm{inc} / \mathrm{d} \log R_p$'


    
def plot_per_rates(fit, logper1, logper2, dlogper, fmtkey, plot_band=True, bandkw={}):
    """
    """
    logper = np.arange(logper1,logper2,dlogper)
    per = 10**logper
    df = dict(per1=per[:-1], per2=per[1:],prad1=fit.y1,prad2=fit.y2)
    df = pd.DataFrame(df)
    df['perc'] = np.sqrt(df.per1 * df.per2)

    rates = []
    for i, row in df.iterrows():
        rates += [fit.occ.occurrence_box(row)]

    rates = pd.DataFrame(rates)    
    rates = pd.concat([df,rates],ignore_index=False,axis=1)
    peri = np.logspace(np.log10(fit.x1),np.log10(fit.x2),300)
    ckscool.plot.occur.plot_rates('perc',rates,fmtkey)
    if hasattr(fit,'mi'):
        f = 10**fit.mi.params['logf']
        rate = f * fit.rate_lambda(peri, fit.mi.params)
        y = peri * rate # planets per e interval dN/dnper
        y = y * dlogper / np.log10(np.e) 
        _ = plot(peri,y,color=ptcolor[fmtkey])

    if plot_band:
        chain = fit.res.flatchain.sample(300)
        chain['f'] = 10**chain['logf']
        rate_lambda = fit.rate_lambda_sample(peri, chain) # nchain x nper
        rate = np.array(chain['f']).reshape(-1,1) * rate_lambda
        y = peri.reshape(1,-1) * rate # planets per e interval dN/dnper
        y = y * dlogper / np.log10(np.e) 
        lo,hi = np.percentile(y,[16,84],axis=0)
        fill_between(peri,lo,hi, color=ptcolor[fmtkey],alpha=0.3)


def plot_per_rates_ratio(fit1, fit2, logper1, logper2, dlogper, fmtkey, plot_band=True, bandkw={}):
    """
    """
    logper = np.arange(logper1,logper2,dlogper)
    per = 10**logper
    peri = np.logspace(np.log10(fit1.x1),np.log10(fit1.x2),300)

    yL = []
    for fit in [fit1, fit2]:
        chain = fit.res.flatchain.sample(300)
        chain['f'] = 10**chain['logf']
        rate_lambda = fit.rate_lambda_sample(peri, chain) # nchain x nper
        rate = np.array(chain['f']).reshape(-1,1) * rate_lambda
        y = peri.reshape(1,-1) * rate # planets per e interval dN/dnper
        y = y * dlogper / np.log10(np.e) 
        yL.append(y)

    y = yL[1] / yL[0]
    lo,hi = np.percentile(y,[16,84],axis=0)
    fill_between(peri,lo,hi, color=ptcolor[fmtkey],alpha=0.3)
        
def plot_sinc_rates(fit, logsinc1, logsinc2, dlogsinc, fmtkey, plot_band=True, bandkw={}):
    """
    """
    logsinc = np.arange(logsinc1,logsinc2,dlogsinc)
    sinc = 10**logsinc
    df = dict(sinc1=sinc[:-1], sinc2=sinc[1:],prad1=fit.y1,prad2=fit.y2)
    df = pd.DataFrame(df)
    df['sincc'] = np.sqrt(df.sinc1 * df.sinc2)

    rates = []
    for i, row in df.iterrows():
        rates += [fit.occ.occurrence_box(row)]

    rates = pd.DataFrame(rates)    
    rates = pd.concat([df,rates],ignore_index=False,axis=1)
    sinci = np.logspace(np.log10(fit.x1),np.log10(fit.x2),300)
    ckscool.plot.occur.plot_rates('sincc',rates,fmtkey)
    if hasattr(fit,'mi'):
        f = 10**fit.mi.params['logf']
        rate = f * fit.rate_lambda(sinci, fit.mi.params)
        y = sinci * rate # planets sinc e interval dN/dnsinc
        y = y * dlogsinc / np.log10(np.e) 
        _ = plot(sinci,y,color=ptcolor[fmtkey])

    if plot_band:
        chain = fit.res.flatchain.sample(300)
        chain['f'] = 10**chain['logf']
        rate_lambda = fit.rate_lambda_sample(sinci, chain) # nchain x nsinc
        rate = np.array(chain['f']).reshape(-1,1) * rate_lambda
        y = sinci.reshape(1,-1) * rate # planets sinc e interval dN/dnsinc
        y = y * dlogsinc / np.log10(np.e) 
        lo,hi = np.percentile(y,[16,84],axis=0)
        fill_between(sinci,lo,hi, color=ptcolor[fmtkey],alpha=0.3)

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
    x = occur[xk]
    y = occur.rate
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw1)
    errorbar(x,y*fac,yerr=yerr*fac, **ebkw2)
    occurul = occur.dropna(subset=['rate_ul'])
    if len(occurul) >0:
        plot(x,occur.rate_ul*fac,**ulkw)

ptcolor = {
    'smass1':'r',
    'smass2':'g',
    'smass3':'b',
    'prad1':'g',
    'prad2':'b',
}

bdcolor = {
    'smass1':'light red',
    'smass2':'light green',
    'smass3':'light blue',
    'prad1':'light green',
    'prad2':'light blue',
}

sns.set_style('ticks')
sns.set_color_codes()


