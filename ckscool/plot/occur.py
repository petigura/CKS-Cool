import lmfit
from matplotlib.pylab import *
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
import pandas as pd
import seaborn as sns
from scipy import ndimage as nd

import ckscool.io
import ckscool.plot.planet
import ckscool.gradient #import R, R_sinc
from .planet import NDPlotter

sns.set_style('ticks')
sns.set_color_codes()

def fig_contour_six_per(gradient=False):
    sns.set_context('paper')
    pl = SixPlotterPerPrad()
    pl.plot()

def fig_contour_six_sinc(gradient=False):
    sns.set_context('paper')
    pl = SixPlotterSincPrad()
    pl.plot()

class SixPlotter(object):
    def __init__(self):
        self.mass1 = [0.5,0.7,1.0]
        self.mass2 = [0.7,1.0,1.4]
        
    def plot(self):
        fig, axL = subplots(nrows=3,ncols=2,figsize=(8.5,9))
        i = 0
        for _mass1, _mass2 in zip(self.mass1,self.mass2):
            key = '{}_smass={}-{}'.format(self.occur_prefix,_mass1,_mass2)
            occ = ckscool.io.load_object(key,cache=1)
            axd = axL[i,0]
            axo = axL[i,1]

            # Detected planets
            sca(axd)
            df = occ.plnt.copy()
            df = df.rename(columns={'prad':'gdir_prad','per':'koi_period','sinc':'giso_sinc'})
            pl = NDPlotter(df,self.xk,zoom=True)
            pl.plot()
            cbar = colorbar(pl.qc,shrink=0.5,format='%.1f')
            cbar.set_label('$dN/d \log P/ d \log R_p$',size='x-small')
            cbar.ax.tick_params(labelsize='xx-small')

            # Occurrence rate 
            sca(axo)
            cp = pl.cp
            pl = self.ORDPlotter(occ, cp)
            pl.plot_ord()
            pl.plot_completeness()
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
        self.occur_prefix = 'occur-per-prad'
        self.ORDPlotter = ORDPlotterPerPrad
        self.xlabel='Orbital Period (days)'
        self.ylabel='Planet Size (Earth-radii)'
        
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

    def plot_ord(self,levels=None):
        """
        plot occurrence rate density
        """
        ds = self.ds
        if levels==None:
            eps = 1e-4
            maxz = ds.occrd.where(ds.ntrial > self.ntrials_min).max()
            maxz = np.round(maxz*1.1,3)
            levels = linspace(0,maxz+eps,14)

        cbarticks = levels[::2]
        cbarticklabels = ["{:.1f}".format(1e2*_yt) for _yt in cbarticks]
        cmap = 'YlGn' 
        kw = dict(levels=levels,extend='neither',cmap=cmap,zorder=0)
        qc = contourf(ds.kxc,ds.kyc,ds.occrd,**kw)
        cbar = colorbar(qc,shrink=0.5,format='%.1f')
        cbar.set_label(self.cbarlabel,size='x-small')
        cbar.ax.tick_params(labelsize='xx-small')

    def plot_completeness(self):
        """
        Gray out region of low completeness
        """
        ds = self.ds
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
        self.cbarlabel = '$df/ d \log P / d \log R_p$'

class ORDPlotterSincPrad(ORDPlotter):
    """
    """
    def __init__(self, occ, cp):
        super(ORDPlotterSincPrad,self).__init__(occ,cp)
        self.xlabel = 'Incident Flux (Earth-units)'
        self.ylabel = 'Planet size (Earth-radii)'
        self.cbarlabel = '$df/ d \log Sinc / d \log R_p$'

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
        #self.plot_best()

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
