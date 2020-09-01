"""
grid
"""
import numpy as np
from numpy import log10, logspace, arange, round, array
from scipy.special import gammaln as gamln
from scipy import special
from sklearn.utils import resample
import pandas as pd 

import calc
import ckscool.io
import ckscool.grid
#import ckscool.gradient

class Occurrence2D(object):
    """
    Base class for occurrence object. Holds generic functionality for
    occurrence objects

    """
    def __init__(self, plnt, comp, nstars):
        self.plnt = plnt
        self.comp = comp
        self.nstars = nstars

    def occurence_box(self, limits):
        """Compute occurrence in a little box

        We make the assumption that the dN/dlogP and dN/lopRp is constant
        within a box.

        Args:
            limits (dict): must contain, per1, per2, prad1, prad2
        """
        
        prad1 = limits['prad1']
        prad2 = limits['prad2']
        per1 = limits['per1']
        per2 = limits['per2']

        # Get planet sample and number of stars
        cut = self.plnt.copy()
        cut = cut[cut.prad.between(prad1,prad2)]
        cut = cut[cut.per.between(per1,per2)]
        nplnt = len(cut)
        prob_trdet_mean, prob_det_mean = self.comp.mean_prob_trdet(
            per1, per2, prad1, prad2
        )
        ntrial = self.nstars * prob_trdet_mean
        rate = nplnt / ntrial

        nsample = int(1e4)
        binom = Binomial(ntrial, nplnt)
        samples = binom.sample(nsample)

        uplim = nplnt==0
        rate = samples_to_rate(samples,uplim=uplim)
        out = {}
        out['ntrial'] = ntrial
        out['nplnt'] = nplnt
        out['prob_trdet_mean'] = prob_trdet_mean
        out['prob_det_mean'] = prob_det_mean
        out = dict(out,**rate)
        return out

    def occurrence_grid(self, per1=None, per2=None, dlogper=None, 
                            prad1=None, prad2=None, dlogprad=None):
        """
        A convenience function to compute occurrence over a grid in period and radius
        """
        eps = 1e-9

        if dlogper is None:
            per = [per1,per2]
        else:
            per = 10**np.arange(log10(per1),log10(per2)+eps,dlogper)
            
        if dlogprad is None:
            prad = [prad1,prad2]
        else:
            prad = 10**np.arange(log10(prad1),log10(prad2)+eps,dlogprad)
        
        df = []
        for i in range(len(per)-1):
            _per1 = per[i]
            _per2 = per[i+1]
            for j in range(len(prad)-1):
                _prad1 = prad[j]
                _prad2 = prad[j+1]
                limits = dict(
                    per1=_per1, per2=_per2, prad1=_prad1, prad2=_prad2
                )
                _out = self.occurence_box(limits)
                _out['perc'] = np.sqrt(_per1*_per2)
                _out['pradc'] = np.sqrt(_prad1*_prad2)
                df.append(_out)
        df = pd.DataFrame(df)
        return df

    def occurrence_rate_density_idem(self, logx, logy):
        """
        Occurrence rate density (ORD) at logx logy

        logx : log10 of the period at which to compute ord
        logy : log10 of the planet radius at which to compute ord
        """
        assert logx.shape ==logy.shape
        plnt = self.plnt.copy()
        plntx = array(plnt[self.xk])
        plnty = array(plnt[self.yk])
        w = 1 / self.comp.prob_trdet_interp(plntx,plnty)
        nplnt = len(plnt)
        logxi = np.log10(plntx)
        logyi = np.log10(plnty)
        occrd = gaussian_2d_kde(
            logx, logy, logxi, logyi, self.xbw, self.ybw ,w=w
        )
        occrd /= self.nstars
        occrd = occrd.reshape(logx.shape)
        return occrd

class OccurrencePerPrad(Occurrence2D):
    """Occurrence

    Args:
        plnt (pandas DataFrame): table with planet parameters must contain
            following keys:
            - prad: planet radius
            - per: planet period
        comp (ckscool.comp.Completeness): see docs
        nstars (int): number of stars in parent stellar population

    """
    xk = 'per'
    yk = 'prad'
    xbw = log10(1 + 1)
    ybw = log10(1 + 0.05)

class OccurrenceSincPrad(Occurrence2D):
    xk = 'sinc'
    yk = 'prad'
    xbw = log10(1 + 1)
    ybw = log10(1 + 0.05)

def gaussian_2d_kde(x, y, xi, yi, xbw, ybw, w=None):
    """
    Convenience function to compute gaussian KDE
    """
    pos = np.vstack([x.flatten(),y.flatten()]).T
    ni = len(xi)
    mu = np.vstack([xi, yi]).T
    cov = np.array(ni * [np.eye(2)])
    cov[:,0,0] *= xbw**2
    cov[:,1,1] *= ybw**2
    return gaussian(pos, mu, cov, w=w) 

class Binomial(object):
    """Class that computes binomial statistics

    Args:
        n (float): number of trials
        k (float): number of sucesses
    """
    def __init__(self, n, k):
        self.n = float(n)
        self.k = float(k)
        self.npdf = 10000
        if k!=0:
            # Detection, calculate pdf of the rate of planets in box
            rate = 1.0 * k / n
            self.maxrate = min(1,10*rate)
        else:
            # Non-detection, upper limit on rate.
            self.maxrate = min(1, 10 / self.n)

    def pdf(self):
        """
        Returns probabilty of rate r
        """
        _rate = np.linspace(0,self.maxrate,self.npdf)
        _pdf = binom_gamma(_rate,self.n,self.k)
        _pdf /= _pdf.sum()
        return _rate, _pdf

    def sample(self, nsamp):
        rate, pdf = self.pdf()
        cdf = pdf.cumsum()
        fp = rate
        x = np.random.random_sample(nsamp)
        return np.interp(x, cdf, rate)

    def hist_samples(self, nsamp, downsamp=10):
        """
        Verify that the sampling working
        """

        rate, pdf = self.pdf()
        samples = self.sample(nsamp)
        weights = ones_like(samples)/nsamp # normalize histogram
        hist(samples,bins=rate[::downsamp],weights=weights)
        plot(rate[::downsamp],pdf[::downsamp]*downsamp)

def binom_gamma(p, n, k):
    """
    Probability that rate is p given there are n trials and k sucesses
    """
    n = 1.0*n
    k = 1.0*k
    p = 1.0*p
    combiln = (gamln(n+1) - (gamln(k+1) + gamln(n - k + 1)))
    _logpdf = combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)
    _pdf = np.exp(_logpdf)
    return _pdf

def sum_cells(ntrial, nplnt):
    """
    Add the occurrence values from different cells.

    Args:
        nplnt (arary): number of planets detected
    """

    ntrial = np.array(ntrial)
    nplnt = np.array(nplnt)
    nsample = int(1e4)

    samplesL = []
    np.random.seed(0)
    for i in range(len(nplnt)):
        binom = Binomial(ntrial[i], nplnt[i])
        samples = binom.sample(nsample)
        samplesL.append(samples)

    samplesL = np.vstack(samplesL)
    isuplim = (nplnt==0) # True if cell yields upper limit

    samples_sum = samplesL[~isuplim].sum(0)
    d = {}
    d = dict(rate=None, rate_err1=None, rate_err2=None, rate_ul=None)

    # All measurements are upper limits
    if (nplnt==0).all():
        samples = samplesL.sum(0)
        d = samples_to_rate(samples,uplim=True)
    else:
        samples = samplesL[~isuplim].sum(0)
        d = samples_to_rate(samples,uplim=False)
    return d

def samples_to_rate(samples, uplim=False):
    d = dict(rate=None, rate_err1=None, rate_err2=None, rate_ul=None)
    if uplim:
        p16, p50, p84, p90 = np.percentile(samples, [16,50,84,90])
        d['rate_ul'] = p90
        d['rate_str'] = "< {rate_ul:.4f} (90%)".format(**d)
    else:
        p16, p50, p84, p90 = np.percentile(samples, [16,50,84,90])
        d['rate'] = p50
        d['rate_err1'] = p84 - p50
        d['rate_err2'] = p16 - p50
        d['rate_str'] = "{rate:.4f} +/- {rate_err1:.4f}/{rate_err2:.4f}".format(**d)
    return d

def load_comp3d():
    """
    Iterate over a bunch of 2D completeness 
    """

def load_occur(objkey, limits, debug=False, sinc=False):
    """
    Constructs occurrence object
    """

    # Derive completeness object
    method = 'fulton-gamma-clip' # treatment for planet detectability
    impact = 0.8 # maximum impact parameter considered.

    field = ckscool.io.load_table('field-cuts',cache=1)
    field = field[~field.isany]
    field = field.rename(columns={'ber19_srad':'srad','ber19_smass':'smass'})
    plnt = ckscool.io.load_table('planets-cuts2')
    plnt = plnt[~plnt.isany]

    namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass',
               'giso_sinc':'sinc'}
    plnt = plnt.rename(columns=namemap)

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']
        field = field[field.smass.between(smass1,smass2)]
        plnt = plnt[plnt.smass.between(smass1,smass2)]

    n1 = len(field)
    field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
    nstars = len(field)
    
    print "{}/{} stars remain after droping nulls ".format(nstars,n1)
    if objkey=='occur-per-prad':
        xbins = np.round(logspace(log10(0.1),log10(1000),65),4)
        ybins = np.round(logspace(log10(0.25),log10(64),51 ),2)
        xk = 'per'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessPerPrad
        Occurrence = ckscool.occur.OccurrencePerPrad

    elif objkey=='occur-sinc-prad':
        xbins = np.round(logspace(log10(0.1),log10(100000),65),4)
        ybins = np.round(logspace(log10(0.25),log10(64),51 ),2)
        xk = 'sinc'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessSincPrad
        Occurrence = ckscool.occur.OccurrenceSincPrad

    else:
        assert False, "{} not supported objkey".format(objkey)

    debug=1
    if debug:
        xbins = xbins[::2]
        ybins = ybins[::2]

    comp_bins_dict = {xk: xbins,yk: ybins}
    spacing_dict = {xk:xscale, yk:yscale}
    grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)
    comp = Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=True)
    comp.compute_grid_prob_tr(verbose=True)
    comp.create_splines()
    occ = Occurrence(plnt, comp, nstars)
    return occ

def load_occur_resample(smass1, smass2, plnt, comp, nstars):
    resample_indices = resample(np.arange(len(plnt)))
    plnt_resampled =  plnt.iloc[resample_indices, :]
    occ = ckscool.occur.Occurrence(plnt_resampled, comp, nstars)
    return occ


def gaussian(pos, mu, cov, w=None):

    """
    Compute the sum of 2D gaussians

    last dimension of x must equal len(mu)
    """
    assert mu.shape[0]==cov.shape[0]
    assert mu.shape[1]==cov.shape[1]

    if type(w) is type(None):
        w = np.ones(mu.shape[0])

    nres = pos.shape[0] # number of resulting points
    ndim = pos.shape[1] # number of dimensions
    npts = mu.shape[0] # number of points used for KDE
    
    invcov = np.linalg.inv(cov) # (npts x ndim)
    detcov = np.linalg.det(cov) # (npts)
    norm = (2 * np.pi)**(-0.5 * ndim) * (detcov**-0.5) # npts
    if type(w) != type(None):
        norm *= w

    res = np.zeros(nres) # provision results array 
    # Iterate over points
    for i in range(npts):
        _invcov = invcov[i]
        _mu = mu[i]
        _norm = norm[i]
        resid = _mu - pos
        temp = np.matmul(_invcov,resid.reshape(-1,ndim,1))
        temp = np.matmul(resid.reshape(-1,1,ndim),temp)
        dist = temp.reshape(-1)
        res += _norm * np.exp(-0.5 * dist )

    return res
