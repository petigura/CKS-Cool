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
from ckscool.comp import CompletenessPerPrad, CompletenessSincPrad

class Occurrence2D(object):
    """
    Base class for occurrence object. Holds generic functionality for
    occurrence objects

    """
    def __init__(self, plnt, comp, nstars):
        self.plnt = plnt
        self.comp = comp
        self.nstars = nstars
        self.plntx = array(plnt[self.xk])
        self.plnty = array(plnt[self.yk])

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


    def planet_weights(self):
        return 1 / self.comp.prob_trdet_interp(self.plntx,self.plnty)

    def occurrence_rate_density_idem(self, logx, logy):
        """
        Occurrence rate density (ORD) at logx logy

        logx : log10 of the period at which to compute ord
        logy : log10 of the planet radius at which to compute ord
        """
        assert logx.shape ==logy.shape
        w = self.planet_weights()
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
    def __init__(self,*args):
        super(OccurrencePerPrad,self).__init__(*args)
        assert isinstance(self.comp,CompletenessPerPrad), \
            "comp must be instance of CompletenessPerPrad"

class OccurrenceSincPrad(Occurrence2D):
    xk = 'sinc'
    yk = 'prad'
    xbw = log10(1 + 1)
    ybw = log10(1 + 0.05)
    def __init__(self,*args):
        super(OccurrenceSincPrad,self).__init__(*args)
        assert isinstance(self.comp,CompletenessSincPrad), \
            "comp must be instance of CompletenessSincPrad"


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


class MeanPlanetSize(object):
    def __init__(self,query):
        smassbins = np.array([0.5,0.7,1.0,1.4])
        self.smass1 = smassbins[:-1]
        self.smass2 = smassbins[1:]
        self.smassc = np.sqrt(self.smass1 * self.smass2)
        self.nsamp = 1000
        self.nbins = len(self.smass1)
        self.w = []
        self.plnt = []
        for i in range(len(self.smass1)):
            smass1, smass2 = self.smass1[i],self.smass2[i]
            k = 'occur-per-prad_smass={}-{}'.format(smass1,smass2)
            occ = ckscool.io.load_object(k,cache=1)
            plnt = occ.plnt.query(query)
            w = 1/ occ.comp.prob_trdet_interp(plnt['per'],plnt['prad'])
            self.plnt.append(plnt)
            self.w.append(w)

    def sample_logprad_single(self):
        """
        sample mean planet size
        """
        mn = []
        for i in range(self.nbins):
            w = self.w[i]
            plnt = self.plnt[i]
            logprad = log10(plnt['prad'])
            val = logprad.sample(len(w),replace=True)
            mn.append(np.average(val,weights=w))

        return mn

    def sample_logprad(self):
        logprad = []
        for i in range(self.nsamp):
            logprad.append(self.sample_logprad_single())

        self.logprad = np.array(logprad)

    def fit_powerlaw(self):
        logprad = np.array(self.logprad)
        logsmassc = log10(self.smassc)
        self.logsmassci = np.linspace(log10(0.5),log10(1.4),100)
        self.mn = np.mean((10**logprad),axis=0)
        self.std = np.std((10**logprad),axis=0)
        self.coeff = np.polynomial.polynomial.polyfit(
            logsmassc,logprad.T,1
        )
        self.logpradf = np.polynomial.polynomial.polyval(
            self.logsmassci,self.coeff
        )

        self.pradf = 10**self.logpradf
        self.q16,self.q84 = np.percentile(self.pradf,[16,84],axis=0)
    
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
