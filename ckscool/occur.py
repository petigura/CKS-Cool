"""
Defines occurrence object
"""
import numpy as np
from numpy import log10, logspace, arange, round
from scipy.special import gammaln as gamln
from scipy import special
from statsmodels.distributions.empirical_distribution import ECDF

import ckscool.io
import ckscool.grid

class Occurrence(object):
    def __init__(self, plnt, comp, nstars, smet_field=None):
        self.plnt = plnt
        self.comp = comp
        self.nstars = nstars
        self.plnt = plnt

    def occurence_box(self, limits):
        """Compute occurrence in a little box

        We make the assumption that the dN/dlogP and dN/lopRp is constant
        within a box.

        Args:
            limits (dict): must contain, per1, per2, prad1, prad2, can
                optionally contain, smet1, smet2
        """
        out = dict()
        prad1 = limits['prad1']
        prad2 = limits['prad2']
        per1 = limits['per1']
        per2 = limits['per2']

        # Get planet sample and number of stars
        cut = self.plnt.copy()
        cut = cut[cut.prad.between(prad1,prad2)]
        cut = cut[cut.per.between(per1,per2)]
        nstars = self.nstars
        nplnt = len(cut)
        prob_trdet_mean, prob_det_mean = self.comp.mean_prob_trdet(
            per1, per2, prad1, prad2
        )
        ntrial = nstars * prob_trdet_mean
        rate = nplnt / ntrial

        nsample = int(1e4)
        binom = Binomial(ntrial, nplnt)
        samples = binom.sample(nsample) 

        uplim = nplnt==0
        rate = samples_to_rate(samples,uplim=uplim)
        out['ntrial'] = ntrial
        out['nplnt'] = nplnt
        out['prob_trdet_mean'] = prob_trdet_mean
        out['prob_det_mean'] = prob_det_mean
        out = dict(out,**rate)
        return out

class Binomial(object):
    def __init__(self, n, k):
        """
        Args:
            n: number of trials
            k: number of sucesses
        """
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

def load_occur(limits, debug=False):
    """
    Constructs occurrence object
    """
    
    # Derive completeness object
    method = 'fulton-gamma-clip' # treatment for planet detectability
    impact = 1 # maximum impact parameter considered.

    field = ckscool.io.load_table('field-cuts',cache=1)
    field = field[~field.isany]
    field = field.rename(columns={'ber18_srad':'srad','m17_smass':'smass'})

    plnt = ckscool.io.load_table('planets-cuts2+iso')
    plnt = plnt[~plnt.isany]
    namemap = {'gdir_prad':'prad','koi_period':'per','m17_smass':'smass'}
    plnt = plnt.rename(columns=namemap)

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']
        field = field[field.smass.between(smass1,smass2)]
        plnt = plnt[plnt.smass.between(smass1,smass2)]

    elif limits.has_key('bmr1'):
        bmr1 = limits['bmr1']  
        bmr2 = limits['bmr2']
        xs = 'gaia2_bpmag - gaia2_rpmag'
        field = field[field.eval(xs).between(bmr1,bmr2)]
        plnt = plnt[plnt.eval(xs).between(bmr1,bmr2)]

    n1 = len(field)
    field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
    n2 = len(field)
    print "{}/{} stars remain after droping nulls ".format(n2,n1)

    # Define grid of period and radius to compute completeness
    comp_per_bins = np.round(logspace(log10(0.1),log10(1000),65),4)
    comp_prad_bins = np.round(logspace(log10(0.25),log10(64),51 ),2)
    
    # debugging
    if debug:
        comp_per_bins = comp_per_bins[:6]
        comp_prad_bins = comp_prad_bins[:6]
        
    comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)

    comp = ckscool.comp.Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=False)
    comp.compute_grid_prob_tr(verbose=False)
    comp.create_splines()
    
    nstars = len(field)
    occ = ckscool.occur.Occurrence(plnt, comp, nstars)
    return occ


