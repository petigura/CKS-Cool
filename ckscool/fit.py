import numpy as np
from lmfit import Parameters 
import lmfit
import cksmet.model
import copy
import corner
from lmfit import Parameters
import numpy as np

fmt = {'kp':"{:.3f}",'logkp':"{:+.2f}",'beta':"{:+.1f}",'alpha':"{:+.1f}",'gamma':"{:.1f}",'gamma+beta':"{:.1f}",'per0':"{:.1f}"}

class Fit(object):
    def __init__(self, x, dx, nplnt, ntrial):
        """
        x (array): centers of bins (can be 1D or 2D)
        dx (array): volume of box
        nplnt (array): number of planets per bin
        ntrial (array): number of trials per bin
        """
        self.x = np.array(x)
        self.dx = dx
        self.nplnt = np.array(nplnt)
        self.ntrial = np.array(ntrial)

    def _loglike(self, params):
        """
        Args:
            params (lmfit.Parameters): parameters
        """
        # Planet distribution function evaluated at the center of each cell
        dfdx = self.model(params, self.x) 
        
        # Number of expected planet detections from model integrated
        # over each cell
        fcell = dfdx * self.dx
        fcell = fcell.reshape(*self.nplnt.shape)
        if (fcell > 1).any():
            print "fcell > 1" 
            return -1e10 # doesn't work for negative likelihood

        nnd = self.ntrial - self.nplnt # number of non detections (bin by bin)
        _loglike = np.sum(
            self.nplnt * np.log(fcell)) + np.sum(nnd * np.log(1 - fcell) 
        )
        if np.isnan(_loglike):
            return -1e10

        return _loglike

    def _negloglike(self,params):
        return -1.0 * self._loglike(params)

    def print_parameters(self,prefix=''):
        lines = self.to_string()
        print "\n".join(lines)
        
    def to_string(self,prefix=''):
        lines = []
        for k in self.pfit.keys():
            s = fmt[k]
            if self.pfit[k].vary:
                chain = self.flatchain[k]
                q = chain.quantile([0.16,0.50,0.84])
                val = s.format(q.ix[0.50])
                err1 = s.format(q.ix[0.84] - q.ix[0.50])
                err2 = s.format(q.ix[0.16] - q.ix[0.50])
                s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
                s = s.replace('++','+')
            else:
                s = "$%s$" % s.format(self.pfit[k].value)

            line = r"{{{}{}}}{{{}}}".format(prefix,k,s)
            lines.append(line)

        return lines
        
    def sample_chain(self, nsamp):
        psamples = self.flatchain.sample(nsamp)
        for k in self.p0.valuesdict().keys():
            if not self.p0[k].vary:
                psamples[k] = self.p0[k]
        return psamples

    def sample(self, x, nsamp,):
        # method only works for 1-D models        
        psamples = self.sample_chain(nsamp)
        
        model_samples = []
        for i, row in psamples.iterrows():
            psamp = copy.copy(self.p0)
            for k in self.p0.valuesdict().keys():
                psamp[k].value = row[k]

            model_samples.append( self.model(psamp, x) ) 

        model_samples = np.vstack(model_samples)
        return model_samples

    def fit(self):
        res = lmfit.minimize(self._negloglike, self.p0, method='Nelder')
        lmfit.report_fit(res)
        self.pfit = res.params
        self.loglike_fit = self._loglike(self.pfit)

    def mcmc(self, **kwargs):
        mini = lmfit.Minimizer(self._loglike, self.pfit, nan_policy='raise')
        res_emcee = mini.emcee(params=self.pfit, seed=1, **kwargs)
        self.flatchain = res_emcee.flatchain 
        
        print "loglike_fit = {:.4f}, max(loglike_mcmc) = {:.4f}".format(self.loglike_fit, np.max(res_emcee.lnprob) )
        if self.loglike_fit < np.max(res_emcee.lnprob):
            print ""
            print "WARNING maxlike from lmfit is smaller than maxlike from emcee"
            print ""
    def corner(self):
        corner.corner(self.flatchain)

class Exponential(Fit):
    def __init__(self, *args, **kwargs):
        super(Exponential, self).__init__(*args, **kwargs)
        self.model = smet_exp
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-10,max=10)
        p.add('beta',value=0.28,vary=True,min=-10,max=10)
        self.p0 = p

class PerPowerLawExpSmetExp(Fit):
    def __init__(self, *args, **kwargs):
        super(PerPowerLawExpSmetExp, self).__init__(*args, **kwargs)
        self.model = per_powerlaw_smet_exp
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-10,max=10)
        p.add('alpha',value=1.8,vary=False,min=-6,max=6)
        p.add('beta',value=0.28,vary=True,min=-6,max=6)
        self.p0 = p

class PowerLawCutoff(Fit):
    def __init__(self, *args, **kwargs):
        super(PowerLawCutoff, self).__init__(*args, **kwargs)
        self.model = powerlaw_cutoff
        p = Parameters()
        p.add('logkp',value=0,vary=True,min=-10,max=10)
        p.add('beta',value=0.0,vary=True)
        p.add('per0',value=7,vary=True, min=0,max=100)
        p.add('gamma',value=2,vary=True,min=0,max=4)
        self.p0 = p

    def to_string(self,prefix=''):
        # Prints the sum of gamma+beta
        lines = super(PowerLawCutoff,self).to_string(prefix=prefix)
        k =  'gamma+beta'
        chain = self.flatchain['gamma'] + self.flatchain['beta'] 
        s = fmt[k]
        q = chain.quantile([0.16,0.50,0.84])
        val = s.format(q.ix[0.50])
        err1 = s.format(q.ix[0.84] - q.ix[0.50])
        err2 = s.format(q.ix[0.16] - q.ix[0.50])
        s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
        s = s.replace('++','+')
        line = r"{{{}{}}}{{{}}}".format(prefix,k,s)
        lines.append(line)
        return lines

class CutoffPowerLaw(Fit):
    def __init__(self, *args, **kwargs):
        super(CutoffPowerLaw, self).__init__(*args, **kwargs)
        self.model = cutoff_powerlaw
        p = Parameters()
        p.add('logkp',value=0,vary=True,min=-10,max=10)
        p.add('beta',value=0.0,vary=False)
        p.add('logsinc0',value=2,vary=True, min=0,max=4)
        p.add('gamma',value=1.6,vary=True,min=0,max=4)
        self.p0 = p

class Gaussian(Fit):
    def __init__(self, *args, **kwargs):
        super(Gaussian, self).__init__(*args, **kwargs)
        self.model = gaussian
        p = Parameters()
        p.add('logkp',value=0,vary=True,min=-10,max=10)
        p.add('mu',value=2.0,vary=True)
        p.add('sig',value=2,vary=True, min=0,max=4)
        self.p0 = p

class LogGaussian(Fit):
    def __init__(self, *args, **kwargs):
        super(LogGaussian, self).__init__(*args, **kwargs)
        self.model = loggaussian
        p = Parameters()
        p.add('logkp',value=0,vary=True,min=-10,max=10)
        p.add('mu',value=.5,vary=True)
        p.add('sig',value=0.1,vary=True, min=0,max=1)
        self.p0 = p


class BrokenPowerLaw(Fit):
    def __init__(self, *args, **kwargs):
        super(BrokenPowerLaw, self).__init__(*args, **kwargs)
        self.model = broken_powerlaw
        p = Parameters()
        p.add('logf0',value=0,vary=True,min=-10,max=10)
        p.add(
            'logper0',value=np.log10(10),vary=True,
            min=np.log10(1),max=np.log10(100)
        )
        p.add('alpha1',value=2,vary=True,min=0,max=4)
        p.add('alpha2',value=0,vary=True,min=-3,max=3)
        self.p0 = p

def powerlaw_cutoff(params, x):
    """
    Equation 8 in Howard 2010
    """
    p = params.valuesdict()
    per = x
    occur = (
        10**p['logkp']  * 
        per**p['beta'] * 
        (1 - np.exp( -1.0 * (per / p['per0'])**p['gamma']))
     ) 
    return occur

def cutoff_powerlaw(params, x):
    """
    Equation 8 in Howard 2010
    """
    p = params.valuesdict()
    sinc = x
    sinc0 = 10**p['logsinc0']
    occur = (
        10**p['logkp']  * 
        sinc**p['beta'] * 
        (1 - np.exp( -1.0 * (sinc0 / sinc)**p['gamma']))
     ) 
    return occur

def broken_powerlaw(params, x):
    """
    2 comp powerlaw

    alpha1 : parameter below the break
    alpha2: parameter above the break
    per0: period of the break
    f0: normalization of the break

    """
    p = params.valuesdict()
    logper = np.log10(x)
    
    b1 = logper <= p['logper0']
    b2 = logper > p['logper0']

    logf = np.zeros_like(x)
    logf[b1] = p['alpha1'] * (logper[b1] - p['logper0']) + p['logf0']
    logf[b2] = p['alpha2'] * (logper[b2] - p['logper0']) + p['logf0']
    f = 10**logf 
    return f

def smet_exp(params, x):
    """
    df/dM = kp * 10**(beta * M)
    """
    smet = x
    p = params.valuesdict()
    occur = 10**p['logkp'] * 10**(p['beta'] * smet)
    return occur

def gaussian(params, x):
    """
    df / dx = kp * e ^ (-0.5 ((x - mu) / sig )**2)
    """
    p = params.valuesdict()
    per = x
    occur = 10**p['logkp']  * np.exp(-0.5 * ((x - p['mu'])/p['sig'])**2)
    return occur

def loggaussian(params, x):
    """
    df / dx = kp * e ^ (-0.5 ((x - mu) / sig )**2)
    """
    p = params.valuesdict()
    per = x

    kp = 10**p['logkp']  
    occur = kp  / x * np.exp(-0.5 * ((np.log(x) - p['mu'])/p['sig'])**2)
    return occur
    

