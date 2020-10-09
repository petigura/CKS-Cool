from lmfit import Parameters, Parameter
import lmfit
import cksmet.model
import copy
import numpy as np
import scipy.integrate
import scipy.special
import scipy.interpolate
import pandas as pd

class Fit(object):
    def __init__(self,occ, per1, per2, prad1, prad2):
        self.occ = occ
        self.prad1 = prad1
        self.prad2 = prad2
        self.per1 = per1
        self.per2 = per2
        b = ( occ.plnt.per.between(per1,per2)
              & occ.plnt.prad.between(prad1,prad2))
        self.plnt = occ.plnt[b]

    def create_completeness_spline(self):
        """
        Evaluate average completeness as a function of period
        """

        x = np.linspace(self.per1,self.per2,100)
        dx = 0.1
        prob_trdet_mean = []
        for _x in x:
            temp, _ = self.occ.comp.mean_prob_trdet(
                _x - dx, _x + dx , self.prad1, self.prad2
            )
            prob_trdet_mean.append(temp)
        prob_trdet_mean = np.array(prob_trdet_mean)
        self.comp = scipy.interpolate.UnivariateSpline(x,prob_trdet_mean,s=0)

    def rate_Lambda(self, p):
        """
        Lambda from Bryson+20
        """
        v = p.valuesdict()
        nknots = 300
        f = 10**v['logf']
        x = np.logspace(np.log10(self.per1),np.log10(self.per2), nknots)
        etalambda = self.comp(x) * self.rate_lambda(x, p) 
        product = scipy.interpolate.UnivariateSpline(x, etalambda, s=0)
        integral =  product.integral(self.per1, self.per2)
        return integral * f * self.occ.nstars

    def logprob(self,p):
        v = p.valuesdict()
        f = 10**v['logf']
        _Lambda = self.rate_Lambda(p)
        _lambda = self.rate_lambda(self.plnt.per, p)
        _sumloglambda = np.sum(np.log(f * _lambda))
        _logprob = - _Lambda + _sumloglambda
        return _logprob

    def fit_max_likelihood(self, p, **kwargs):
        mi = lmfit.minimize(lambda x : -self.logprob(x), p, **kwargs)
        mi.params.pretty_print()
        self.mi = mi

    def run_mcmc(self,short=False):
        burn = 1000
        steps = 5000 
        if short:
            burn = 100
            steps = 500 
        thin = 1 
        nwalkers = 8
        res = lmfit.minimize(
            self.logprob, method='emcee', nan_policy='omit', burn=burn,
            steps=steps, thin=thin, nwalkers=nwalkers, params=self.mi.params,
            progress=True, float_behavior='posterior'
        )
        res.params.pretty_print()
        self.res = res

    def rate_lambda_sample(self, per, chain):
        """
        sample the normalize occurrence rate density from mcmc
        
        Returns:
            array: shape = (len(chain), len(nper))
        """
        rate = []
        for i,row in chain.iterrows():
            params = Parameters()
            for k,param0 in self.res.params.iteritems():
                if param0.vary:
                    params[k] = Parameter(k,value=row[k])
                else:
                    params[k] = Parameter(k,value=param0.value)
            rate.append( self.rate_lambda(per,params) )
        rate = np.vstack(rate)
        return rate

class FitSmoothBrokenPowerLaw(Fit):
    def IndefInt(self,k1,k2,u,x0):
        k1p1 = k1+1
        dk =  k1 - k2
        out = ( x0 * u**k1p1/k1p1
                * scipy.special.hyp2f1(1, k1p1/dk, k1p1/dk + 1, -u**dk) )
        return out

    def DefInt(self,k1,k2,x0,x1,x2):
        return self.IndefInt(k1,k2,x2/x0,x0) - self.IndefInt(k1,k2,x1/x0,x0)

    def rate_lambda(self, per, p):
        """
        must integrate to 1
        """
        v = p.valuesdict()
        per0 = 10**v['logper0']
        u = per / per0
        denom = u**(-v['k1']) + u**(-v['k2'])
        norm = self.DefInt(v['k1'],v['k2'],per0,self.per1,self.per2)
        _lambda_per = 1.0 / denom / norm
        return _lambda_per

class FitBrokenPowerLaw(Fit):
    def InDefInt1(self, per, per0, k1):
        return per**(k1 + 1) / (k1 + 1)

    def InDefInt2(self, per, per0, k1, k2):
        return per0**(k1 - k2) / (k2 + 1) * per**(k2 + 1)

    def rate_lambda(self, per, p):
        """
        must integrate to 1
        """
        if type(per) is float:
            per = np.array([per])
        v = p.valuesdict()
        per0 = 10**v['logper0']
        _lambda_per = np.zeros_like(per)
        b = per <= per0 
        _lambda_per[b] = per[b]**v['k1']
        b = per > per0 
        _lambda_per[b] = per0**(v['k1'] - v['k2']) * per[b]**v['k2']

        norm = ( self.InDefInt1(per0, per0, v['k1'])
                 - self.InDefInt1(self.per1, per0, v['k1'])
                 + self.InDefInt2(self.per2, per0, v['k1'], v['k2'])
                 - self.InDefInt2(per0, per0, v['k1'], v['k2']) ) 
        
        _lambda_per = _lambda_per / norm
        return _lambda_per
