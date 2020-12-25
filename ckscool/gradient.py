from __future__ import division
import time

import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from lmfit import minimize, Parameters
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import ckscool.comp
import ckscool.io
import ckscool.occur
import ckscool.plot.occur

class Gradient(object):
    """
    Gradient 

    base class to compute gradients
    """
    
    def __init__(self, objkey, smass_lims):
        self.objkey = objkey
        self.smass_lims = smass_lims

    def resampled_gradient(self, plnt_full, comp, nstars, iteration):

        # first calculate detection gradient
        # resample with replacement (seeded)
        nplnt = len(plnt_full)
        plnt = plnt_full.sample(nplnt, replace=True, random_state=iteration)
        if self.objkey=='grad-per-prad':
            Occurrence = ckscool.occur.OccurrencePerPrad
        elif self.objkey=='grad-sinc-prad':
            Occurrence = ckscool.occur.OccurrenceSincPrad
        occ = Occurrence(plnt, comp, nstars)
        df = occ.plnt.copy()
        df = df.rename(
            columns={
                'prad':'gdir_prad','per':'koi_period','sinc':'giso_sinc'
            }
        )
        # produce arrays for gradient calculations
        if self.objkey=='grad-per-prad':
            pl = ckscool.plot.planet.NDPlotter(df,'koi_period',zoom=False)
        elif self.objkey=='grad-sinc-prad':
            pl = ckscool.plot.planet.NDPlotter(df,'giso_sinc',zoom=False)
        x, y, Z_det = pl.plot(gradient_array=True)

        # find line of least detection
        sol_det = self.find_gradient(x, y, Z_det)

        # then calculate occurrence gradient
        cp = pl.cp
        pl = ckscool.plot.occur.ORDPlotter(occ, cp)
        x_occ, y_occ, Z_occ = pl.plot_ord(gradient_array=True)
        Z_comp, ntrials_min = pl.plot_completeness(gradient_array=True)

        # set low completeness area to null occurrence
        Z_comp[Z_comp <= ntrials_min] = 0.0
        Z_comp[Z_comp >  ntrials_min] = 1.0
        Z_occ = Z_occ * Z_comp

        # find line of least occurrence
        sol_occ = self.find_gradient(x, y, Z_occ)
        print('{0} complete...'.format(iteration))

        return [sol_det, sol_occ]

    def gradient_chain(self, N_cores, N_iter):

        [smass1, smass2] = self.smass_lims

        # load planet sample
        field = ckscool.io.load_table('field-cuts',cache=1)
        field = field[~field.isany]
        field = field.rename(columns={'ber20_srad':'srad','ber20_smass':'smass'})
        plnt = ckscool.io.load_table('planets-cuts2')
        plnt = plnt[~plnt.isany]

        namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass',
                    'giso_sinc':'sinc'}
        plnt = plnt.rename(columns=namemap)
        field = field[field.smass.between(smass1,smass2)]
        plnt = plnt[plnt.smass.between(smass1,smass2)]
        field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
        nstars = len(field)

        # load completeness object
        if self.objkey=='grad-per-prad':
            key = 'comp-per-prad_smass={0}-{1}'.format(smass1,smass2)
        elif self.objkey=='grad-sinc-prad':
            key = 'comp-sinc-prad_smass={0}-{1}'.format(smass1,smass2)
        comp = ckscool.io.load_object(key,cache=1)

        # compute gradient chain (parallelised)
        self.grad_chain = Parallel(n_jobs=N_cores)(delayed(self.resampled_gradient)(plnt, comp, nstars, i) for i in np.arange(N_iter))

    def line_integral(self, params, KDE):
        """Return the line integral for a given line

        Arguments:

            params: array of parameters that will be passed to self.func

            KDE: interpolated occurence / detection map, must be from
                        the following -
                        scipy.interpolate.RectBivariateSpline(S_array,
                        R_array, occ)

        Returns:
            Total line integral from P[min] to P[max] with contributions coming from
            the occurence map.

        """

        nbins = 30
        x0, x1 = np.log10(self.x0), np.log10(self.x1)
        xbins = np.linspace(x0, x1, nbins+1)
        xmid = 0.5 * ( xbins[1:] + xbins[:-1] )
        y0 = self.func(params, x0)
        y1 = self.func(params, x1)
        ymid = self.func(params, xmid)
        L = np.sqrt( (x1 - x0)**2 + (y1 - y0)**2 )
        dl = L / nbins
        points = np.vstack([xmid, ymid]).T
        integral = np.sum(KDE(points) * dl)
        nintegral = integral / L # normalized integral
        return nintegral

    def find_gradient(self, x, y, z):
        """
        """
        # KDE interpolation
        map_interp = RegularGridInterpolator((x, y), z)

        def obj(params):
            return self.line_integral(params, map_interp)

        out = minimize(obj, self.params, method='Nelder-Meade')
        return out
    
class GradientPerPrad(Gradient):
    def __init__(self, *args):
        super(GradientPerPrad, self).__init__(*args)
        params = Parameters()
        params.add('m', value=-0.01,min=-0.15,max=0.15)
        params.add('logR_10', value=np.log10(1.7),min=np.log10(1.4),max=np.log10(2.3))
        self.params = params
        self.x0 = 1
        self.x1 = 100
        
    
    def func(self,params, logP):
        """Returns the value of the radius versus period line.
        """

        p = params.valuesdict()
        return p['m'] * logP + p['logR_10'] - p['m']

class GradientSincPrad(Gradient):
    def __init__(self, *args):
        super(GradientSincPrad, self).__init__(*args)
        params = Parameters()
        params.add('m', value=0.05, min=-0.05, max=0.2)
        params.add('logR_100', value=np.log10(2.0), min=np.log10(1.5), max=np.log10(2.5))
        self.params = params
        self.x0 = 10
        self.x1 = 1000
    
    def func(self,params, logS):
        """Returns the value of the radius versus flux line.
        """

        p = params.valuesdict()
        return p['m'] * logS + p['logR_100'] - 2 * p['m']

def construct_grad(objkey, limits, N_cores=4, N_sample=10000):
    # constructs Gradient object, calculates chains and makes output
    # dataframe

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']

    if objkey=='grad-per-prad':
        grad = GradientPerPrad(objkey, [smass1, smass2])
    if objkey=='grad-sinc-prad':
        grad = GradientSincPrad(objkey, [smass1, smass2])
        if [smass1, smass2] == [0.5,0.7]:
            x0 = 3
            x1 = 100
        elif [smass1, smass2] == [0.7,1.0]:
            x0 = 10
            x1 = 300
        elif [smass1, smass2] == [1.0,1.4]:
            x0 = 30
            x1 = 1000
        elif [smass1, smass2] == [0.5,1.4]:
            x0 = 3
            x1 = 1000
        else:
            assert False, "mass limits not supported"

        grad.x0 = x0 
        grad.x1 = x1
        
    grad.gradient_chain(N_cores, N_sample)
    return grad

