from __future__ import division
import time

import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from joblib import Parallel, delayed

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

    def line_integral(self, params, KDE, width=0):
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
        '''
        L = np.sqrt( (x1 - x0)**2 + (y1 - y0)**2 )
        dl = L / nbins
        points = np.vstack([xmid, ymid]).T
        integral = np.sum(KDE(points) * dl)
        nintegral = integral / L # normalized integral
        '''
        if width==0:
            points = np.vstack([xmid, ymid]).T
        if width > 0:
            nbinsy = 5
            dy = np.linspace(-0.5*width, 0.5*width, nbinsy)
            ymid = ymid.reshape(-1,1) + dy.reshape(1,-1)
            xmid = np.vstack([xmid] * nbinsy)
            points = np.array([xmid.flatten(),ymid.flatten()]).T
            
        nintegral = np.sum(KDE(points) )
        return nintegral

    def compute_gradient(self, x, y, z):
        """
        """
        # KDE interpolation
        map_interp = RegularGridInterpolator((x, y), z)

        def obj(params):
            return self.line_integral(params, map_interp)

        self.out = minimize(obj, self.params, method='Nelder-Meade')
        #self.out = minimize(obj, self.params, method='brute')
    
class GradientPerPrad(Gradient):
    def __init__(self, *args):
        super(GradientPerPrad, self).__init__(*args)
        params = Parameters()
        params.add('m', value=-0.01,min=-0.25,max=0.15)
        params.add('logR_10', value=np.log10(1.7),min=np.log10(1.4),max=np.log10(2.3))
        self.params = params
        self.x0 = 3
        self.x1 = 30
    
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
        params.add('logR_100', value=np.log10(1.7), min=np.log10(1.5), max=np.log10(2.5))
        self.params = params
        self.x0 = 10
        self.x1 = 1000
    
    def func(self,params, logS):
        """Returns the value of the radius versus flux line.
        """

        p = params.valuesdict()
        return p['m'] * logS + p['logR_100'] - 2 * p['m']

def load_gradient(objkey, seed=None, full_output=False):
    # constructs Gradient object, calculates chains and makes output
    # dataframe

    s1, s2 = objkey.split('_')
    name, xk, yk, mode = s1.split('-')
    smass1, smass2 = s2.replace('smass=','').split('-')
    smass1, smass2 = float(smass1),float(smass2)
    occkey = ( objkey.replace('grad','occur')
               .replace('-det','')
               .replace('-occ','') )
                
    occ = ckscool.io.load_object(occkey,cache=1,verbose=0)
    plnt = occ.plnt.copy()
    nplnt = len(plnt)
    if seed is not None:
        print("seed = {}".format(seed))
        plnt = plnt.sample(nplnt, replace=True, random_state=seed)
        
    if xk=='per':
        ndplotterkey = 'koi_period'
        Occurrence = ckscool.occur.OccurrencePerPrad
        grad = GradientPerPrad(objkey, [smass1, smass2])

    elif xk=='sinc':
        ndplotterkey = 'giso_sinc'
        Occurrence = ckscool.occur.OccurrenceSincPrad
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

    namemap = {'prad':'gdir_prad','per':'koi_period','sinc':'giso_sinc'}
    pl = ckscool.plot.planet.NDPlotter(
        plnt.rename(columns=namemap),ndplotterkey,zoom=False
    )
    if mode=='det':
        x, y, Z_det = pl.plot(gradient_array=True)
        grad.compute_gradient(x, y, Z_det)

        
    if mode=='occ':
        occ = Occurrence(plnt, occ.comp, occ.nstars)
        cp = pl.cp
        pl = ckscool.plot.occur.ORDPlotter(occ, cp)
        x_occ, y_occ, Z_occ = pl.plot_ord(gradient_array=True)
        Z_comp, ntrials_min = pl.plot_completeness(gradient_array=True)

        # set low completeness area to null occurrence
        Z_comp[Z_comp <= ntrials_min] = 0.0
        Z_comp[Z_comp >  ntrials_min] = 1.0
        Z_occ = Z_occ * Z_comp

        # find line of least occurrence
        grad.compute_gradient(x, y, Z_occ)

    if full_output:
        return grad, pl    

    return grad


def load_gradient_chain(key):
    niter = 128
    Ncores = 8
    func = lambda i : load_gradient(key, seed=i)
    obj = Parallel(n_jobs=Ncores)(delayed(func)(i) for i in range(niter))
    return obj
