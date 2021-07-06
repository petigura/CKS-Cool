from __future__ import division
import time

import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from joblib import Parallel, delayed
import os

import ckscool.comp
import ckscool.io
import ckscool.occur
import ckscool.plot.occur

import pandas as pd
import copy 
from matplotlib import pylab as plt
from numpy import log10

class Gradient(object):

    def __init__(self, objkey):
        s1, s2 = objkey.split('_')
        self.name, self.xk, self.yk, self.mode = s1.split('-')
        smass1, smass2 = s2.replace('smass=','').split('-')
        smass1, smass2 = float(smass1),float(smass2)
        self.occkey = ( objkey.replace('grad','occur')
                   .replace('-det','')
                   .replace('-occ','') )

        self.outdir = os.path.join(ckscool.io.ANALYSISDIR, objkey)
        os.system('mkdir -p {}'.format(self.outdir)) 
        self.objkey = objkey
        self.smass1 = smass1
        self.smass2 = smass2

        fn = self.objkey+'.csv'
        fn = os.path.join(self.outdir,fn)
        self.csvfn = fn

        if self.xk=='per':
            self.ndplotterkey = 'koi_period'
            self.xlim = (0.5,1.5) # compute minima only over a certain range
            self.ylim = (0.15,0.35)
            self.x0 = 1 # log of anchor point
            self.ORDPlotter = ckscool.plot.occur.ORDPlotterPerPrad
            self.Occurrence2D = ckscool.occur.OccurrencePerPrad
            
        if self.xk=='sinc':
            self.ndplotterkey = 'giso_sinc'
            if smass1==0.5 and smass2==0.7:
                self.xlim = (log10(3),log10(100)) 
            if smass1==0.7 and smass2==1.0:
                self.xlim = (log10(10),log10(300)) 
            if smass1==1.0 and smass2==1.4:
                self.xlim = (log10(30),log10(1000)) 
            if smass1==0.5 and smass2==1.4:
                self.xlim = (log10(10),log10(1000)) 
            
            self.ORDPlotter = ckscool.plot.occur.ORDPlotterSincPrad

            self.ylim = (0.15,0.35)
            self.x0 = 2

        if self.xk=='smass':
            self.ndplotterkey = 'giso_smass'
            self.xlim = (np.log10(0.6),np.log10(1.4))
            self.ylim = (0.15,0.35)
            self.x0 = 0

            #ndplotkey = ['koi_period','giso_sinc','giso_smass','cks_smet']

        if self.xk=='smet':
            self.ndplotterkey = 'cks_smet'
            self.xlim = (-0.3,0.3) # compute minima only over a certain range
            self.ylim = (0.15,0.35)
            self.x0 = 0

        self.xi = np.linspace(self.xlim[0],self.xlim[1],100)

    def sample(self, nsamples, nplots=20):
        seeds = [None] + range(nsamples)
        fits = []
        i = 0
        plnt = ckscool.io.load_table('planets-cuts2',cache=1)
        plnt = plnt[~plnt.isany]
        plnt = plnt[plnt.giso_smass.between(self.smass1,self.smass2)]
        plnt0 = plnt.copy()
        nplnt = len(plnt)

        if self.mode=='occ':
           occ0 = ckscool.io.load_object(self.occkey,cache=1)

        for seed in seeds:
            outfile = self.objkey+'_seed={}.pdf'.format(seed)
            outfile = os.path.join(self.outdir,outfile)
            if seed is not None:
                print("seed = {}".format(seed))
                plnt = plnt0.sample(nplnt, replace=True, random_state=seed)

            # NDPlotter needed for both detections and occurrence
            # gradient calculations. 
            ndplotter = ckscool.plot.planet.NDPlotter(
                plnt,self.ndplotterkey,zoom=True)
            
            x, y, Z = ndplotter.plot(gradient_array=True)

            if self.mode=='occ':
                occ = copy.deepcopy(occ0)
                namemap = {
                    'koi_period':'per',
                    'giso_sinc':'sinc',
                    'gdir_prad':'prad'
                }
                occ = self.Occurrence2D(
                    plnt.rename(columns=namemap),occ0.comp, occ0.nstars
                )
                ordplotter = self.ORDPlotter(occ, ndplotter.cp)

                # we need xy from ndplotter, but we overwrite Z with
                # the gradient array
                _, _, Z = ordplotter.plot_ord(gradient_array=True)

            Z = Z.T # need to transpose array here

            # generate a 2D meshgrids with all rows except top and
            # bottom and all columns.
            irows = np.arange(1, Z.shape[0] - 1) 
            irows2, icols2 = np.meshgrid(
                irows, range(Z.shape[1]), indexing='ij'
            )

            # Identfiy all points lower than the row directly above and below 
            b = ( (Z[irows,:] < Z[irows+1,:])
                  & (Z[irows,:] < Z[irows-1,:]) )


            df = pd.DataFrame(dict(row=irows2[b],col=icols2[b]))
            df['x'] = x[df.col]
            df['y'] = y[df.row]
            df['z'] = Z[df.row,df.col]

            # only consider the minima winthin the prespecified range
            df = df[df.x.between(*self.xlim) & df.y.between(*self.ylim)]
            dfmin  = df.sort_values(by=['x','z']) \
                    .groupby('x',as_index=False).first()
            
            pfit = np.polyfit(dfmin.x - self.x0, dfmin.y, 1)
            fits.append(dict(seed=seed,m=pfit[0],y0=pfit[1]) )

            i += 1
            if i > nplots:
                continue
            
            plt.figure(figsize=(6,4))

            if self.mode=='det':
                ndplotter.plot()

            if self.mode=='occ':
                ordplotter.plot_ord()
            
            plt.plot(df.x,df.y,'.b')
            plt.plot(dfmin.x,dfmin.y,'.r')
            plt.plot(self.xi,np.polyval(pfit,self.xi-self.x0))
            plt.tight_layout()
            plt.gcf().savefig(outfile)
                
        fits = pd.DataFrame(fits)
        fits.to_csv(self.csvfn)

    def load_csv(self):
        self.fits = pd.read_csv(self.csvfn)

    def plot_gradients(self, mode, xi=None):
        if xi==None:
            xi = self.xi

        yi = []
        for i, row in self.fits.iterrows():
            pfit = np.array([row.m, row.y0])
            yi += [np.polyval(pfit,xi-self.x0)]
        yi = np.vstack(yi)

        if mode=='band':
            lo, hi = np.percentile(yi,[16,84],axis=0)
            plt.fill_between(xi, lo, hi, color='b', alpha=0.4)
            return 

        if mode=='samples':
            plt.plot(xi, yi.T, 'r')
            return

        assert False, "mode = {} not supported ".format(mode)
            
        
'''    
class Gradient(object):
    """
    Gradient 

    base class to compute gradients

    Implements the following methods
        - line-integral which minimizes the line integral
        - chisq which minimizes chisq along the lowest minimum between two specified bounds. 

    """
    
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

    
    def compute_gradient(self, x, y, z, method):
        """
        x: grid of x values
        y: grid of y values
        z: grid of z values
        method: fitting method
        """

        if method=='line-integral':
            # KDE interpolation
            map_interp = RegularGridInterpolator((x, y), z)

            def obj(params):
                return self.line_integral(params, map_interp)

            self.out = minimize(obj, self.params, method='Nelder-Meade')

        if method=='valley-chisq':
            # have to transpose the array to make indexing work
            Z = z.T

            nrows = Z.shape[0]
            irows = np.arange(1,nrows-1)
            b = (Z[irows,:] < Z[irows+1,:]) & (Z[irows,:] < Z[irows-1,:])

            irows2,icols2 = np.meshgrid(irows,range(Z.shape[1]),indexing='ij')
            lmink_ = pd.DataFrame(dict(row=irows2[b],col=icols2[b]))
            lmin['x'] = x[lmin.col]
            lmin['y'] = y[lmin.row]
            lmin['z'] = Z[lmin.row,lmin.col]

            # only consider the minima winthin the prespecified range
            lmin = lmin.query('0.5 < x  < 1.5 and 0.15 < y <  0.35') 
            lmin  = lmin.sort_values(by=['x','z']).groupby('x',as_index=False).first()

            def obj(params):
                ymod = self.func(params, lmin.x)
                return np.sum((lmin.y - ymod)**2)

            self.out = minimize(obj, self.params, method='Nelder-Meade')
            
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

'''    
def load_z(objkey, seed=None, full_output=False):
    """
    Loads x, y, z arrays that are used in gradient calculations
    """
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

    elif xk=='sinc':
        ndplotterkey = 'giso_sinc'
        Occurrence = ckscool.occur.OccurrenceSincPrad
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


    namemap = {'prad':'gdir_prad','per':'koi_period','sinc':'giso_sinc'}
    pl = ckscool.plot.planet.NDPlotter(
        plnt.rename(columns=namemap),ndplotterkey,zoom=False
    )
    if mode=='det':
        x, y, Z_det = pl.plot(gradient_array=True)
        z = Z_det
        
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
        z = Z_occ
        
    return x, y, z

def load_gradient(objkey, seed=None, full_output=False):
    x, y, z = load_z(objkey, seed=seed)
    s1, s2 = objkey.split('_')
    name, xk, yk, mode = s1.split('-')
    if xk=='per':
        grad = GradientPerPrad()
    grad.compute_gradient(x,y,z,'valley-chisq')
    return grad
    
def load_gradient_chain(key):
    niter = 100
    Ncores = 1
    func = lambda i : load_gradient(key, seed=i)
    obj = Parallel(n_jobs=Ncores)(delayed(func)(i) for i in range(niter))
    return obj
