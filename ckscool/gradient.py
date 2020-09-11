from __future__ import division
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy import interpolate
from scipy import optimize
from numpy import random as rand
import corner
from numpy import logspace, log10, arange, sqrt
from joblib import Parallel, delayed
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import cPickle as pickle
import pandas as pd


import ckscool.comp
import ckscool.io
import ckscool.grid
import ckscool.occur
import ckscool.plot.occur
import ckscool.gradient
from ckscool.plot.planet import NDPlotter
from ckscool.plot.occur import ORDPlotter


class Gradient(object):

    def __init__(self, objkey, smass_lims):
        self.objkey = objkey
        self.smass_lims = smass_lims

    def resampled_gradient(self, plnt_full, comp, nstars, iteration):

        # first calculate detection gradient
        # resample with replacement (seeded)
        plnt = plnt_full.sample(len(plnt_full),replace=True, random_state=iteration)
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
            pl = NDPlotter(df,'koi_period',zoom=False)
        elif self.objkey=='grad-sinc-prad':
            pl = NDPlotter(df,'giso_sinc',zoom=False)
        x, y, Z_det = pl.plot(gradient_array=True)

        # find line of least detection
        sol_det = self.find_gradient(x, y, Z_det)


        # then calculate occurrence gradient
        cp = pl.cp
        pl = ORDPlotter(occ, cp)
        x_occ, y_occ, Z_occ = pl.plot_ord(gradient_array=True)
        Z_comp, ntrials_min = pl.plot_completeness(gradient_array=True)

        # set low completeness area to null occurrence
        Z_comp[Z_comp <= ntrials_min] = 0.0
        Z_comp[Z_comp >  ntrials_min] = 1.0
        Z_occ = Z_occ * Z_comp

        # find line of least occurrence
        sol_occ = self.find_gradient(x, y, Z_occ)

        print('{0} complete...'.format(iteration))

        return [sol_det[0], sol_det[1], sol_occ[0], sol_occ[1]]


    def gradient_chain(self, N_cores, N_iter):

        [smass1, smass2] = self.smass_lims

        # load planet sample
        field = ckscool.io.load_table('field-cuts',cache=1)
        field = field[~field.isany]
        field = field.rename(columns={'ber19_srad':'srad','ber19_smass':'smass'})
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
            comp = ckscool.io.load_object('comp-per-prad_smass={0}-{1}'.format(smass1,smass2),cache=1)
        elif self.objkey=='grad-sinc-prad':
            comp = ckscool.io.load_object('comp-sinc-prad_smass={0}-{1}'.format(smass1,smass2),cache=1)

        # compute gradient chain (parallelised)
        self.grad_chain = Parallel(n_jobs=N_cores)(delayed(self.resampled_gradient)(plnt, comp, nstars, i) for i in np.arange(N_iter))
        self.grad_chain = np.array(self.grad_chain)


    def gradient_file(self):

        """
        Writes gradient values (to be quoted in paper) to pandas dataframe
        """

        grad_chain = self.grad_chain

        m_det_bounds = np.percentile(grad_chain[:,0], [16,50,84])
        m_det_err = [m_det_bounds[2]-m_det_bounds[1], m_det_bounds[1]-m_det_bounds[0]]

        m_occ_bounds = np.percentile(grad_chain[:,2], [16,50,84])
        m_occ_err = [m_occ_bounds[2]-m_occ_bounds[1], m_occ_bounds[1]-m_occ_bounds[0]]

        R0_det_bounds = np.percentile(grad_chain[:,1], [16,50,84])
        R0_det_err = [10**R0_det_bounds[2]-10**R0_det_bounds[1], 10**R0_det_bounds[1]-10**R0_det_bounds[0]]

        R0_occ_bounds = np.percentile(grad_chain[:,3], [16,50,84])
        R0_occ_err = [10**R0_occ_bounds[2]-10**R0_occ_bounds[1], 10**R0_occ_bounds[1]-10**R0_occ_bounds[0]]

        if self.objkey=='grad-per-prad':
            data = {'m'        :  [m_det_bounds[1], m_occ_bounds[1]],
                    'm_err1'   :  [m_det_err[0], m_occ_err[0]],
                    'm_err2'   :  [-m_det_err[1], -m_occ_err[1]],
                    'Rp10'     :  [10**R0_det_bounds[1], 10**R0_occ_bounds[1]],
                    'Rp10_err1':  [R0_det_err[0], R0_occ_err[0]],
                    'Rp10_err2':  [-R0_det_err[1], -R0_occ_err[1]],
                    'map'      :  ['Detection', 'Occurrence']
                    }

        elif self.objkey=='grad-sinc-prad':
            data = {'m'         :  [m_det_bounds[1], m_occ_bounds[1]],
                    'm_err1'    :  [m_det_err[0], m_occ_err[0]],
                    'm_err2'    :  [-m_det_err[1], -m_occ_err[1]],
                    'Rp10'      :  [R0_det_bounds[1], R0_occ_bounds[1]],
                    'Rp100_err1':  [10**R0_det_err[0], 10**R0_occ_err[0]],
                    'Rp100_err2':  [-R0_det_err[1], -R0_occ_err[1]],
                    'map'       :  ['Detection', 'Occurrence']
                    }


        self.output_data = pd.DataFrame(data=data)

            



class GradientPerPrad(Gradient):


    def prad_per(self,logP,m,logR_10):
        """
        Returns the value of the radius versus period line.

        Arguments:
        logP    : log period (days) value
        m       : gradient
        logR_10 : log intercept at 10 day period

        Returns:
        log R: log radius (R_earths)
        """

        return m * logP + logR_10 - m


    def line_integral_per(self, theta, KDE):

        """
        Return the line integral for a given straight line


        Arguments:
            theta:       array of [m,logR_10] where m is gradient and logR_10 is
                        intercept at 10 days
            occurence:   interpolated occurence / detection map, must be from the following -
                        scipy.interpolate.RectBivariateSpline(S_array, R_array, occ)

        Returns:
            Total line integral from P[min] to P[max] with contributions coming from
            the occurence map.
        """

        m, logR_10 = theta
        P_domain = np.linspace(log10(3),log10(80),30)

        
        logR_10_lim = [log10(1.3), log10(3.0)]
        m_lim = [-0.15,0.15]

        integral = 0.0

        if  logR_10_lim[0] <= logR_10 <= logR_10_lim[1] \
        and m_lim[0] <= m <= m_lim[1]:

            for i in range(len(P_domain)-1):

                L = sqrt( (P_domain[-1] - P_domain[0])**2 + (self.prad_per(P_domain[-1],m,logR_10) - self.prad_per(P_domain[0],m,logR_10))**2)
                dl = L / len(P_domain)

                R_mid = ( self.prad_per(P_domain[i+1],m,logR_10) + self.prad_per(P_domain[i],m,logR_10) ) / 2
                P_mid = ( P_domain[i+1] + P_domain[i] ) / 2

                contribution = KDE([P_mid, R_mid]) * dl

                if contribution <= 1e-5: # spline sometimes <= 0 for masked region
                    contribution = 0.0
                    integral += contribution
                else:
                    integral += contribution

            return integral
        else:
            return 1e10


    def find_gradient(self, logP_array, logR_array, KDE):
        """
        Runs a minimisation algorithm to find the gradient and intercept that
        minimise the line integral through the occurence map.

        Arguments:
            P_array:        Array of per values used to create meshgrid
            R_array:        Array of radii values used to create meshgrid
            occurence:      Occurence defined on meshgrid of np.meshgrid(P_array, R_array)

        Returns:
            Solution to minimisation problem:  [m, logR_10]
        """
        # KDE interpolation
        map_interp = RegularGridInterpolator((logP_array, logR_array), KDE)

        sol = optimize.minimize(self.line_integral_per, x0=[0.01,log10(1.7)], args=(map_interp),
                                method='Nelder-Mead')
        return sol.x





class GradientSincPrad(Gradient):


    def prad_sinc(self,logS,m,logR_100):
        """
        Returns the value of the radius versus sinc line.

        Arguments:
        logS     : log bolometric flux (Earth units) value
        m        : gradient
        logR_100 : log intercept at 100 x Earth incident bolometric flux

        Returns:
        log R: log radius (R_earths)
        """

        return m * logS + logR_100 - 2*m




    def line_integral_sinc(self, theta, KDE):
        """
        Return the line integral for a given straight line


        Arguments:
            theta:       array of [m,logR_100] where m is gradient and logR_100 is
                         intercept at 100 x earth incident flux
            occurence:   interpolated occurence / detection map, must be from the following -
                         scipy.interpolate.RectBivariateSpline(S_array, R_array, occ)

        Returns:
            Total line integral from S[min] to S[max] with contributions coming from
            the occurence map.
        """

        m, logR_100 = theta
        S_domain = np.linspace(log10(10),log10(1000),30)

        logR_100_lim = [log10(1.3), log10(3.0)]
        m_lim = [-0.15,0.15]

        integral = 0.0

        if  logR_100_lim[0] <= logR_100 <= logR_100_lim[1] \
        and m_lim[0] <= m <= m_lim[1]:

            for i in range(len(S_domain)-1):

                L = sqrt( (S_domain[-1] - S_domain[0])**2 + (self.prad_sinc(S_domain[-1],m,logR_100) - self.prad_sinc(S_domain[0],m,logR_100))**2)
                dl = L / len(S_domain)

                R_mid = ( self.prad_sinc(S_domain[i+1],m,logR_100) + self.prad_sinc(S_domain[i],m,logR_100) ) / 2
                S_mid = ( S_domain[i+1] + S_domain[i] ) / 2

                contribution = KDE([S_mid, R_mid]) * dl

                if contribution <= 1e-5: # spline sometimes <= 0 for masked region
                    contribution = 0.0
                    integral += contribution
                else:
                    integral += contribution

            return integral
        else:
            return 1e10


    def find_gradient(self, logS_array, logR_array, KDE):
        """
        Runs a minimisation algorithm to find the gradient and intercept that
        minimise the line integral through the occurence map.

        Arguments:
            S_array:        Array of sinc values used to create meshgrid
            R_array:        Array of radii values used to create meshgrid
            occurence:      Occurence defined on meshgrid of np.meshgrid(P_array, R_array)

        Returns:
            Solution to minimisation problem:  [m, logS_100]
        """
        # KDE interpolation
        map_interp = RegularGridInterpolator((logS_array, logR_array), KDE)

        sol = optimize.minimize(self.line_integral_sinc, x0=[0.01,log10(1.7)], args=(map_interp),
                                method='Nelder-Mead')
        return sol.x




def construct_grad(objkey, limits, N_cores=4, N_sample=10000):

    # constructs Gradient object, calculates chains 
    # and makes output dataframe

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']

    if objkey=='grad-per-prad':
        Gradient = GradientPerPrad(objkey, [smass1, smass2])
    if objkey=='grad-sinc-prad':
        Gradient = GradientSincPrad(objkey, [smass1, smass2])

    Gradient.gradient_chain(N_cores, N_sample)
    Gradient.gradient_file()

    return Gradient

