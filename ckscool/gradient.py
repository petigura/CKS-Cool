from __future__ import division
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy import interpolate
from scipy import optimize
from numpy import random as rand
import corner
from numpy import logspace, log10, arange
from joblib import Parallel, delayed
from scipy import interpolate


import ckscool.comp
import ckscool.io
import ckscool.grid
import ckscool.occur
import ckscool.plot.occur
import ckscool.gradient


from sklearn.utils import resample


# ---------------------------------------------------------------------------- #

def R(P,m,Rp_10):
    """
    Returns the value of the radius versus period line.
    As this straight line is in log-space, it takes the
    form of a power law in linear space.

    Arguments:
    P     : period value (days)
    m     : gradient
    Rp_10 : intercept at 10 day period

    Returns:
    R: Radius (M_earths)
    """

    return Rp_10 * (P/10)**m



def R_sinc(F,m,Rp_100):
    """
    Returns the value of the radius versus period line.
    As this straight line is in log-space, it takes the
    form of a power law in linear space.

    Arguments:
    P     : period value (days)
    m     : gradient
    Rp_10 : intercept at 10 day period

    Returns:
    R: Radius (M_earths)
    """

    return Rp_100 * (F/100)**m



# ---------------------------------------------------------------------------- #

def line_integral(theta, map, maptype, sinc=False):
    """
    Return the line integral for a given straight line


    Arguments:
        theta:       array of [m,c] where m is gradient and c is
                     intercept
        occurence:   interpolated occurence map, must be from the following -
                     scipy.interpolate.RectBivariateSpline(P_array, R_array, occ)

    Returns:
        Total line integral from P[min] to P[max] with contributions coming from
        the occurence map.
    """

    m, Rp_10 = theta
    if sinc:
        if maptype == 'occurrence':
            P_domain = np.logspace(np.log10(1),np.log10(80),100)
        elif maptype == 'detections':
            P_domain = np.logspace(np.log10(1),np.log10(300),100)

        Rp10_lim = [1.5, 2.3]
        m_lim = [-0.1,0.1]
    else:
        if maptype == 'occurrence':
            P_domain = np.logspace(np.log10(1),np.log10(200),100)
        elif maptype == 'detections':
            P_domain = np.logspace(np.log10(1),np.log10(200),100)
        Rp10_lim = [1.7, 2.4]
        m_lim = [-0.2,0.2]

    # constraints for intercepts with P(R_min) and P(R_max)
    if m_lim[0] <= m <= m_lim[1] and Rp10_lim[0] <= Rp_10 <= Rp10_lim[1]:

        # sum contributions to integrand
        integral = 0
        for i in range(len(P_domain)-1):

            if sinc:
                dl = np.sqrt((P_domain[i+1]-P_domain[i])**2 + (R_sinc(P_domain[i],m,Rp_10)-R_sinc(P_domain[i+1],m,Rp_10))**2)
                R_mid = ( R_sinc(P_domain[i+1],m,Rp_10) + R_sinc(P_domain[i],m,Rp_10) ) / 2

            else:
                dl = np.sqrt((P_domain[i+1]-P_domain[i])**2 + (R(P_domain[i],m,Rp_10)-R(P_domain[i+1],m,Rp_10))**2)
                R_mid = ( R(P_domain[i+1],m,Rp_10) + R(P_domain[i],m,Rp_10) ) / 2

            P_mid = ( P_domain[i+1] + P_domain[i] ) / 2
            contribution = map(P_mid,R_mid)[0,0] * dl

            if contribution <= 1e-5: # spline sometimes <= 0 for masked region
                contribution = 0.0
                integral += contribution
            else:
                integral += contribution



        return integral
    else:
        return 1e10

# ---------------------------------------------------------------------------- #

def gradient(P_array, R_array, map, maptype, sinc=False):
    """
    Runs a minimisation algorithm to find the gradient and intercept that
    minimise the line integral through the occurence map.

    Arguments:
        P_array:        Array of period values used to create meshgrid
        R_array:        Array of radii values used to create meshgrid
        occurence:      Occurence defined on meshgrid of np.meshgrid(P_array, R_array)

    Returns:
        Solution to minimisation problem:  [m, Rp_10]
    """
    # KDE interpolation
    map_interp = interpolate.RectBivariateSpline(P_array,R_array,map)

    if sinc:
        sol = optimize.minimize(line_integral, x0=[0.0,1.6], args=(map_interp, maptype, sinc),
                                method='Nelder-Mead')
    else:
        sol = optimize.minimize(line_integral, x0=[-0.08,1.8], args=(map_interp, maptype, sinc),
                                method='Nelder-Mead')

    return sol

# ---------------------------------------------------------------------------- #

def bootstrap_occurrence(smass1, smass2, plnt, comp, nstars, x_range, y_range, seed, sinc=False, debug=False):

    """
    resamples planet population (with replacement) and finds best fit gradient
    of Fulton gap.

    Args:
        smass1, smass2     : lower and upper stellar mass for sample
        plnt               : entire sample of planets
        comp               : completeness object
        nstars             : number of stars in field after cuts
    """

    # parallelised so need to reimport on each cpu
    import ckscool.occur
    import ckscool.plot.occur
    import ckscool.gradient

    # seed random number generator
    np.random.seed(seed)
    if sinc:
        ntrials_min = 100
    else:
        ntrials_min = 100

    if sinc:
        # resample and calculate occurence
        occ_i = ckscool.occur.load_occur_resample(smass1, smass2, plnt, comp, nstars, sinc=True)
        cp_i = ckscool.plot.occur.load_contour_plotter_sinc(occ_i)
        # get occurence meshgrid and values
        X, Y, Z, ntrial = ckscool.plot.occur.gradient_arrays_sinc(cp_i)
        if debug:
            ckscool.plot.occur.contour_sinc(cp_i, ntrials_min=ntrials_min)
    else:
        # resample and calculate occurence
        occ_i = ckscool.occur.load_occur_resample(smass1, smass2, plnt, comp, nstars)
        cp_i = ckscool.plot.occur.load_contour_plotter(occ_i)
        # get occurence meshgrid and values
        X, Y, Z, ntrial = ckscool.plot.occur.gradient_arrays(cp_i)
        if debug:
            ckscool.plot.occur.contour(cp_i, ntrials_min=ntrials_min)


    # mask occurence with low ntrial
    ntrial[ntrial <= ntrials_min] = 0.0
    ntrial[ntrial >  ntrials_min] = 1.0
    Z = Z * ntrial

    # find gradient and intercept
    if sinc:
        sol_i = ckscool.gradient.gradient(x_range, y_range, Z, 'occurrence', sinc=True)
        if debug:
            sinc_domain = np.log10(np.logspace(0,3,100))
            R_sinc = log10([ckscool.gradient.R_sinc(10**i,sol_i.x[0], sol_i.x[1]) for i in sinc_domain])
            plt.plot(sinc_domain, R_sinc,label=r'm={0:.2f}, R$_p$(100)={1:.2f}'.format(sol_i.x[0], sol_i.x[1]))
            plt.legend()
            plt.savefig('./test_{}'.format(seed))
            plt.clf()
            plt.close()
    else:
        sol_i = ckscool.gradient.gradient(x_range, y_range, Z, 'occurrence', sinc=False)
        if debug:
            per_domain = np.log10(np.logspace(0,2,100))
            R_per = log10([ckscool.gradient.R(10**i,sol_i.x[0], sol_i.x[1]) for i in per_domain])
            plt.plot(per_domain, R_per,label=r'm={0:.2f}, R$_p$(10)={1:.2f}'.format(sol_i.x[0], sol_i.x[1]))
            plt.legend()
            plt.savefig('./test_{}'.format(seed))
            plt.clf()
            plt.close()
    return list(sol_i.x)

# ---------------------------------------------------------------------------- #

def bootstrap_detection(plnt, seed, sinc=False, debug=False):

    """
    resamples planet population (with replacement) and finds best fit gradient
    of Fulton gap.

    Args:
        smass1, smass2     : lower and upper stellar mass for sample
        plnt               : entire sample of planets
        comp               : completeness object
        nstars             : number of stars in field after cuts
    """

    # parallelised so need to reimport on each cpu
    import ckscool.occur
    import ckscool.plot.occur
    import ckscool.gradient

    # seed random number generator
    np.random.seed(seed)

    resample_indices = resample(np.arange(len(plnt)))
    plnt =  plnt.iloc[resample_indices, :]

    if sinc:
        X,Y,Z = ckscool.plot.planet._sinc_prad(plnt, for_gradient=True)
        sinc_range = np.array([10**i for i in X])
        prad_range = np.array([10**i for i in Y])
        sol_i = ckscool.gradient.gradient(sinc_range, prad_range, Z, 'detections', sinc=True)
        if debug:
            ckscool.plot.planet._sinc_prad(plnt, zoom=True)
            R_sinc = log10([ckscool.gradient.R_sinc(10**i,sol_i.x[0], sol_i.x[1]) for i in sinc_range])
            plt.plot(sinc_range, R_sinc)
            plt.savefig('./test_{}'.format(seed))
            plt.clf()
            plt.close()

    else:
        X,Y,Z = ckscool.plot.planet._per_prad(plnt, for_gradient=True)
        per_range = np.array([10**i for i in X])
        prad_range = np.array([10**i for i in Y])
        sol_i = ckscool.gradient.gradient(per_range, prad_range, Z)
        if debug:
            ckscool.plot.planet._per_prad(plnt, zoom=True)
            R_per = log10([ckscool.gradient.R(10**i,sol_i.x[0], sol_i.x[1]) for i in per_range])
            plt.plot(per_range, R_per)
            plt.savefig('./test_{}'.format(seed))
            plt.clf()
            plt.close()

    return list(sol_i.x)

# ---------------------------------------------------------------------------- #

def bootstrap_chain(smass_bins, map, n_iter=10000, n_cores=1, debug=False):

    for k in range(len(smass_bins)-1):

        smass1, smass2 = smass_bins[k], smass_bins[k+1]

        if map=="detections":

            plnt = ckscool.io.load_table('planets-cuts2+iso')
            plnt = plnt[~plnt.isany]
            plnt = plnt[plnt.giso_smass.between(smass1,smass2)]
            chain = Parallel(n_jobs=n_cores)(delayed(bootstrap_detection)(plnt, seed=i, debug=debug) for i in np.arange(n_iter))
            np.savetxt("./data/chain_detv2_{0}-{1}-smass.csv".format(smass1, smass2), chain, delimiter=',')

        elif map=="occurrence":

            # Derive completeness object
            method = 'fulton-gamma-clip' # treatment for planet detectability
            impact = 0.8 # maximum impact parameter considered.

            field = ckscool.io.load_table('field-cuts',cache=1)
            field = field[~field.isany]
            field = field.rename(columns={'ber19_srad':'srad','ber19_smass':'smass'})
            field = field[field.smass.between(smass1,smass2)]
            n1 = len(field)
            field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
            n2 = len(field)
            print "{}/{} stars remain after droping nulls ".format(n2,n1)
            nstars = n2


            # Define grid of period and radius to compute completeness
            comp_per_bins = np.round(logspace(log10(0.1),log10(1000),65),4)
            comp_prad_bins = np.round(logspace(log10(0.25),log10(64),51 ),2)

            comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
            spacing_dict = {'per':'log','prad':'log'}
            grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)

            comp = ckscool.comp.Completeness(field, grid, method, impact)
            comp.compute_grid_prob_det(verbose=False)
            comp.compute_grid_prob_tr(verbose=False)
            comp.create_splines()

            # load planet population
            plnt = ckscool.io.load_table('planets-cuts2+iso')
            plnt = plnt[~plnt.isany]
            namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass'}
            plnt = plnt.rename(columns=namemap)
            plnt = plnt[plnt.smass.between(smass1,smass2)]

            # ranges must match those defined in ckscool.plot.occur.load_contour_plotter
            per_range = np.logspace(np.log10(0.1),np.log10(300),80)
            prad_range = np.logspace(np.log10(0.5),np.log10(4),80)


            chain = Parallel(n_jobs=n_cores)(delayed(bootstrap_occurrence)(smass1, smass2, plnt, comp, nstars, per_range, prad_range, seed=i, debug=debug) for i in np.arange(n_iter))
            np.savetxt("./data/chain_occv2_{0}-{1}-smass.csv".format(smass1, smass2), chain, delimiter=',')

        else:
            raise NameError(' "map" argument  must be one of the following: "occurrence" or "detections" ')

# ---------------------------------------------------------------------------- #

def bootstrap_chain_sinc(smass_bins, map, n_iter=10000, n_cores=1, debug=False):

    for k in range(len(smass_bins)-1):

        smass1, smass2 = smass_bins[k], smass_bins[k+1]

        if map=="detections":

            plnt = ckscool.io.load_table('planets-cuts2+iso')
            plnt = plnt[~plnt.isany]
            plnt = plnt[plnt.giso_smass.between(smass1,smass2)]
            chain = Parallel(n_jobs=n_cores)(delayed(bootstrap_detection)(plnt, seed=i, sinc=True, debug=debug) for i in np.arange(n_iter))
            np.savetxt("./data/chain_det_SincPradv3_{0}-{1}-smass.csv".format(smass1, smass2), chain, delimiter=',')


        elif map=="occurrence":
            # Derive completeness object
            method = 'fulton-gamma-clip' # treatment for planet detectability
            impact = 0.8 # maximum impact parameter considered.

            field = ckscool.io.load_table('field-cuts',cache=1)
            field = field[~field.isany]
            field = field.rename(columns={'ber19_srad':'srad','ber19_smass':'smass'})
            field = field[field.smass.between(smass1,smass2)]
            n1 = len(field)
            field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
            n2 = len(field)
            print "{}/{} stars remain after droping nulls ".format(n2,n1)
            nstars = n2

            # Define grid of period and radius to compute completeness
            comp_sinc_bins = np.round(logspace(log10(0.1),log10(100000),65),4)
            comp_prad_bins = np.round(logspace(log10(0.25),log10(64),51 ),2)

            comp_bins_dict = {'sinc':comp_sinc_bins, 'prad': comp_prad_bins}
            spacing_dict = {'sinc':'log','prad':'log'}
            grid = ckscool.grid.Grid(comp_bins_dict, spacing_dict)
            comp = ckscool.comp.Completeness_SincPrad(field, grid, method, impact)
            comp.compute_grid_prob_det_sinc(verbose=False)
            comp.compute_grid_prob_tr_sinc(verbose=False)
            comp.create_splines_sinc()

            # load planet population
            plnt = ckscool.io.load_table('planets-cuts2+iso')
            plnt = plnt[~plnt.isany]
            namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass','giso_sinc':'sinc'}
            plnt = plnt.rename(columns=namemap)
            plnt = plnt[plnt.smass.between(smass1,smass2)]

            # ranges must match those defined in ckscool.plot.occur.load_contour_plotter
            sinc_range = logspace(log10(0.1),log10(30000),80)
            prad_range = logspace(log10(0.5),log10(4),80)


            chain = Parallel(n_jobs=n_cores)(delayed(bootstrap_occurrence)(smass1, smass2, plnt, comp, nstars, sinc_range, prad_range, seed=i, sinc=True, debug=debug) for i in np.arange(n_iter))
            np.savetxt("./data/chain_occ_SincPradv3_{0}-{1}-smass.csv".format(smass1, smass2), chain, delimiter=',')

        else:
            raise NameError(' "map" argument  must be one of the following: "occurrence" or "detections" ')
