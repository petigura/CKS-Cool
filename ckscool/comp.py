"""
Completeness class
"""
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import scipy.integrate
from scipy.stats import gamma
from astropy import constants as c
from astropy import units as u

TDUR_EARTH_SUN_HRS = (
    ((4 * c.R_sun**3 * 1.0*u.yr / np.pi / c.G / (1.0*c.M_sun))**(1.0/3.0)).to(u.hr)).value

DEPTH_EARTH_SUN = ((c.R_earth / c.R_sun)**2).cgs.value

SINC_EARTH = 1.3608e3

__STARS_REQUIRED_COLUMNS__ = (
    "logcdpp3 logcdpp6 logcdpp12 tobs smass srad".split()
)

class CompletenessPerPrad(object):
    """Class that compute completeness using the noise properties of a
    representive ensemble of stars.

    Examples:

        # 1. Define sample of stars.
        >>> stars = ckscool.io.load_table('field-cuts',cache=1)
        >>> stars = stars[~stars.isany]
        >>> stars = stars.rename(
               columns={'ber18_srad':'srad','m17_smass':'smass'}
            )
        >>> stars = stars.query('smass > 1.3')

        # 2. Define grid to compute completeness
        comp_per_bins = round(logspace(log10(0.1),log10(1000),33),4)
        comp_prad_bins = round(logspace(log10(0.25),log10(64),25 ),2)
        comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
        spacing_dict = {'per':'log','prad':'log'}
        grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)

        # 3. Compute completeness grid
        comp = ckscool.comp.Completeness(stars, grid, method, impact)
        comp.compute_grid_prob_det(verbose=True)
        comp.compute_grid_prob_tr(verbose=True)

        # 4. Compute spline interpolation over completeness grid
        comp.create_splines()
        comp.mean_prob_trdet(1,10,1,2)

    """

    def __init__(self, stars, grid, method, impact, mesfac=None):
        """Args:
            stars (pandas.DataFrame): Sample of stars from which planets
                are detected. Must be as close as possible to be sample of
                of stars used in the planet search. Must contain the following
                keys.
                - logcdpp3: three-hour CDPP noise
                - logcdpp6: six-hour CDPP noise
                - logcdpp12: twelve-hour CDPP noise
                - tobs: days in which target was observed
                - smass: Stellar mass (solar masses)
                - srad: Stellar radius (solar-radii)
            grid (Grid): Grid object that contains boundaries of bins.
            method (str): method for converting per, prad into prob_det either:
                - mes-step
                - fulton-gamma: Fulton et al. fit to gamma function.
            impact: maximum impact parameter considered for our sample
            mesfactor: conversion between theoretical SNR and actual MES

        Notes:

            The Kepler pipeline detects planets by computing the
            multiple event statistic (MES) over a range of periods and
            epochs. The MES is like a signal to noise ratio, but
            beacuse the pipeline computes MES over a finite grid, the
            MES usually lower than the SNR computed after fitting a
            model to the light curve. When we compute what the SNR
            would be of a putative planet at (per,prad), we are
            implicitly assuming that we would have computed it at the
            correct period, epoch. So to convert the SNR to an
            expected MES we have to apply a correction of ~20%. Note
            this is only relevant to the mes-step method of converting
            computing the completeness. The Fulton-gamma perscription
            fit the recovered planets as a function of the theoretical
            SNR so it already incorporates this bias.

        """

        for col in __STARS_REQUIRED_COLUMNS__:
            assert list(stars.columns).count(col)==1,\
                "required column {} missing".format(col)

            assert stars[col].notnull().all(),\
                "column {} cannot contain nulls".format(col)

        self.stars = stars

        # Define the regular grid for interpolation
        x0 = np.array(stars.index) # (nstar)
        x1 = np.log10([3,6,12]) # (3)
        points  = (x0, x1)
        cols = 'logcdpp3 logcdpp6 logcdpp12'.split()
        values = stars[cols] # (nstar, 3)
        values = np.array(values)
        logcdpp_interp = RegularGridInterpolator(
            points, values, bounds_error=False, fill_value=None
        )

        self.logcdpp_interp = logcdpp_interp
        self.x0 = x0
        self.x1 = x1
        self.values = values
        self.nstars = len(stars)
        self.grid = grid
        self.method = method
        self.impact = impact
        self.mesfac = mesfac

    def snr(self, per, prad):
        """
        Calculate expected transit SNR

        For a planet of a given period and size, calculate the
        expected SNR for every star in the sample.

        Args:
            per (float) : orbital period of planet days
            prad (float): size of planet in Earth-radii

        Returns:
            array: SNR

        """
        depth = self._depth(prad)
        tdur = self._tdur(per)
        cdpp = self._cdpp(tdur) * 1e-6
        num_transits = self._num_transits(per)
        snr = depth / cdpp * np.sqrt(num_transits)
        return snr

    def mes_scaled(self, per, prad):
        """
        Calculate Multiple Event Statistic and apply scaling
        """
        return self.mesfac * self.mes(per, prad)

    def _tdur(self, per):
        """
        Compute duration for a putative transit of a given orbital period

        Args:
            per (float): orbital period

        Returns:
            pandas.Series: transit duration for each star in the sample.
        """
        per_yrs = per / 365.25
        srad = self.stars['srad']
        smass = self.stars['smass']
        tdur =  TDUR_EARTH_SUN_HRS * srad * (per_yrs / smass )**0.33
        return tdur

    def _smax(self, per):
        """
        Compute semi-major axis putative transit of a given orbital period

        Args:
            per (float): orbital period

        Returns:
            pandas.Series: transit duration for each star in the sample.
        """
        per_yrs = per / 365.25
        smax = self.stars['smass']**(1.0/3.0) * per_yrs**(2.0/3.0)
        return smax

    def _depth(self, prad):
        """
        Compute transit depth for a putative planet of a given size.

        Args:
            prad (float): planet radius (Earth-radii)

        Returns:
            pd.Series: the depth of a `prad` planet around each star

        """

        _depth = (prad / self.stars.srad)**2 * DEPTH_EARTH_SUN
        return _depth

    def _num_transits(self, per):
        """
        Compute number of transits for a planet having a given orbital
        period, factoring in duty cycle.

        Args:
            per (float): orbital period (days)

        Returns:
            pd.Series: the number of transits for each star in the sample
        """

        return self.stars.tobs / per

    def _cdpp(self, tdur):
        """
        Calculate noise (CDPP) over a specified duration

        Args:
            tdur (pandas.Series or array): the transit duration for each star

        Retruns:
            pd.Series: The CDPP for each star

        """
        tdur = np.array(tdur)
        logtdur = np.log10(tdur)

        x0i = self.x0
        x1i = logtdur
        pointsi = np.vstack([x0i,x1i]).T
        logcdpp = self.logcdpp_interp(pointsi)
        cdpp = pd.Series(index=self.stars.index,data=10**logcdpp)
        return cdpp

    def prob_det(self, per, prad, interp=False):
        """Probability that a planet would be detectable

        Probability that transiting planet with orbital period `per`
        and planet radius `prad` would be detectable around a
        randomly-drawn star from the stellar sample.

        Args:
            per (float): orbital period
            prad (float): planet size
            method (str): One of the following:
                - direct: compute prob of detectability by counting stars
                  where the MES is sufficient to detect.
                - interp: use the interpolator

        Returns:
            float: fraction of stars in sample where we could have detected
                planet

        """
        if interp==False:
            if self.method.count('step'):
                mes = self.mes_scaled(per, prad)
                _prob_det = 1.0*self.prob_det_mes(mes).sum() / self.nstars

            elif self.method=='fulton-gamma':
                snr = self.snr(per, prad)
                _prob_det = fulton_gamma(snr).sum() / self.nstars

            elif self.method=='fulton-gamma-clip':
                snr = self.snr(per, prad)
                _prob_det = fulton_gamma_clip(snr).sum() / self.nstars
        else:
            _prob_det = self.prob_det_interp(per, prad)

        return _prob_det

    def prob_tr(self, per):
        """
        Probability that a given planet would transit
        """
        a = self._smax(per)
        srad = (np.array(self.stars['srad']) * u.R_sun).to(u.AU).value
        _prob_tr = srad * self.impact / a
        return _prob_tr.mean()

    def compute_grid_prob_det(self, verbose=0):
        """Compute a grid of detection probabilities"""

        print "Computing grid of detection probabilities"
        def rowfunc(row):
            return self.prob_det(row.perc, row.pradc)

        def callback(row):
            s ="{perc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_det'] = self._grid_loop(rowfunc, callback)

    def compute_grid_prob_tr(self,verbose=0):
        """Compute a grid of transit probabilities"""

        print "Computing grid of transit probabilities"
        def rowfunc(row):
            return self.prob_tr(row.perc)

        def callback(row):
            s ="{perc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_tr'] = self._grid_loop(rowfunc, callback)

    def _grid_loop(self, rowfunc, callback):
        df = self.grid.ds.to_dataframe()
        i = 0
        for idx, row in df.iterrows():
            df.ix[idx,'temp'] = rowfunc(row)
            if i%100==0:
                callback(df.ix[idx])
            i+=1
        return df['temp'].to_xarray()

    def create_splines(self):
        """Intialize detection probability interpolator
        """
        # Define the regular grid for interpolation
        x0 = np.log(np.array(self.grid.ds.per))
        x1 = np.log(np.array(self.grid.ds.prad))
        points = (x0, x1)

        prob_det = self.grid.ds['prob_det'].transpose('per','prad')
        prob_trdet = self.grid.ds['prob_tr'] * self.grid.ds['prob_det']
        prob_trdet = prob_trdet.transpose('per','prad')
        self._prob_trdet_spline = RectBivariateSpline(x0, x1, prob_trdet)
        self._prob_det_spline = RectBivariateSpline(x0, x1, prob_det)

    def prob_det_interp(self, per, prad):
        _prob_det =  self._prob_det_interp_spline(per, prad)
        assert np.isfinite(_prob_det), " error in computing transit prob"
        return _prob_det

    def mean_prob_trdet(self, per1, per2, prad1, prad2):
        """Mean probability of transiting and being detected"""
        x1 = np.log(per1)
        x2 = np.log(per2)
        y1 = np.log(prad1)
        y2 = np.log(prad2)

        area = (x2-x1) * (y2-y1)
        integral = self._prob_trdet_spline.integral(x1, x2, y1, y2)
        prob_trdet_mean = integral / area

        integral = self._prob_det_spline.integral(x1, x2, y1 ,y2)
        prob_det_mean = integral / area
        return prob_trdet_mean, prob_det_mean

    def prob_trdet_interp(self, per, prad):
        x = np.log(per)
        y = np.log(prad)
        return self._prob_trdet_spline(x,y,grid=False)
    
class CompletenessSincPrad(CompletenessPerPrad):
    def _a(self, sinc, srad, steff, smass):
        """Returns semimajor axis for a given incident flux,
        stellar mass, stellar radius and stellar effective temperature
        """
        _a = np.sqrt( ( (c.R_sun*srad)**2 * c.sigma_sb * (steff**4) ) / ( SINC_EARTH*sinc ) )
        return _a.value

    def _per(self, a, smass):
        """ Returns the orbital period for a given semimajor axis and
        stellar mass
        """
        _per = 2.0 * np.pi * np.sqrt(a**3 / (c.G * smass * c.M_sun) )
        _per_days = _per / (60.0 * 60.0 * 24.0)
        return _per_days.value

    def prob_det(self, sinc, prad, interp=False):
        """Probability that a planet would be detectable

        Probability that transiting planet with incident flux `sinc`
        and planet radius `prad` would be detectable around a
        randomly-drawn star from the stellar sample.

        Args:
            sinc (float): stellar incident flux
            prad (float): planet size
            method (str): One of the following:
                - direct: compute prob of detectability by counting stars
                  where the MES is sufficient to detect.
                - interp: use the interpolator

        Returns:
            float: fraction of stars in sample where we could have detected
                planet

        """
        srad = self.stars['srad']
        steff = self.stars['gaia2_steff']
        smass = self.stars['smass']

        a = self._a(sinc, srad, steff, smass)
        per = self._per(a, smass)

        if interp==False:
            if self.method.count('step'):
                mes = self.mes_scaled(per, prad)
                _prob_det = 1.0*self.prob_det_mes(mes).sum() / self.nstars

            elif self.method=='fulton-gamma':
                snr = self.snr(per, prad)
                _prob_det = fulton_gamma(snr).sum() / self.nstars

            elif self.method=='fulton-gamma-clip':
                snr = self.snr(per, prad)
                _prob_det = fulton_gamma_clip(snr).sum() / self.nstars
        else:
            _prob_det = self.prob_det_interp(per, prad)

        return _prob_det

    def compute_grid_prob_det(self,verbose=0):
        """Compute a grid of detection probabilities"""

        print "Computing grid of detection probabilities"
        def rowfunc(row):
            return self.prob_det_sinc(row.sincc, row.pradc)

        def callback(row):
            s ="{sincc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_det'] = self._grid_loop(rowfunc, callback)

    def compute_grid_prob_tr(self,verbose=0):
        """Compute a grid of transit probabilities"""

        print "Computing grid of transit probabilities"
        def rowfunc(row):
            srad = self.stars['srad']
            steff = self.stars['gaia2_steff']
            smass = self.stars['smass']

            a_c = self._a(row.sincc, srad, steff, smass)
            perc = self._per(a_c, smass)
            return self.prob_tr(perc)

        def callback(row):
            s ="{sincc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_tr'] = self._grid_loop(rowfunc, callback)

    def create_splines(self):
        """Intialize detection probability interpolator
        """
        # Define the regular grid for interpolation
        x0 = np.log(np.array(self.grid.ds.sinc))
        x1 = np.log(np.array(self.grid.ds.prad))
        points = (x0, x1)

        prob_det = self.grid.ds['prob_det'].transpose('sinc','prad')
        prob_trdet = self.grid.ds['prob_tr'] * self.grid.ds['prob_det']
        prob_trdet = prob_trdet.transpose('sinc','prad')
        self._prob_trdet_spline = RectBivariateSpline(x0, x1, prob_trdet)
        self._prob_det_spline = RectBivariateSpline(x0, x1, prob_det)

    def prob_det_interp(self, sinc, prad):
        _prob_det =  self._prob_det_interp_spline(per, prad)
        assert np.isfinite(_prob_det), " error in computing transit prob"
        return _prob_det

    def mean_prob_trdet(self, sinc1, sinc2, prad1, prad2):
        """Mean probability of transiting and being detected"""
        x1 = np.log(sinc1)
        x2 = np.log(sinc2)
        y1 = np.log(prad1)
        y2 = np.log(prad2)
        area = (x2-x1) * (y2-y1)
        integral = self._prob_trdet_spline.integral(x1, x2, y1, y2)
        prob_trdet_mean = integral / area

        integral = self._prob_det_spline.integral(x1, x2, y1 ,y2)
        prob_det_mean = integral / area
        return prob_trdet_mean, prob_det_mean

import copy
class Completeness3D_PerPradSmass(object):
    """
    
    """
    def __init__(self, stars, grid, massc, masswid, method, impact, mesfac=None):
        # Derive completeness object
        method = 'fulton-gamma-clip' # treatment for planet detectability
        impact = 0.8 # maximum impact parameter considered.

        self.stars = stars

        # default values. can be changed
        self.massc = massc
        self.masswid = masswid
        self.grid = grid
        self.method = method
        self.impact = impact

    def get_comp2d(self, smass1, smass2):
        stars = self.stars[self.stars.smass.between(smass1,smass2)].copy()
        n1 = len(stars)
        stars = stars.dropna(subset=__STARS_REQUIRED_COLUMNS__)
        n2 = len(stars)
        print "{}/{} stars remain after droping nulls ".format(n2,n1)
        comp = Completeness(stars, self.grid, self.method, self.impact)
        comp.compute_grid_prob_det(verbose=False)
        comp.compute_grid_prob_tr(verbose=False)
        comp = copy.deepcopy(comp) # don't understand why I need this but I d
        return comp

    def compute_completeness_grid(self):
        compL = []
        for i in range(self.massc.shape[0]):
            smass1 = self.massc[i] / 10**(0.5*self.masswid)
            smass2 = self.massc[i] * 10**(0.5*self.masswid)
            comp = self.get_comp2d(smass1,smass2)
            compL.append(comp)
        self.compL = compL

    def create_interpolators(self):
        comp0 = self.compL[0]
        xi = np.log10(self.grid.binsc['per'])
        yi = np.log10(self.grid.binsc['prad'])
        zi = np.log10(self.massc)

        prob_tr = []
        prob_det = []
        for comp in self.compL:
            prob_tr.append(np.array(comp.grid.ds.prob_tr))
            prob_det.append(np.array(comp.grid.ds.prob_det))
        
        prob_tr = np.dstack(prob_tr)
        prob_det = np.dstack(prob_det)
        self._prob_tr = RegularGridInterpolator((xi,yi,zi), prob_tr)
        self._prob_det = RegularGridInterpolator((xi,yi,zi), prob_det)

    def _transform_to_grid(self, per, prad, smass):
        points = np.vstack([np.array(per),np.array(prad),np.array(smass)]).T
        points = np.log10(points)
        return points

    def prob_det(self, per, prad, smass):
        return self._prob_det(self._transform_to_grid(per, prad, smass))

    def prob_tr(self, per, prad, smass):
        return self._prob_tr(self._transform_to_grid(per, prad, smass))
        
def load_occur(limits, debug=False, sinc=False):
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
    namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass','giso_sinc':'sinc'}
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
    field = field.dropna(subset=__STARS_REQUIRED_COLUMNS__)
    n2 = len(field)
    print "{}/{} stars remain after droping nulls ".format(n2,n1)

    if sinc:
        comp_sinc_bins = np.round(logspace(log10(0.1),log10(100000),65),4)
        comp_prad_bins = np.round(logspace(log10(0.25),log10(64),51 ),2)

        # debugging
        if debug:
            comp_sinc_bins = comp_sinc_bins[:6]
            comp_prad_bins = comp_prad_bins[:6]

        comp_bins_dict = {'sinc':comp_sinc_bins, 'prad': comp_prad_bins}
        spacing_dict = {'sinc':'log','prad':'log'}
        grid = ckscool.grid.Grid(comp_bins_dict, spacing_dict)
        comp = ckscool.comp.Completeness_SincPrad(field, grid, method, impact)
        comp.compute_grid_prob_det_sinc(verbose=False)
        comp.compute_grid_prob_tr_sinc(verbose=False)
        comp.create_splines_sinc()
        nstars = len(field)
        occ = ckscool.occur.Occurrence_SincPrad(plnt, comp, nstars)

    else:
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



    
# ---------------------------------------------------------------------------- #

def fulton_gamma(snr):
    k = 17.56
    l = 1.00
    theta = 0.49
    return gamma.cdf(snr, k, l, theta)

def fulton_gamma_clip(snr):
    return fulton_gamma(snr) * (snr > 10)

