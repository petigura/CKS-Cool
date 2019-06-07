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


__STARS_REQUIRED_COLUMNS__ = (
    "logcdpp3 logcdpp6 logcdpp12 tobs smass srad".split()
)

class Completeness(object):
    """
    Class for that can compute completeness using the noise properties
    of a representive ensemble of stars.  
    """

    def __init__(self, stars, grid, method, impact):
        """
        Args:
            stars (pandas.DataFrame): Sample of stars from which planets 
                are detected. Must be as close as possible to be sample of 
                of stars used in the planet search. Must contain the following
                keys.
                - logcdpp3: three hour CDPP
                - logcdpp6: six
                - logcdpp12: twelve
                - tobs: days in which target was observed
                - smass: Stellar mass (solar masses) 
                - srad: Stellar radius (solar-radii) 
            grid (Grid): Grid object that contains boundaries of bins.
            method (str): method for converting per, prad into prob_det either:
                - mes-step
                - fulton-gamma: Fulton et al. fit to gamma function.
            impact: maximum impact parameter considered for our sample
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
        tdur = ( 
            TDUR_EARTH_SUN_HRS * self.stars['srad'] * 
            (per_yrs / self.stars['smass'] )**(1.0/3.0)
        )
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

    def prob_det_mes(self, mes):
        """Recovery rate vs. multiple event statistic

        Args:
            mes (array): Multiple event statistic
            
        """
        out = np.zeros_like(mes) 
        mes = np.array(mes) 
        if self.prob_det_mes_name.count('step')==1:
            min_mes = self.prob_det_mes_name.split('-')[1]
            min_mes = float(min_mes)
            out[mes > min_mes] = 1.0

        return out

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
                _prob_det = 1.0*self.prob_det_mes(mes).sum()
                _prob_det /= self.nstars

            elif self.method.count('fulton-gamma'):
                snr = self.snr(per, prad)
                _prob_det = fulton_gamma(snr).sum() 
                _prob_det/=self.nstars
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

    def compute_mes_factor(self, plnt_mes):
        """
        Using a simple SNR calculation, we tend to over estimate the MES
        that the pipeline would have found.

        Args:
            plnt_mes (pd.DataFrame): Planets used to calibrate MES. Must have 
                the fllowing keys:
                - per
                - prad
                - id_kic
                - id_koid
             
        """
        plnt_mes = plnt_mes.copy()
        self.plnt_mes = plnt_mes
        self.plnt_mes['mes_pipeline'] = self.plnt_mes.mes
        self._compute_mes_factor_loop(False)

        logmes_pipeline = np.log10(self.plnt_mes.mes_pipeline)
        logmes_formula = np.log10(self.plnt_mes.mes_formula)
        logmes_diff_med = np.nanmedian(logmes_formula - logmes_pipeline)
        logmes_diff_std = np.nanstd(logmes_formula - logmes_pipeline)
        print "med(log(mes_formula/mes_pipeline)) {:.2f} (dex)".format(
            logmes_diff_med
        )
        self.logmes_diff_std = logmes_diff_med
        self.logmes_diff_med = logmes_diff_std 
        self.mesfac = 10**(-1.0 * logmes_diff_med)
        self._compute_mes_factor_loop(True)
 
    def _compute_mes_factor_loop(self, scaled):
        i = 0
        for id_koicand, row in self.plnt_mes.iterrows():
            try:
                if scaled:
                    fmes = self.mes_scaled
                    key = 'mes_formula_scaled' 
                else:
                    fmes = self.mes
                    key = 'mes_formula' 

                mes_formula = fmes(row.per,row.prad).ix[row.id_kic]
                mes_pipeline = row.mes
                self.plnt_mes.ix[id_koicand,key] = mes_formula

                if i< 20:
                    s =  "{} {:.2f} {:.2f}".format(
                        id_koicand, mes_pipeline, mes_formula
                    )
                    print s
            except KeyError:
                pass
        
            i+=1

    def compute_grid_prob_det(self,verbose=0):
        """Compute a grid of detection probabilities"""

        print "Compute a grid of detection probabilities"
        def rowfunc(row):
            return self.prob_det(row.perc, row.pradc)

        def callback(row):
            s ="{perc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_det'] = self._grid_loop(rowfunc, callback)

    def compute_grid_prob_tr(self,verbose=0):
        """Compute a grid of transit probabilities"""

        print "Compute a grid of transit probabilities"
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

def fulton_gamma(snr):
    k = 17.56 
    l = 1.00 
    theta = 0.49
    return gamma.cdf(snr, k, l, theta)



