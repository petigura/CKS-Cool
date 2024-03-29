"""
Defining cuts for the occurrence analysis
"""

import numpy as np
import sys
import inspect
from collections import OrderedDict
from numpy import log10

cd = OrderedDict()
cd['max-kepmag'] = 16.0
cd['max-steff'] = 5000
cd['min-steff'] = 3500
cd['max-srad'] = 1.5
cd['max-score'] = 0.75


texdict ={
    'teff':'\\teff',
    'logg':'\\logg'
}

plotdict ={
    'teff':'\mathregular{T}_\mathregular{eff}',
    'logg':'\log g'
}

class CutBase(object):
    cuttype = 'base'
    plotstr = 'base'
    texstr = 'base'

    def __init__(self, df, sample):
        self.df = df
        self.sample = sample 
        self.nrows = len(df)

    def allpass(self):
        """A noop that lets all values pass
        """
        return np.zeros(self.nrows).astype(bool)

class CutNone(CutBase):
    cuttype = 'none'
    plotstr = 'Full sample'
    texstr = 'Full sample'
    def cut(self):
        return self.allpass()

##############################################################
# Cuts applied to both parent stellar population and planets #
##############################################################

class CutDiluted(CutBase):
    cuttype = 'diluted'
    plotstr = 'Not diluted'
    texstr = 'Not diluted'
    def cut(self):
        b = self.df['gaia2_gflux_dilution'] > 1.1
        return b

class CutRizzuto(CutBase):
    cuttype = 'rizzuto'
    plotstr = 'No excess astrometric noise'
    texstr = plotstr
    def cut(self):
        b = self.df['gaia2_astrometric_excess_noise_sig'] > 10
        return b

class CutRUWE(CutBase):
    cuttype = 'ruwe'
    plotstr = 'RUWE < 1.2'
    texstr = plotstr
    def cut(self):
        b = self.df['gaia2_ruwe'] > 1.2
        return b

class CutGiant(CutBase):
    cuttype = 'giant'
    plotstr = 'Not giant'
    texstr = 'Not giant'
    def cut(self):
        srad = self.df['gaia2_srad']
        steff = self.df['m17_steff']
        b = np.log10(srad) > 0.00015 * (steff - 5500) + 0.20
        return b

class CutGiantCMD(CutBase):
    """
    Cut out giants based on CMD
    """
    cuttype = 'giantcmd'
    plotstr = 'Not giant'
    texstr = 'Not giant'
    


    points1 = np.array(
        [[0.5,2],
         [0.8,2],
         [1.25,5.7],
         [2.5,9.25]]
    )

    points2 = np.array(
        [[.5,3.25],
         [0.8,5.5],
         [2.5,10.4]]
    )

    assert points1[0][0]==points2[0][0]
    assert points1[-1][0]==points2[-1][0]

    def fabove(self,x):
        return np.interp(x,self.points1[:,0],self.points1[:,1])

    def fbelow(self,x):
        return np.interp(x,self.points2[:,0],self.points2[:,1])

    def cut(self):
        m = self.df
        m['dmod'] = 5 * log10(m.bj_dist/10)
        xs = 'gaia2_bpmag - gaia2_rpmag'
        ys =  'gaia2_gmag - dmod'
        babove = m.eval(ys) < self.fabove(m.eval(xs))
        bbelow = m.eval(ys) > self.fbelow(m.eval(xs))
        btoored = m.eval(xs) > self.points1[-1][0]
        btooblue = m.eval(xs) < self.points1[0][0]
        return babove | bbelow | btoored | btooblue

class CutFaint(CutBase):
    cuttype = 'faint'
    plotstr = cuttype
    texstr = cuttype
    def cut(self):
        steff = self.df['m17_steff']
        kepmag = self.df['m17_kepmag']
        b1 = steff.between(5000,6500) & kepmag.between(0,14.2)
        b2 = steff.between(3000,5000) & kepmag.between(0,16)
        b = ~(b1 | b2)
        return b



##########################################
# Cuts only applied to planet population #
##########################################

class CutSB2(CutBase):
    """Remove stars that are SB2s"""
    cuttype = 'sb2'
    texstr = r'Removing SB2s'
    plotstr = r'Removing SB2s'
    def cut(self):
        b = self.df.rm_sb2.isin([4,5])
        return b

class CutPrad(CutBase):
    """Remove planets with large radius uncertainties"""
    cuttype = 'badprad'
    texstr = r'$R_p < 20 R_E$'
    plotstr = texstr
    def cut(self):
        b = ~self.df.gdir_prad.between(0,20)
        return b

class CutNotReliable(CutBase):
    """Remove stars Remove stars with a FP designation from a number of catalogs

    koi_disposition is the disposition on the exoplanet archive.
    koi_disposition_sup is the disposition on the exoplanet archive supplemental
    koi_pdisposition is the disposition based on kepler data

    must not be listed as a false positive based on kepler project,
    with a score of > 0.75. Must also not be listed as a false positive on exoplanet archive

    """
    cuttype = 'notreliable'
    plotstr = 'Reliable KOI'
    texstr = plotstr
    def cut(self):
        if self.sample=='koi-thompson18':
            b1 = self.df.koi_pdisposition.str.contains('FALSE')  
            b2 = self.df.koi_score < cd['max-score']
            b3 = self.df.koi_disposition_sup.str.contains('FALSE')
            return b1 |  b2 | b3
        elif self.sample=='koi-mullally15':
            b1 = self.df.koi_disposition.str.contains('FALSE')  
            return b1
        else:
            assert False," error"

class CutSNR(CutBase):
    cuttype = 'lowsnr'
    plotstr = 'snr'
    texstr = 'snr'
    def cut(self):
        b = self.df['koi_max_mult_ev'] < 10
        return b

class CutMCMC(CutBase):
    cuttype = 'nomcmc'
    plotstr = 'mcmc'
    texstr = 'mcmc'
    def cut(self):
        b = ~self.df.koi_fittype.str.contains('MCMC')
        return b

class CutCROWDSAP(CutBase):
    cuttype = 'dilutedsap'
    plotstr = 'sap'
    texstr = 'sap'
    def cut(self):
        b = ~self.df.koi_fittype.str.contains('MCMC')
        return b

class CutImpact(CutBase):
    """Remove planets with high impact parameters"""
    cuttype = 'badimpact'
    texstr = r'$b$ < 0.8'
    plotstr = texstr
    def cut(self):
        b = self.df.dr25_b > 0.8
        return b

class CutImpactTau(CutBase):
    """Remove planets with high impact parameters"""
    cuttype = 'badimpacttau'
    texstr = r'$\tau / \tau_0$ > 0.6'
    plotstr = texstr
    def cut(self):
        b = self.df.eval('dr25_period > 3 and dr25_tau / giso_tau0 < 0.6')
        return b

class CutPradPrecision(CutBase):
    """Remove planets with large radius uncertainties"""
    cuttype = 'badpradprec'
    texstr = r'$\sigma(R_p) / R_p $ < 0.2'
    plotstr = texstr
    def cut(self):
        b = self.df.eval('gdir_prad_err1/gdir_prad > 0.2') 
        return b

# Not really necessary
class CutSpecParallax(CutBase):
    """Remove stars the spectroscopic parallax does not agree with the trig parallax
    """
    cuttype = 'badspecparallax'
    texstr = r'Delta Parallax < 4 sigma'
    plotstr = r'$|\pi_{trig} - \pi_{spec}| < 4 \sigma$'
    def cut(self):

        diff = self.df.gaia2_sparallax - self.df.giso2_sparallax
        sigmadiff = np.sqrt(
            self.df.giso2_sparallax_err1**2 + self.df.gaia2_sparallax_err**2 
        ) 
        ndiff = diff/sigmadiff
        teff = self.df['m17_steff']
        b = np.abs(ndiff) > 4
        return b

class CutVsinI(CutBase):
    """Remove stars the spectroscopic parallax does not agree with the trig parallax
    """
    cuttype = 'badvsini'
    texstr = r'VsinI > 10'
    plotstr = r'$v \sin i > 20$ km/s'
    def cut(self):
        b = self.df['cks_svsini'] > 20 
        return b


class CutFPP(CutBase):
    """Remove planets with high false positive probability"""
    cuttype = 'largefpp'
    texstr = r'FPP < 0.9'
    plotstr = texstr
    def cut(self):
        b = self.df.fpp_prob > 0.9
        return b

clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

def get_cut(cuttype):
    for name, obj in clsmembers:
        if name.count('Cut'):
            if getattr(obj,'cuttype')==cuttype:
                return obj

    assert False, "{} not defined".format(cuttype)

def add_cuts(df, cuttypes, sample):
    isany = np.zeros(len(df))
    for cuttype in cuttypes:
        obj = get_cut(cuttype)
        cut = obj(df,sample)
        key = 'is'+cuttype
        df[key] = cut.cut()
        isany += df[key].astype(int)
    
    df['isany'] = isany > 0 
    df.cuttypes=cuttypes
    return df 
