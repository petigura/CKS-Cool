"""
Module for dealing with cuts.
"""
import numpy as np
import sys
import inspect
from collections import OrderedDict

texdict ={
    'teff':'\\teff',
    'logg':'\\logg'
}

plotdict ={
    'teff':'\mathregular{T}_\mathregular{eff}',
    'logg':'\log g'
}

#samples = 'koi-mullally15 koi-thompson18 ckscool-planets'.split()

cd = OrderedDict()
cd['max-kepmag'] = 16.0
cd['max-steff'] = 5000
cd['min-steff'] = 3500
cd['max-srad'] = 1.5
cd['max-score'] = 0.75

class CutBase(object):
    cuttype = 'base'
    plotstr = 'base'
    texstr = 'base'

    def __init__(self, df, sample):
        #assert samples.count(sample)==1,'{} not supported sample'.format(sample)
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

class CutFaint(CutBase):
    cuttype = 'faint'
    plotstr = '$Kp$ < {} mag'.format(cd['max-kepmag'])
    texstr = '$Kp$ < {} mag'.format(cd['max-kepmag'])
    def cut(self):
        kepmag = self.df['kic_kepmag']
        b = kepmag > cd['max-kepmag']
        return b

class CutGiant(CutBase):
    cuttype = 'giant'
    plotstr = 'Not giant'.format(cd['max-kepmag'])
    texstr = 'Not giant'.format(cd['max-kepmag'])
    def cut(self):
        b = self.df['gaia2_srad'] > cd['max-srad']
        return b

class CutBadParallax(CutBase):
    cuttype = 'badparallax'
    plotstr = '$\sigma(\pi) / \pi < 10$'
    texstr = '$\sigma(\pi) / \pi < 10$'
    def cut(self):
        b1 = self.df['gaia2_sparallax'].isnull()
        b2 = self.df['gaia2_sparallax_over_err'] < 10 
        b = b1 | b2
        return b

class CutDiluted(CutBase):
    cuttype = 'diluted'
    plotstr = 'Not diluted'
    texstr = 'Not diluted'
    def cut(self):
        b = self.df['gaia2_gflux_ratio'] > 1.1
        return b

class CutFurlan(CutBase):
    cuttype = 'dilutedfurlan'
    plotstr = 'furlan RCF < 1.05'
    texstr = 'furlan RCF < 1.05'
    def cut(self):
        b = self.df['f17_rcf_avg'] > 1.05
        return b

class CutAO(CutBase):
    cuttype = 'dilutedao'
    plotstr = 'No companions with $\Delta K < 2.5$'
    texstr = 'No companions with $\Delta K < 2.5$'
    def cut(self):
        b1 = self.df['f17_rcf_avg'] > 1.05

        # K4 star has M_K = 4.31 and M = 0.73 Msun. Want to find the
        # mass of a star with 10% of flux in So looking for star with
        # M_K 4.31 + 2.5 4.31 + 2.5 = 6.81 7.55 M4 star has M_K 7.55
        # and M = 0.24 Msun. So using a mass ratio of 0.3 will remove
        # stars with a dilutoin factor of 10% in Kmag.
        b2 = self.df['k16_max_massratio'] > 0.3
        return b1 | b2 

class CutTeffPhot(CutBase):
    """Remove stars where Teff falls outside allowed range
    """
    cuttype = 'badteffphot'
    texstr = r'${teff:}$ = ${min-steff:}$-${max-steff:}$ K'.format(**dict(texdict, **cd))
    plotstr = r'${teff:}$ = ${min-steff:}$-${max-steff:}$ K'.format(**dict(plotdict, **cd))
    def cut(self):
        teff = self.df['m17_steff']
        b = ~teff.between(cd['min-steff'],cd['max-steff'])
        return b

class CutNotReliable(CutBase):
    """Remove stars Remove stars with a FP designation from a number of catalogs
    """
    cuttype = 'notreliable'
    plotstr = 'Reliable KOI'
    texstr = plotstr
    def cut(self):
        if self.sample=='koi-thompson18':
            b1 = self.df.koi_disposition.str.contains('FALSE')  
            b2 = self.df.koi_score < cd['max-score']
            return b1 | b2
        elif self.sample=='koi-mullally15':
            b1 = self.df.koi_disposition.str.contains('FALSE')  
            return b1
        else:
            assert False," error"

####################################
# Cuts for filtered stellar sample #
####################################

class CutSpecParallax(CutBase):
    """Remove stars the spectroscopic parallax does not agree with the trig parallax
    """
    cuttype = 'sb2'
    texstr = r'Removing SB2s'
    plotstr = r'Removing SB2s'
    def cut(self):

        diff = self.df.gaia2_sparallax - self.df.giso2_sparallax
        sigmadiff = np.sqrt(
            self.df.giso2_sparallax_err1**2 + self.df.gaia2_sparallax_err**2 
        ) 
        ndiff = diff/sigmadiff
        teff = self.df['m17_steff']
        b = np.abs(ndiff) > 3
        return b

class CutSB2(CutBase):
    """Remove stars that are SB2s"""
    cuttype = 'badspecparallax'
    texstr = r'Delta Parallax < 3 sigma'
    plotstr = r'$|\pi_{trig} - \pi_{spec}| < 3 \sigma$'
    def cut(self):
        b = self.df.rm_sb2.isin([4,5])
        return b

###################################
# Cuts for filtered planet sample #
###################################

class CutImpact(CutBase):
    """Remove planets with high impact parameters"""
    cuttype = 'badimpact'
    texstr = r'$b$ < 0.8'
    plotstr = texstr
    def cut(self):
        b = self.df.koi_impact > 0.8
        return b

class CutPradPrecision(CutBase):
    """Remove planets with large radius uncertainties"""
    cuttype = 'badpradprecision'
    texstr = r'$\sigma(R_p) / R_p $ < 0.2'
    plotstr = texstr
    def cut(self):
        b = self.df.eval('gdir_prad_err1/gdir_prad') > 0.2
        return b

class CutPrad(CutBase):
    """Remove planets with large radius uncertainties"""
    cuttype = 'badprad'
    texstr = r'$R_p < 20 R_E$'
    plotstr = texstr
    def cut(self):
        b = ~self.df.gdir_prad.between(0,20)
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
    return df 
