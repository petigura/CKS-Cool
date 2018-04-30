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

samples = 'koi-mullally15 koi-thompson18'.split()

cd = OrderedDict()
cd['max-kepmag'] = 15.8
cd['max-steff'] = 5200
cd['min-steff'] = 3500
cd['max-score'] = 0.75

class CutBase(object):
    cuttype = 'base'
    plotstr = 'base'
    texstr = 'base'

    def __init__(self, df, sample):
        assert samples.count(sample)==1,'{} not supported sample'.format(sample)
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
        b = self.df['gaia2_srad'] > 1.5
        return b

class CutBadParallax(CutBase):
    cuttype = 'badparallax'
    plotstr = 'Bad Parallax'
    texstr = 'Bad parallax'
    def cut(self):
        b = self.df['gaia2_sparallax_over_err'] < 25
        return b

class CutDiluted(CutBase):
    cuttype = 'diluted'
    plotstr = 'Diluted'
    texstr = 'Diluted'.format(cd['max-kepmag'])
    def cut(self):
        b = self.df['gaia2_gflux_ratio'] > 1.1
        return b

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
