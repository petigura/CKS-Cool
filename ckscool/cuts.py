"""
Module for dealing with cuts.
"""
import numpy as np
import sys
import inspect

texdict ={
    'teff':'\\teff',
    'logg':'\\logg'
}

plotdict ={
    'teff':'\mathregular{T}_\mathregular{eff}',
    'logg':'\log g'
}

samples = 'koi-thompson18'.split()

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
    plotstr = '$Kp$ < 16.0 mag'
    texstr = '$Kp$ < 16.0 mag'
    def cut(self):
        kepmag = self.df['m17_kepmag']
        b = kepmag > 16.0
        return b

class CutTeffPhot(CutBase):
    """Remove stars where Teff falls outside allowed range
    """
    cuttype = 'badteffphot'
    texstr = r'${teff:}$ < $4900$ K'.format(**texdict)
    plotstr = r'${teff:}$ < 4900 K'.format(**plotdict)
    def cut(self):
        teff = self.df['m17_steff']
        b = ~teff.between(0,4900)
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
            b2 = self.df.koi_score < 0.75 
            return b1 | b2
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
