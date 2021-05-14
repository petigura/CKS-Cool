from collections import OrderedDict

import numpy as np
import pandas as pd

TAU_CONST = 2.036 # duration [hr] of a P = 1 day planet orbiting a 1 g/cc star


def update_planet_parameters(df):
    """Update planet parameters

    Args:
         df: DataFrame

    Returns
        pandas.DataFrame: with following columns added

         -    gdir_prad and [_err1 and _err2]
         -    giso_prad
         -    giso_sinc
         -    giso_sma
         -    giso_tau0

    """
    nsamp = 10000

    # load up chains:
    namemap = {'RD1':'dr25_ror','RHO':'dr25_rho','BB1':'dr25_b','PE1':'dr25_period','DIL':'dr25_dil'}
    samp = dict(dr25_ror=[],dr25_rho=[],dr25_b=[],dr25_period=[],dr25_dil=[]) 
    print("loading up {} chains".format(len(df)))
    drop = []
    for i, row in df.iterrows():
        if i % 100 == 0:
            print(i)
        fn = 'data/dr25-chains_trimmed-thinned.hdf'
        try:
            with pd.HDFStore(fn) as store:
                chain = pd.read_hdf(store,row.id_koicand)
                chain2 = chain.dropna(subset=namemap.keys())
                diff = len(chain) - len(chain2) 
                if diff > 0:
                    print("dropped {} nan rows from {} chains"
                          .format(diff,row.id_koicand))
                chain = chain2

        except KeyError:
            drop.append(row.id_koicand)
            print("no chains for {}".format(row.id_koicand))
            continue

        for key,val in namemap.iteritems():
            samp[val].append(chain[key].sample(nsamp,replace=True))

    df = df.set_index('id_koicand').drop(drop).reset_index()
            
    for k in samp.keys():
        samp[k] = np.vstack(samp[k]).astype(float)
        
    samp['dr25_tau'] = (
        TAU_CONST
        * (1 - samp['dr25_b']**2)**0.5
        * samp['dr25_period']**0.33
        * samp['dr25_rho']**-0.33
    )
    def sample(loc, scale, seed):
        np.random.RandomState(seed)
        loc = np.array(loc).reshape(-1, 1)
        scale = np.array(scale).reshape(-1, 1)
        return loc + scale * np.random.randn(1, nsamp)


    # stellar radius direct method
    loc = df.gdir_srad
    scale = df.eval('0.5 * (gdir_srad_err1 - gdir_srad_err2)')
    samp['gdir_srad'] = sample(loc, scale, 0)

    # stellar radius grid method
    loc = df.giso_srad
    scale = df.eval('0.5 * (giso_srad_err1 - giso_srad_err2)')
    samp['giso_srad'] = sample(loc, scale, 1)

    # stellar mass
    loc = df.giso_smass
    scale = df.eval('0.5 * (giso_smass_err1 - giso_smass_err2)')
    samp['smass'] = sample(loc, scale, 2)

    # Teff
    loc = df.cks_steff
    scale = df.cks_steff_err
    samp['steff'] = sample(loc, scale, 3)

    # density
    loc = df.giso_srho
    scale = df.eval('0.5 * (giso_srho_err1 - giso_srho_err2)')
    samp['giso_srho'] = sample(loc, scale, 4)
    
    # Planet radius
    samp['gdir_prad'] = samp['dr25_ror'] * samp['gdir_srad'] * 109.245
    samp['giso_prad'] = samp['dr25_ror'] * samp['giso_srad'] * 109.245
    
    # Semi-major axis
    period = np.array(df.koi_period).reshape(-1, 1)
    samp['giso_sma'] = (samp['smass'] * (period/365)**2) ** 0.33

    # insolation flux
    samp['giso_sinc'] = ((samp['steff']/5778)**4
                         * (samp['giso_srad']/samp['giso_sma'])**2)

    # transit duration
    period = np.array(df[['koi_period']])
    samp['giso_tau0'] = (TAU_CONST * period**0.33 * samp['giso_srho']**-0.33)

    keys = 'dr25_period dr25_ror dr25_rho dr25_b dr25_tau dr25_dil gdir_prad giso_prad giso_sinc giso_sma giso_tau0'.split()
    for k in keys:
        kerr1 = k+'_err1'
        kerr2 = k+'_err2'
        df[k] = np.median(samp[k], axis=1)
        df[kerr1] = np.std(samp[k], axis=1)
        df[kerr2] = -1.0 * df[kerr1]

    df['dr25_fgraz'] = ((samp['dr25_b'] > 1.0).sum(axis=1).astype(float)
                        / samp['dr25_b'].shape[1])
    return df, samp


