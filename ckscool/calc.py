from collections import OrderedDict

import numpy as np
import pandas as pd

SEED = 0  # reproducability
np.random.RandomState(SEED)

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
    samp = dict(dr25_ror=[],dr25_rho=[],dr25_b=[],dr25_period=[]) 
    print("loading up {} chains".format(df))
    for i, row in df.iterrows():
        fn = '../Kepler-Radius-Ratio/test_trimmed-thinned_comp.hdf'
        chain = pd.read_hdf(fn,row.id_koicand)
        samp['dr25_ror'].append(chain['RD1'].sample(nsamp,replace=True))
        samp['dr25_rho'].append(chain['RHO'].sample(nsamp,replace=True))
        samp['dr25_b'].append(chain['BB1'].sample(nsamp,replace=True))
        samp['dr25_period'].append(chain['PE1'].sample(nsamp,replace=True))
        #if i % 100 == 0:
        print(i)
            
    for k in samp.keys():
        samp[k] = np.vstack(samp[k]).astype(float)
        
    samp['dr25_tau'] = (
        TAU_CONST
        * (1 - samp['dr25_b']**2)**0.5
        * samp['dr25_period']**0.33
        * samp['dr25_rho']**-0.33
    )

    def sample(loc, scale):
        loc = np.array(loc).reshape(-1, 1)
        scale = np.array(scale).reshape(-1, 1)
        return loc + scale * np.random.randn(1, nsamp)

    # stellar radius direct method
    loc = df.gdir_srad
    scale = df.eval('0.5 * (gdir_srad_err1 - gdir_srad_err2)')
    samp['gdir_srad'] = sample(loc, scale)

    # stellar radius grid method
    loc = df.giso_srad
    scale = df.eval('0.5 * (giso_srad_err1 - giso_srad_err2)')
    samp['giso_srad'] = sample(loc, scale)

    # stellar mass
    loc = df.giso_smass
    scale = df.eval('0.5 * (giso_smass_err1 - giso_smass_err2)')
    samp['smass'] = sample(loc, scale)

    # Teff
    loc = df.cks_steff
    scale = df.cks_steff_err
    samp['steff'] = sample(loc, scale)

    # density
    loc = df.giso_srho
    scale = df.eval('0.5 * (giso_srho_err1 - giso_srho_err2)')
    samp['giso_srho'] = sample(loc, scale)
    
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

    keys = 'dr25_ror dr25_rho dr25_b gdir_prad giso_prad giso_sinc giso_sma giso_tau0'.split()
    for k in keys:
        kerr1 = k+'_err1'
        kerr2 = k+'_err2'
        df[k] = np.median(samp[k], axis=1)
        df[kerr1] = np.std(samp[k], axis=1)
        df[kerr2] = -1.0 * df[kerr1]

    df['dr25_fgraz'] = ((samp['dr25_b'] > 1.0).sum(axis=1).astype(float)
                        / samp['dr25_b'].shape[1])
    return df
