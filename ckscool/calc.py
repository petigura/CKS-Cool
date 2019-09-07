from collections import OrderedDict

import numpy as np
import pandas as pd

SEED = 0  # reproducability
np.random.RandomState(SEED)

def update_planet_parameters(df):
    """Update planet parameters

    Args:
         df: DataFrame

    Returns
        pandas.DataFrame: with following columns added

         -    gdir_prad
         -    gdir_prad_err1
         -    gdir_prad_err2
         -    giso_prad
         -    giso_prad_err1
         -    giso_prad_err2
         -    giso_sinc
         -    giso_sinc_err1
         -    giso_sinc_err2
         -    giso_sma
         -    giso_sma_err1
         -    giso_sma_err2

    """
    nsamp = 1000
    def sample(loc, scale):
        loc = np.array(loc).reshape(-1, 1)
        scale = np.array(scale).reshape(-1, 1)
        return loc + scale * np.random.randn(1, nsamp)

    # stellar radius direct method
    loc = df.gdir_srad
    scale = df.eval('0.5 * (gdir_srad_err1 - gdir_srad_err2)')
    gdir_srad_samp = sample(loc, scale)

    # stellar radius grid method
    loc = df.giso_srad
    scale = df.eval('0.5 * (giso_srad_err1 - giso_srad_err2)')
    giso_srad_samp = sample(loc, scale)

    # stellar mass
    loc = df.giso_smass
    scale = df.eval('0.5 * (giso_smass_err1 - giso_smass_err2)')
    smass_samp = sample(loc, scale)

    # Teff
    loc = df.cks_steff
    scale = df.cks_steff_err
    steff_samp = sample(loc, scale)

    # Radius ratio
    loc = df.dr25_ror
    scale = df.eval('0.5 * (dr25_ror_err1 - dr25_ror_err2)')
    empty = scale.isnull()
    scale = scale.fillna(0)
    ror_samp = sample(loc, scale)
    ror_samp[empty, :] = np.nan # fill with null

    # Planet radius
    gdir_prad_samp = ror_samp * gdir_srad_samp * 109.245
    giso_prad_samp = ror_samp * giso_srad_samp * 109.245

    # Semi-major axis
    period = np.array(df.koi_period).reshape(-1, 1)
    sma_samp = (smass_samp * (period/365)**2) ** 0.33

    # insolation flux
    sinc_samp = (steff_samp/5778)**4 * (giso_srad_samp/sma_samp)**2

    df['gdir_prad'] = np.median(gdir_prad_samp, axis=1)
    df['gdir_prad_err1'] = np.std(gdir_prad_samp, axis=1)
    df['gdir_prad_err2'] = -1.0 * df['gdir_prad_err1']
    df['giso_prad'] = np.median(giso_prad_samp)
    df['giso_prad_err1'] = np.std(giso_prad_samp)
    df['giso_prad_err2'] = -1.0 * df['giso_prad_err1']
    df['giso_sinc'] = np.median(sinc_samp, axis=1)
    df['giso_sinc_err1'] = np.std(sinc_samp, axis=1)
    df['giso_sinc_err2'] = -1.0 * df['giso_sinc_err1']
    df['giso_sma'] = np.median(sma_samp, axis=1)
    df['giso_sma_err1'] = np.std(sma_samp, axis=1)
    df['giso_sma_err2'] = -1.0 * df['giso_sma_err1']
    return df
