import pandas as pd
import numpy as np
from collections import OrderedDict

SEED = 0  # reproducability
rand = np.random.RandomState(SEED)

albedo = 0.30

def update_planet_parameters(df):
    newkeys = [
        'gdir_prad','gdir_prad_err1','gdir_prad_err2',
        'giso_sinc','giso_sinc_err1','giso_sinc_err2',
        'giso_sma','giso_sma_err1','giso_sma_err2',
    ]

    for key in newkeys:
        df[key] = np.nan

    for i, row in df.iterrows():

        # stellar radius
        srad = row.gdir_srad
        srad_err = 0.5 * (row.gdir_srad_err1 - row.gdir_srad_err2)
        srad_samp = rand.normal(row.gdir_srad, srad_err, size=10000)

        # stellar mass
        smass = row.giso_srad
        smass_err = 0.5 * (row.giso_smass_err1 - row.giso_smass_err2)
        smass_samp = rand.normal(row.giso_smass, smass_err, size=10000)

        # Teff
        steff = row.cks_steff
        steff_err = row.cks_steff_err
        steff_samp = rand.normal(row.cks_steff, steff_err, size=10000)

        # Radius ratio
        ror = row.dr25_ror
        ror_err = 0.5 * (row.dr25_ror_err1 - row.dr25_ror_err2)
        ror_samp = rand.normal(row.dr25_ror, ror_err, size=10000)

        # Planet radius
        prad_samp = ror_samp * srad_samp * 109.245

        # Semi-major axis
        sma_samp = (smass_samp * (row.koi_period / 365.) ** 2.) ** (1. / 3.)

        # insolation flux
        sinc_samp = (steff_samp / 5778.) ** 4.0 * (srad_samp / sma_samp) ** 2.0


        df.ix[i,'gdir_prad'] = np.median(prad_samp)
        df.ix[i,'gdir_prad_err1'] = np.std(prad_samp)
        df.ix[i,'gdir_prad_err2'] = -1.0 * np.std(prad_samp)
        df.ix[i,'giso_sinc'] = np.median(sinc_samp)
        df.ix[i,'giso_sinc_err1'] = np.std(sinc_samp)
        df.ix[i,'giso_sinc_err2'] = -1.0 * np.std(sinc_samp)
        df.ix[i,'giso_sma'] = np.median(sma_samp)
        df.ix[i,'giso_sma_err1'] = np.std(sma_samp)
        df.ix[i,'giso_sma_err2'] = -1.0 * np.std(sma_samp)

    return df
