import pandas as pd
import numpy as np
from collections import OrderedDict
import cksgaia.io
from cksgaia.errors import equad

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
        steff = row.sm_steff
        steff_err = row.sm_steff_err
        steff_samp = rand.normal(row.sm_steff, steff_err, size=10000)

        # Radius ratio
        ror = row.koi_ror
        ror_err = 0.5 * (row.koi_ror_err1 - row.koi_ror_err2)
        ror_samp = rand.normal(row.koi_ror, ror_err, size=10000)

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

def table_statistics():
    d = OrderedDict()
    cache = 1
    df = cksgaia.io.load_table('cks', cache=cache)
    d['nstars-cks'] = "{}".format(len(df))

    df = cksgaia.io.load_table('cks+nea', cache=cache)
    d['ncand-cks'] = "{}".format(len(df))

    # Model dep errors
    for k in 'smass srad-dw srad-gi sage'.split():
        d['equad-' + k] = "{:.0f}\%".format(equad[k] * 100)

    # Error summaries with and without floor
    for table in ['iso', 'iso-floor']:
        df = cksgaia.io.load_table(table, cache=cache)
        for k in 'smass srad slogage'.split():
            if k == 'slogage':
                err = df['giso_' + k + '_err1']
                fmt = ".2f"
                unit = "~dex"
            else:
                fmt = ".1f"
                err = df['giso_' + k + '_frac_err'] * 100
                unit = "\%"
            for p in [5, 50, 95]:
                fmt2 = "{:%s}%s" % (fmt, unit)
                percentile = fmt2.format(np.percentile(err, p))
                d["{}-{}-err-{:02d}".format(table, k, p, unit)] = percentile

    df = cksgaia.io.load_table('cks+nea+iso-floor', cache=cache)
    d['cks-rp-frac-err-median'] = "{:.0f}\%".format(
        (df['gdir_prad_err1'] / df['gdir_prad']).median() * 100)
    d['cks-sinc-frac-err-median'] = "{:.0f}\%".format(
        (df['giso_insol_err1'] / df['giso_insol']).median() * 100)
    d['cks-sma-frac-err-median'] = "{:.1f}\%".format(
        (df['giso_sma_err1'] / df['giso_sma']).median() * 100)

    for k, v in d.iteritems():
        print r"{{{}}}{{{}}}".format(k, v)


def kdeslice(x, xvalue, kde):
    """
    Return a slice through a KDE at the specificed x value.
    :param x: array of x values corresponding to the columns in the KDE array
    :param xvalue: x location to slice
    :param kde: 2D KDE to be sliced
    :return: 1d array of KDE values

    """

    pos = np.argmin(np.abs(x - xvalue))
    print pos
    sl = kde[pos,:].flatten()

    return sl


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values - average) ** 2, weights=weights)
    return (average, np.sqrt(variance / len(values)))


def average_in_box(sample, box, col1='gdir_prad', col2='koi_period', logparam=True):
    radlim_low = box[0][1]
    radlim_high = box[1][1]
    plim_low = box[0][0]
    plim_high = box[1][0]
    q = sample.query('@radlim_low <= gdir_prad < @radlim_high & @plim_low < koi_period <= @plim_high')
    n = len(q)

    q = q[np.isfinite(q[col1]) & np.isfinite(q[col2]) & np.isfinite(q['weight'])]

    if logparam:
        logcol1 = np.log10(q[col1])
        logcol2 = np.log10(q[col2])
        avg_rad, err_rad = weighted_avg_and_std(logcol1, q['weight'])
        avg_per, err_per = weighted_avg_and_std(logcol2, q['weight'])
        return [(10 ** avg_per, 10 ** avg_rad), (err_per * (10 ** avg_per), err_rad * (10 ** avg_rad))]
    else:
        avg_rad, err_rad = weighted_avg_and_std(q[col1], q['weight'])
        avg_per, err_per = weighted_avg_and_std(q[col2], q['weight'])
        return [(avg_per, avg_rad), (err_per, err_rad)]


