"""
Module for the reading and writing of gaia data 
"""

import numpy as np
import pandas as pd
import ckscool.io

def foverlap(x0):
    step = 0.02
    x = np.arange(-15,15,step)
    area = step*step
    x,y = np.meshgrid(x,x)

    FWHM = 2.5
    sigma = FWHM / 2.355

    A = 1.0 / 2 / np.pi / sigma**2
    a = A * np.exp(-((x - x0)**2 / 2 / sigma**2 + y**2 / 2 / sigma**2))
    a.sum() * area
    R = 4
    return a[np.sqrt(x**2 + y**2) < R].sum() * area

xi = np.arange(0,8,0.2)
overlap = np.array([foverlap(x) for x in xi])

def read_xmatch_gaia2(fn):
    df = pd.read_csv(fn)
    namemap = {
        'source_id':'id_gaia2',
        'ra':'ra', 
        'dec':'dec',
        'parallax':'sparallax', 
        'parallax_error':'sparallax_err', 
        'phot_g_mean_flux':'gflux',
        'phot_g_mean_flux_error':'gflux_err',
        'phot_g_mean_mag':'gmag',

        'phot_bp_mean_flux':'bpflux',
        'phot_bp_mean_flux_error':'bpflux_err',
        'phot_bp_mean_mag':'bpmag',
        'phot_rp_mean_flux':'rpflux',
        'phot_rp_mean_flux_error':'rpflux_err',
        'phot_rp_mean_mag':'rpmag',
        'parallax_over_error':'sparallax_over_err',
        'astrometric_excess_noise':'astrometric_excess_noise',
        'astrometric_excess_noise_sig':'astrometric_excess_noise_sig',
        'id_kic':'id_kic',
        'dist':'angdist',
        'ruwe':'ruwe',
    }
    
    df['steff'] = df['teff_val']
    df['steff_err1'] = df.eval('teff_percentile_upper - teff_val')
    df['steff_err2'] = df.eval('teff_percentile_lower - teff_val')
    df['srad'] = df['radius_val']
    df['srad_err1'] = df.eval('radius_percentile_upper - radius_val')
    df['srad_err2'] = df.eval('radius_percentile_lower - radius_val')
    df = df.rename(columns=namemap)
    df['angdist'] *= 60*60
    cols = namemap.values() + 'steff steff_err1 steff_err2 srad srad_err1 srad_err2 '.split()
    df = df[cols]
    df = ckscool.io.add_prefix(df, 'gaia2_')
    return df


def xmatch_gaia2(df, gaia):
    """
    Crossmatch the sources in Gaia 2

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR2 table
    """

    assert len(df.id_kic)==len(df.id_kic.drop_duplicates()), "No duplicate stars"
    gaia['id_gaia2'] = gaia['id_gaia2'].astype(str)

    # Estimate Kmag and kmag flux from gaia
    p = [1.7,0]
    bp_rp = gaia.eval('gaia2_bpmag-gaia2_rpmag')
    g_k = np.polyval(p,bp_rp)
    gaia['gaia2_kmag'] = gaia['gaia2_gmag'] - g_k
    gaia['gaia2_kflux'] = 2e10*10**(-0.4*gaia.eval('gaia2_kmag'))

    foverlap = np.interp(gaia['gaia2_angdist'], xi, overlap)

    gaia['gaia2_gflux_aper'] = gaia['gaia2_gflux'] * foverlap
    gaia['gaia2_kflux_aper'] = gaia['gaia2_kflux'] * foverlap

    # just want the stars
    m = pd.merge(df,gaia,on='id_kic',how='left')
    m['id_gaia2']= m['id_gaia2'].fillna('-99').astype(np.int64)
    ndf = len(df)
    print "max(gaia_angdist) = {} (arcsec)".format(m['gaia2_angdist'].max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    m = m.set_index('id_kic')
    g = m.groupby('id_kic')
    m['gaia2_gflux_sum'] = g['gaia2_gflux_aper'].sum()
    m['gaia2_kflux_sum'] = g['gaia2_kflux_aper'].sum()
    m['gaia2_absdiff_gkep'] = np.abs(m['gaia2_gmag'] - m['m17_kepmag'])
    m['gaia2_n_8arcsec'] = g.size()
    m = m.reset_index()

    # Match candidate
    #mbest = m.query('gaia2_angdist < 0.5 and -0.1 < m17_kepmag - gaia2_gmag < 0.15')
    mbest = m.query('gaia2_angdist < 1 and -0.1 < m17_kepmag - gaia2_gmag < 0.15')
    mbest['gaia2_n_1arcsec'] = mbest.groupby('id_kic').size()

    print "{} gaia sources within 1 arcsec of {} target sources".format(
        len(mbest),mbest.id_kic.drop_duplicates().count()
    )

    mbest = mbest.sort_values(by=['id_kic','gaia2_absdiff_gkep'])
    g = mbest.groupby('id_kic',as_index=False)
    mbest['gaia2_n_1arcsec'] = g.size()
    mbest = g.nth(0) 
    mbest['gaia2_gflux_ratio'] = mbest.eval('gaia2_gflux_sum / gaia2_gflux')
    mbest['gaia2_kflux_ratio'] = mbest.eval('gaia2_kflux_sum / gaia2_kflux')
    return mbest, m 
