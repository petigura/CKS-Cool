from collections import OrderedDict

import ckscool.io
import ckscool.cuts
import numpy as np
from time import gmtime, strftime
from astropy.time import Time
import pandas as pd


def val_stat(return_dict=False):
    """
    Statistics of sample
    """
    d = OrderedDict()


    cand = ckscool.io.load_table('planets-cuts1',cache=2)
    star = cand.groupby('id_kic',as_index=False).nth(0)
    field = ckscool.io.load_table('field-cuts',cache=2)

    d['n-koi'] = "{}".format( len(cand) )
    d['n-star'] = "{}".format( len(star) )

    d['n-star-field'] = "{}".format( len(field) )
    d['n-star-field-pass0'] = "{}".format( len(field[~field.isany]) )

    df = ckscool.io.load_table('planets-cuts1')
    df = ckscool.cuts.occur.add_cuts(df, field.cuttypes, 'field')
    d['n-cand'] = "{}".format( len(df) )

    cut = df[~df.isany]
    d['n-cand-pass0'] = "{}".format( len(cut ) )
    d['n-star-pass0'] = "{}".format( len(cut.id_koi.drop_duplicates()) )

    df = ckscool.io.load_table('planets-cuts1')
    cuttypes1 = df.cuttypes
    cand = df[~df.isany]
    d['n-cand-pass1'] = "{}".format( len(cand) )
    star1 = cand.groupby('id_koi',as_index=False).nth(0)
    d['n-star-pass1'] = "{}".format( len(star1) )

    df = ckscool.io.load_table('planets-cuts2+iso')
    star2 = df.groupby('id_koi',as_index=False).nth(0)
    m = pd.merge(star1,star2[['id_koi','cks_sprov']],how='left')
    m['cks_sprov'].fillna('none',inplace=True)
    s = m.groupby('cks_sprov').size()
    for sprov, count in s.iteritems():
        d['n-star-pass1-{}'.format(sprov,count)] = count

    d['n-star-pass1-newobs'] = s['smemp'] + s['smsyn']
    d['n-star-pass1-lt4700'] = len(m.query('ber18_steff < 4700'))
    
    df[df[["is"+c for c in cuttypes1]].sum(axis=1)==0]
    d['sparallax-med'] = "{:.1f}".format(star.eval('gaia2_sparallax').median())
    d['sparallax-ferr-med'] = "{:.1f}".format(star.eval('gaia2_sparallax_err / gaia2_sparallax').median() * 100)
    d['sdistance-med'] = "{:.0f}".format(star.eval('1000 / gaia2_sparallax').median())

    d['kmag-med'] = "{:.2f}".format(star.m17_kmag.median())
    d['kmag-err-med'] = "{:.2f}".format(star.m17_kmag_err.median())
    d['kepmag-med'] = "{:.1f}".format(star.m17_kepmag.median())
    d['kepmag-err-med'] = "{:.1f}".format(star.m17_kepmag.median())
    
    av_to_ak = (0.161+0.063)/(2.9197679+0.063) # extinction law

    df = ckscool.io.load_table('planets+iso',cache=1)
    d['av-med'] = "{:.03f}".format(df.gdir_avs.median())
    d['av-min'] = "{:.03f}".format(df.gdir_avs.min())
    d['av-max'] = "{:.03f}".format(df.gdir_avs.max())
    d['ak-med'] = "{:.03f}".format(float(d['av-med']) * av_to_ak)
    d['ak-min'] = "{:.03f}".format(float(d['av-min']) * av_to_ak)
    d['ak-max'] = "{:.03f}".format(float(d['av-max']) * av_to_ak)


    ferr = df.eval('0.5*(gdir_srad_err1 - gdir_srad_err2) / gdir_srad')
    d['srad-med'] = "{:.1f}".format(df.gdir_srad.median())
    d['srad-ferr-med'] = "{:.1f}".format(100*ferr.median())

    '''
    df = ckscool.io.load_table('ckscool-stars-cuts',cache=0)
    cuttypes = df.cuttypes
    df = df.groupby('id_koi').nth(0)

    for cut in cuttypes:
        d['cut-nstars-'+cut] = df['is'+cut].sum()

    d['cut-nstars-any'] = df['isany'].sum()

    obs = ckscool.observe.ObserveNIRC2()
    obs.label_observed()

    d['nirc2-k18-observed'] = obs.sample['k18_observed'].sum()
    d['nirc2-f17-observed'] = obs.sample['f17_observed'].sum()
    d['nirc2-p18-observed'] = obs.sample['p18_observed'].sum()
    d['nirc2-tot-observed'] = obs.sample['isobserved'].sum()


    df = ckscool.io.load_table('cksgaia-planets')
    df = df.groupby('id_koi',as_index=False).nth(0)
    d['cks1-sdistance-med'] = "{:.0f}".format(df.query('kic_kepmag < 14.2').eval('1000 / gaia2_sparallax').median())

    cd = ckscool.cuts.cd
    for key in cd.keys():
        d[key] = cd[key]
    

    '''
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

def val_sample(return_dict=False):
    d = OrderedDict()

    m = ckscool.io.load_table('DR1+DR2+CXM-overlap')
    d['nstars dr1'] = m.in_dr1.sum()
    d['nstars dr2'] = m.in_dr2.sum()
    d['nstars cxm'] =  m.in_cxm.sum()
    d['nstars dr1 & cxm'] = (m.in_dr1 & m.in_cxm).sum()
    d['nstars dr1 & ~cxm'] = (m.in_dr1 & ~m.in_cxm).sum()
    d['nstars dr2 & cxm'] = (m.in_dr2 & m.in_cxm).sum()
    d['nstars dr2 & ~cxm'] = (m.in_dr2 & ~m.in_cxm).sum()
    d['nstars ~dr1 & ~dr2 & cxm'] = (~m.in_dr1 & ~m.in_dr2 & m.in_cxm).sum()
    d['nstars (dr1 | dr2) & cxm'] = ((m.in_dr1 | m.in_dr2) & m.in_cxm).sum()
    d['nstars dr1 & dr2'] = (m.in_dr1 & m.in_dr2).sum()
    d['nstars dr1 | dr2'] = (m.in_dr1 | m.in_dr2).sum()
    d['nstars dr1 | dr2 | cxm'] = (m.in_dr1 | m.in_dr2 | m.in_cxm).sum()

    d['nstars dr2 & cxm & pre-2018'] = (
        m.in_dr2
        & m.in_cxm
        & (m.obs_bjd < Time('2018-01-01',format='iso').jd)
    ).sum()

    d['nstars dr2 & cxm & post-2018'] = (
        m.in_dr2
        & m.in_cxm
        & (m.obs_bjd > Time('2018-01-01',format='iso').jd)
    ).sum()

    d['nstars dr2 & post-2018'] = (
        m.in_dr2
        & (m.obs_bjd > Time('2018-01-01',format='iso').jd)
    ).sum()

    b = ( m.in_dr2
          & m.in_cxm
          & (m.obs_bjd < Time('2018-01-01',format='iso').jd)
          & (m.obs_counts < 1500) )

    d['nstars dr2 & cxm & pre-2018 counts < 1500'] = b.sum()
    kois = ", ".join(m[b].id_koi.astype(int).astype(str))
    d['kois dr2 & cxm & pre-2018 counts < 1500'] = kois

    table = 'planets-cuts1'
    df = ckscool.io.load_table(table,cache=2)
    i = 0 
    bpass = np.zeros(len(df))
    for cuttype in df.cuttypes:
        key = 'is'+cuttype
        obj = ckscool.cuts.occur.get_cut(cuttype)
        cut = obj(df,table)
        bpass += df[key].astype(int)
        dfcut = df[bpass==0]
        key = 'nplanets {}-cut-{}-{}'.format(table,i,cuttype)
        d[key] = len(dfcut)

        key = 'nstars {}-cut-{}-{}'.format(table,i,cuttype)
        d[key] = len(dfcut.id_kic.drop_duplicates())
        i+=1


    table = 'field-cuts'
    # needs to be freshly generated to get cuttypes
    df = ckscool.io.load_table(table,cache=2) 
    i = 0 
    bpass = np.zeros(len(df))
    for cuttype in df.cuttypes:
        key = 'is'+cuttype
        obj = ckscool.cuts.occur.get_cut(cuttype)
        cut = obj(df,table)
        bpass += df[key].astype(int)
        dfcut = df[bpass==0]
        key = 'nstars {}-cut-{}-{}'.format(table,i,cuttype)
        d[key] = len(dfcut.id_kic.drop_duplicates())
        i+=1

    table = 'planets-cuts2'
    # needs to be freshly generated to get cuttypes
    df = ckscool.io.load_table(table,cache=2) 
    i = 0 
    bpass = np.zeros(len(df))
    for cuttype in df.cuttypes:
        key = 'is'+cuttype
        obj = ckscool.cuts.occur.get_cut(cuttype)
        cut = obj(df,table)
        bpass += df[key].astype(int)
        dfcut = df[bpass==0]
        key = 'nstars {}-cut-{}-{}'.format(table,i,cuttype)
        d[key] = len(dfcut.id_kic.drop_duplicates())
        key = 'nplanets {}-cut-{}-{}'.format(table,i,cuttype)
        d[key] = len(dfcut.id_kic)
        i+=1

    key = 'nstars {}-cut-all'.format(table)
    d[key] = len(dfcut.id_kic.drop_duplicates())
    key = 'nplanets {}-cut-all'.format(table)
    d[key] = len(dfcut.id_kic)

    
    df = ckscool.io.load_table('star')
    s = df.groupby('cks_sprov').size()
    d['nstars smemp'] = s.emp
    d['nstars smsyn'] = s.syn

    faint  = df.query('m17_kepmag>14.2 ')
    bright  = df.query('m17_kepmag<14.2 ')
    d['kmag-err-med kepmag<14.2'] = "{:.2f}".format(
        bright.m17_kmag_err.median()
    )
    d['kmag-err-med kepmag>14.2'] = "{:.2f}".format(
        faint.m17_kmag_err.median()
    )

    _eval = '100 * gaia2_sparallax_err / gaia2_sparallax'
    d['parallax-ferr-med kepmag<14.2'] = "{:.1f}".format(
        bright.eval(_eval).median()
    )
    d['parallax-ferr-med kepmag>14.2'] = "{:.1f}".format(
        faint.eval(_eval).median()
    )

    av_to_ak = (0.161+0.063)/(2.9197679+0.063) # extinction law
    d['ak-med kepmag<14.2'] = "{:.2f}".format(bright.gdir_avs.median()
                                              * av_to_ak)
    d['ak-med kepmag>14.2'] = "{:.2f}".format(faint.gdir_avs.median()
                                              * av_to_ak)

    d['ak-max kepmag<14.2'] = "{:.2f}".format(bright.gdir_avs.max()
                                              * av_to_ak)
    d['ak-max kepmag>14.2'] = "{:.2f}".format(faint.gdir_avs.max()
                                              * av_to_ak)


    _eval = '100 * 0.5*(gdir_srad_err1 - gdir_srad_err2) / gdir_srad'
    d['gdir_srad-ferr-med kepmag<14.2'] = "{:.1f}".format(
        bright.eval(_eval)
        .median()
    )
    d['gdir_srad-ferr-med kepmag>14.2'] = "{:.1f}".format(
        faint.eval(_eval)
        .median()
    )

    
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)


    if return_dict:
        return d

    return lines

def val_fit(return_dict=False):
    d = OrderedDict()
    keys = [
        'fitdetected_sn-m',
        'fitdetected_se-m',
        'fitdetected_sn-mm',
        'fitdetected_se-mm',
        'fitdetected_sn-mma',
        'fitdetected_se-mma'
    ]
    for key in keys:
        fitter = ckscool.io.load_object(key,cache=1)
        mode = fitter.mode
        samp = fitter.samples
        samp['R0'] = 10**(fitter.samples.logR_0)
        desc = samp.describe()
        d[mode+'-R0'] = "{:.2f} \pm {:.2f}".format(*desc.loc[['mean','std'],'R0'].tolist())
        d[mode+'-alpha-mass'] = "{:.2f} \pm {:.2f}".format(*desc.loc[['mean','std'],'alpha_mass'].tolist())
        d[mode+'-alpha-met'] = "{:.2f} \pm {:.2f}".format(*desc.loc[['mean','std'],'alpha_met'].tolist())
        d[mode+'-alpha-age'] = "{:.2f} \pm {:.2f}".format(*desc.loc[['mean','std'],'alpha_age'].tolist())

    keys = [
        'mps_size-se',
        'mps_size-sn',
    ]

    for key in keys:
        mps = ckscool.io.load_object(key,cache=1)
        d[key+'-alpha'] = "{:.2f} \pm {:.2f}".format(np.mean(mps.coeff[1]), np.std(mps.coeff[1]))
        d[key+'-R0'] = "{:.2f} \pm {:.2f}".format(np.mean(10**mps.coeff[0]), np.std(10**mps.coeff[0]))
    
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)


    if return_dict:
        return d

    return lines

