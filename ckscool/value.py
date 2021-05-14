from collections import OrderedDict

import ckscool.io
import ckscool.cuts
import numpy as np
from time import gmtime, strftime
from astropy.time import Time
import pandas as pd

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

    b = ~m.in_dr1 & ~m.in_dr2 & m.in_cxm
    kois = ", ".join(m[b].id_koi.astype(int).astype(str))
    d['kois ~dr1 & ~dr2 & cxm'] = kois

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

    # gradient slope
    smass = [0.5,0.7,1.0,1.4]
    for i in range(3):
        smass1, smass2 = smass[i],smass[i+1]
        key = 'grad-per-prad_smass={0}-{1}'.format(
            smass1, smass2
        )            
        grad = ckscool.io.load_object(key,cache=1)
        m = grad.grad_chain[:,0]
        Rp_10 = 10**grad.grad_chain[:,1]
        
        d[key+'m'] = "{:.2f} \pm {:.2f}".format(m.mean(),m.std())
        d[key+'Rp10'] = "{:.2f} \pm {:.2f}".format(Rp_10.mean(),Rp_10.std())


        
    smass = [0.5,0.7,1.0,1.4]
    for i in range(3):
        smass1, smass2 = smass[i],smass[i+1]
        prad = [1.0,1.7,4.0]
        for j in range(2):
            prad1, prad2 = prad[j],prad[j+1]
            key = 'fitper_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            chain = fit.res.flatchain
            chain['x0'] = 10**chain['logx0']
            chain['f'] = 10**chain['logf']
            q = chain.quantile([0.14,0.5,0.86])
            q.loc['up'] = q.loc[0.86] - q.loc[0.50]
            q.loc['lo'] = q.loc[0.14] - q.loc[0.50]
            d[key+'-logf'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'logf'].tolist())
            d[key+'-f'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'f'].tolist())
            d[key+'-k1'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'k1'].tolist())
            #d[key+'-k2'] = "{:.1f}^{{ {:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'k2'].tolist())
            d[key+'-logx0'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'logx0'].tolist())
            d[key+'-x0'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'x0'].tolist())

            
    smass = [0.5,0.7,1.0,1.4]
    for i in range(3):
        smass1, smass2 = smass[i],smass[i+1]
        prad = [1.0,1.7,4.0]
        for j in range(2):
            prad1, prad2 = prad[j],prad[j+1]
            key = 'fitsinc_smass={}-{}-prad={}-{}'.format(
                smass1, smass2, prad1, prad2,
            )
            fit = ckscool.io.load_object(key,cache=1)
            chain = fit.res.flatchain
            chain['x0'] = 10**chain['logx0']
            chain['f'] = 10**chain['logf']
            q = chain.quantile([0.14,0.5,0.86])
            q.loc['up'] = q.loc[0.86] - q.loc[0.50]
            q.loc['lo'] = q.loc[0.14] - q.loc[0.50]
            d[key+'-logf'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'logf'].tolist())
            d[key+'-f'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'f'].tolist())
            #d[key+'-k1'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'k1'].tolist())
            d[key+'-k2'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'k2'].tolist())
            d[key+'-logx0'] = "{:.1f}^{{ +{:.1f} }}_{{ {:.1f} }}".format(*q.loc[[0.5,'up','lo'],'logx0'].tolist())
            d[key+'-x0'] = "{:.0f}^{{ +{:.0f} }}_{{ {:.0f} }}".format(*q.loc[[0.5,'up','lo'],'x0'].tolist())
            
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

        
    if return_dict:
        return d

    return lines

def val_grad(return_dict=False):
    d = OrderedDict()
    smass = ['0.5-0.7','0.7-1.0','1.0-1.4','0.5-1.4']
    outdir = 'gapfitting-comparison/gradient-column/'
    for _smass in smass:
        fn = '{}/grad-per-prad-det_smass={}.csv'.format(outdir, _smass)
        key = 'per-prad-det_smass={}'.format(_smass)
        df = pd.read_csv(fn)
        df['Rp10'] = 10**df.y0
        q = df.quantile([0.14,0.5,0.86])
        q.loc['up'] = q.loc[0.86] - q.loc[0.50]
        q.loc['lo'] = q.loc[0.14] - q.loc[0.50]
        d[key+'-Rp10'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'Rp10'].tolist())
        d[key+'-m'] = "{:.2f}^{{ +{:.2f} }}_{{ {:.2f} }}".format(*q.loc[[0.5,'up','lo'],'m'].tolist())
        
        
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)
        
    if return_dict:
        return d

    return lines

