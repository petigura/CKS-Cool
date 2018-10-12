from collections import OrderedDict

import ckscool.io
import ckscool.observe
import numpy as np
from time import gmtime, strftime


def val_stat(return_dict=False):
    """
    Statistics of sample
    """
    d = OrderedDict()
    cand = ckscool.io.load_table('ckscool-cuts',cache=1)
    stars = cand.groupby('id_kic',as_index=False).first()

    d['n-koi'] = "{}".format( len(cand) )
    d['n-stars'] = "{}".format( len(stars) )

    cand = cand[~cand.isany]
    stars = cand.groupby('id_kic',as_index=False).first()

    d['n-koi-pass'] = "{}".format( len(cand))
    d['n-stars-pass'] = "{}".format( len(stars))

    obs = ckscool.observe.ObserveHIRES()
    d['hires-min-counts'] = obs.counts
    d['hires-min-snr'] = "{:.0f}".format(np.sqrt(obs.counts / 10.0)  * 45 )
    obs.label_observed('1980-01-01','2018-01-01')
    _d = obs.expected_observing_time()
    for key in 'hires-nstars-notobserved hires-time-notobserved'.split():
        d[key] = _d[key]

    # drop stars that have been observed prior to 2018
    sample = obs.sample
    sample.index = sample.id_kic
    sample = sample.drop(sample[sample['isobserved']].index)
    sample = sample.drop(['isobserved'],axis=1)
    obs.sample = sample
    obs.label_observed('2018-01-01','2019-01-10')
    _d = obs.expected_observing_time()

    d['date'] = strftime("%Y-%m-%d", gmtime())

    for key in 'hires-nstars-notobserved hires-time-notobserved'.split():
        d[key+'-2018'] = _d[key]


    obs = ckscool.observe.ObserveNIRC2()
    obs.label_observed()
    d['n-koi-sig-diluted'] = len(obs.sample[obs.sample.f17_sig_diluted])
    d['n-koi-k18-observed'] = len(obs.sample[obs.sample.k18_observed])

    _d = obs.expected_observing_time()
    for key in 'nirc2-nstars-notobserved nirc2-time-notobserved'.split():
        d[key] = _d[key]
    

    df = ckscool.io.load_table('ckscool-planets',cache=1)
    d['sparallax-med'] = "{:.1f}".format(df.eval('gaia2_sparallax').median())
    d['sparallax-ferr-med'] = "{:.1f}".format(df.eval('gaia2_sparallax_err / gaia2_sparallax').median() * 100)

    
    av_to_ak = (0.161+0.063)/(2.9197679+0.063) # extinction law
    
    d['av-med'] = "{:.03f}".format(df.gdir_avs.median())
    d['av-min'] = "{:.03f}".format(df.gdir_avs.min())
    d['av-max'] = "{:.03f}".format(df.gdir_avs.max())
    d['ak-med'] = "{:.03f}".format(float(d['av-med']) * av_to_ak)
    d['ak-min'] = "{:.03f}".format(float(d['av-min']) * av_to_ak)
    d['ak-max'] = "{:.03f}".format(float(d['av-max']) * av_to_ak)

    d['mk-med'] = "{:.02f}".format(df.m17_kmag.median())
    d['mk-err-med'] = "{:.02f}".format(df.m17_kmag_err.median())


    ferr = df.eval('0.5*(gdir_srad_err1 - gdir_srad_err2) / gdir_srad')
    d['srad-med'] = "{:.1f}".format(df.gdir_srad.median())
    d['srad-ferr-med'] = "{:.1f}".format(100*ferr.median())


    '''
    d['ror-med'] = "{:.1f}".format(100*df.koi_ror.median())
    d['ror-ferr-med'] = "{:.1f}".format(100*ferr.median())
    ferr= cut.eval('0.5*(gdir_prad_err1 - gdir_prad_err2) / gdir_prad')
    d['prad-ferr-med'] = "{:.1f}".format(100*ferr.median())
    d['prad-med'] = "{:.1f}".format(cut.gdir_prad.median())


    d['cks-gaia-distmod-err-med'] = "{:.2f}".format(cut.eval('gaia2_sparallax_err / gaia2_sparallax').median())

    dist = 1 / (1e-3  * df.gaia2_sparallax)
    mu = 5 * np.log10(dist) - 5 
    d['cks-gaia-distmod-med'] = "{:.2f}".format(mu.median())




    d['steff-med'] = "{:.0f}".format(df.cks_steff.median())



    # Properties of cks+gaia2 table
    df = cksgaia.io.load_table('cksgaia-planets',cache=1)
    cut = df
    '''


    cd = ckscool.cuts.cd
    for key in cd.keys():
        d[key] = cd[key]
    

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

