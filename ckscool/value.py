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
    cand = ckscool.io.load_table('ckscool-targets-cuts',cache=2)
    star = cand.groupby('id_kic',as_index=False).nth(0)

    d['n-koi'] = "{}".format( len(cand) )
    d['n-star'] = "{}".format( len(star) )

    cand = cand[~cand.isany]
    star = cand.groupby('id_kic',as_index=False).nth(0)

    d['n-koi-pass'] = "{}".format(len(cand))
    d['n-star-pass'] = "{}".format(len(star))
    d['sparallax-med'] = "{:.1f}".format(star.eval('gaia2_sparallax').median())
    d['sparallax-ferr-med'] = "{:.1f}".format(star.eval('gaia2_sparallax_err / gaia2_sparallax').median() * 100)

    d['sdistance-med'] = "{:.0f}".format(star.eval('1000 / gaia2_sparallax').median())

    d['kmag-med'] = "{:.2f}".format(star.m17_kmag.median())
    d['kmag-err-med'] = "{:.2f}".format(star.m17_kmag_err.median())
    d['kepmag-med'] = "{:.1f}".format(star.kic_kepmag.median())
    d['kepmag-err-med'] = "{:.1f}".format(star.kic_kepmag.median())
    
    av_to_ak = (0.161+0.063)/(2.9197679+0.063) # extinction law
    
    df = ckscool.io.load_table('ckscool-stars',cache=1)
    df = df.groupby('id_starname',as_index=False).nth(0)

    d['av-med'] = "{:.03f}".format(df.gdir_avs.median())
    d['av-min'] = "{:.03f}".format(df.gdir_avs.min())
    d['av-max'] = "{:.03f}".format(df.gdir_avs.max())
    d['ak-med'] = "{:.03f}".format(float(d['av-med']) * av_to_ak)
    d['ak-min'] = "{:.03f}".format(float(d['av-min']) * av_to_ak)
    d['ak-max'] = "{:.03f}".format(float(d['av-max']) * av_to_ak)


    ferr = df.eval('0.5*(gdir_srad_err1 - gdir_srad_err2) / gdir_srad')
    d['srad-med'] = "{:.1f}".format(df.gdir_srad.median())
    d['srad-ferr-med'] = "{:.1f}".format(100*ferr.median())
    
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
    
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

