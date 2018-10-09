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

