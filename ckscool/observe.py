"""
Module to plan CKS-Cool observations.

- Calculate needed observing time
- Calculate total open shutter time for observations (anything observed before 2018) doesn't count

- Generate HIRES scripts
- Generate NIRC2 scripts

"""

import ckscool.io
import cpsutils.kbc
import astropy.time
import pandas as pd
import cpsutils.hires
import cpsutils.hires.exposure
from collections import OrderedDict
import numpy as np

class Observe(object):
    def __init__(self):
        df = ckscool.io.load_table('ckscool-cuts')
        df = df[~df.isany]
        df = df.groupby('id_kic',as_index=False).first()
        self.sample = df

class ObserveNIRC2(Observe):
    texp = 10.0/60.0 # ten minutes per target.
    def label_observed(self, start, stop):
        """
        Label stars as observed fold in Adam Kraus's
        """
        k18 = ckscool.io.load_table('nrm-previous')
        #start_jd = astropy.time.Time(start,format='iso').jd
        #stop_jd = astropy.time.Time(stop,format='iso').jd
        #prev = prev[prev.jd.between(start_jd,stop_jd)]
        k18['k18_observed'] = True
        self.sample = pd.merge(self.sample,k18,on='id_koi',how='left')
        self.sample['k18_observed'] = self.sample['k18_observed'].fillna(False)
        
        f17 = ckscool.io.load_table('furlan17-tab9')
        f17 = f17['id_koi f17_rcorr_avg'.split()]
        f17['f17_diluted'] = True
        f17['f17_sig_diluted'] = (f17.f17_rcorr_avg - 1) > 0.05
        self.sample = pd.merge(self.sample, f17, on='id_koi',how='left')

        # Is observed by Keck?
        f17 = ckscool.io.load_table('furlan17-tab2')
        f17 = f17['id_koi f17_ao_obs'.split()]
        #f17 = f17[f17.f17_ao_obs.str.contains('Keck')]
        #f17['f17_observed_keck'] = True
        self.sample = pd.merge(self.sample, f17, on='id_koi',how='left')
        self.sample['f17_ao_obs'] = self.sample['f17_ao_obs'].fillna('none')
        self.sample['f17_sig_diluted'] = self.sample['f17_sig_diluted'].fillna(False)

    def expected_observing_time(self):
        """
        Return amount of on sky NIRC2 LGS time we would need.
        """

        # Only count stars that haven't been observed
        df = self.sample
        vmag = self.sample.m17_kepmag + 0.4 # 
        df['texp'] = self.texp

        d = OrderedDict()
        d['nirc2-nstars'] = len(df)
        d['nirc2-time'] = (df.texp).sum()
        
        bobs = ~df.k18_observed & ~df.f17_sig_diluted & ~df.f17_ao_obs.str.contains('Keck|Gem')
        df = df[bobs]
        d['nirc2-nstars-notobserved'] = len(df)
        d['nirc2-time-notobserved'] = np.round((df.texp).sum(),1)
        return d

    def remaining_observing_time():
        """
        Return amount of on sky NIRC2 LGS time we would need.
        """
        pass 

class ObserveHIRES(Observe):
    tover = 1.0/60. # one minute of overhead perstar
    counts = 3 # target time on the exposure meter

    def label_observed(self, start, stop):
        """
        Label stars as observed
        """
        prev = cpsutils.kbc.loadkbc()
        prev = prev[
            prev.name.str.contains('^K\d{5}$|^CK\d{5}$') & (prev['type']=='t')
        ]

        prev['id_koi'] = prev.name.str.slice(start=-5).astype(int)
        prev['jd'] += 2440000
        start_jd = astropy.time.Time(start,format='iso').jd
        stop_jd = astropy.time.Time(stop,format='iso').jd
        prev = prev[prev.jd.between(start_jd,stop_jd)]
        prev['isobserved'] = True
        self.sample = pd.merge(self.sample,prev,on='id_koi',how='left')
        self.sample['isobserved'] = self.sample['isobserved'].fillna(False)

    def expected_observing_time(self):
        """
        Expected amount of HIRES time
        """
        # Only count stars that haven't been observed
        df = self.sample
        vmag = self.sample.m17_kepmag + 0.4 # 
        
        df['texp'] = cpsutils.hires.exposure.exposure_time(vmag, self.counts) 
        df['texp'] /= 3600

        d = OrderedDict()
        d['hires-nstars'] = len(df)
        d['hires-time'] = (df.texp + self.tover).sum()

        df = df[~df.isobserved]
        d['hires-nstars-notobserved'] = len(df)
        d['hires-time-notobserved'] = np.round((df.texp + self.tover).sum(),1)
        return d
        
