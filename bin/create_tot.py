#!/usr/bin/env python
import pandas as pd
import glob
import os

sources='cks1 smsyn smemp'.split()
suffixes='direct grid-parallax-yes grid-parallax-no'.split()
modes='direct grid grid'.split()

for source in sources:
    for suffix,mode in zip(suffixes,modes):
        fmt = dict(source=source,suffix=suffix,mode=mode)
        # all stars
        df = pd.read_csv('data/isoclassify-{source:}-{suffix:}.csv'.format(**fmt))
        # existing output
        fmt['outdir'] = 'isoclassify/{source:}/{suffix:}'.format(**fmt)
        files = '{outdir:}/*/*.csv'.format(**fmt)
        files = glob.glob(files)
        files = pd.DataFrame(files,columns=['starname'])
        files['id_starname'] = files.starname.apply(lambda x : x.split('/')[-2])
        files['existing'] = True
        df = pd.merge(df,files,how='left')
        df['existing'].fillna(False,inplace=True)
        df = df[~df['existing']]
        csvfn = 'isoclassify-{source:}-{suffix:}-new.csv'.format(**fmt)
        fmt['csvfn'] = csvfn
        df.to_csv(csvfn)
        print "{} new files to process".format(len(df))
        
        cmd = """isoclassify batch {mode:} {csvfn:} -o {outdir:} --plot none | grep mkdir  > isoclassify-{source:}-{suffix:}.tot """.format(**fmt)
        print cmd
        os.system(cmd)

