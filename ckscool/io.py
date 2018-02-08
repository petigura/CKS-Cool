import os
import cPickle as pickle

import pandas as pd
import numpy as np
import cksspec.io

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')

def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksmet

    Args:
        table (str): name of table. must be one of
            - nea 


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table)
            print "read table {} from {}".format(table,cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table)
        return df

    elif table=='koi-thompson18':
        df = load_table_koi(table)

    elif table=='koi-coughlin16':
        df = load_table_koi(table)

    elif table=='koi-mullally15':
        df = load_table_koi(table)

    elif table=='mathur17':
        df = cksspec.io.load_table('stellar17')
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)

    elif table=='cks-physical-merged':
        df = pd.read_csv('data/cks_physical_merged.csv',index_col=0)

    elif table=='cks-physical-merged+mathur17':
        df = load_table('cks-physical-merged')
        m17 = load_table('mathur17')
        df = pd.merge(df,m17, on='id_kic')
    else:
        assert False, "table {} not valid table name".format(table)
    return df

def load_table_koi(table):
    """
    Helper function to load KOI tables
    """
    if table=='koi-thompson18':
        csvfn = os.path.join(DATADIR,'q1_q17_dr25_koi.csv')
        df = pd.read_csv(csvfn,comment='#',index_col=0)

    elif table=='koi-coughlin16':
        csvfn = os.path.join(DATADIR,'q1_q17_dr24_koi.csv')
        df = pd.read_csv(csvfn,comment='#',index_col=0)

    elif table=='koi-mullally15':
        csvfn = os.path.join(DATADIR,'q1_q16_koi.csv')
        df = pd.read_csv(csvfn, comment='#', skipinitialspace=True)

    namemap = {
        'kepid':'id_kic',
        'kepoi_name':'id_koicand',
        'kepler_name':'id_kepler_name',
    }
    
    df = df.rename(columns=namemap)
    df['id_koi'] = df.id_koicand.str.slice(start=1,stop=6).astype(int)
    star = load_table('mathur17',cache=1)
    df = pd.merge(df,star)
    return df 

def add_prefix(df,prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = prefix + col 
    df = df.rename(columns=namemap)
    return df

def sub_prefix(df, prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = col.replace(prefix,'') 
    df = df.rename(columns=namemap)
    return df

