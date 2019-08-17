"""
Module for CKS-Cool I/O


"""
import os
import cPickle as pickle
import warnings

import pandas as pd
import numpy as np
from astropy.io import ascii
from scipy.io import idl
import tables
from astropy.table import Table

import ckscool.cuts.occur
import ckscool.calc
import ckscool.comp
from ckscool.pdplus import LittleEndian
import ckscool.gaia
import ckscool.occur
import ckscool.plot.occur
import ckscool._isoclassify

# Ignore the Natural name warning
warnings.simplefilter('ignore', tables.NaturalNameWarning)
warnings.simplefilter('ignore', pd.errors.PerformanceWarning)

# Define paths to various cache files. DATADIR stores data tables that
# should not change with different code runs. CACHEDIR stores
# temporary files that should be cleared out when the analysis
# changes. The CACHEDIR can be set for different git branches
FILE = os.path.dirname(__file__)
DATADIR = os.path.join(FILE, '../data/')
CACHEDIR = os.path.join(FILE, '../cache/field-m17-plnt-giso/')
CACHEFN = os.path.join(CACHEDIR, 'load_table_cache.hdf')
os.system('mkdir -p {}'.format(CACHEDIR)) # creates CACHEDIR if doesn't exist
KBCFN = os.path.join(DATADIR,'kbcvel.csv')

def load_table(table, cache=0, verbose=False, cachefn=None):
    """Load tables used in ckscool

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
    if cachefn is None:
        cachefn = CACHEFN

    if cache == 1:
        try:
            df = pd.read_hdf(cachefn, table, mode='a')
            print "read table {} from {}".format(table, cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache = 2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache = 2

    if cache == 2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn, table)
        return df

    if table == 'coldefs':
        tablefn = os.path.join(DATADIR, 'column-definitions.txt')
        colspecs = [(0, 1), (3, 4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[20, 100],
            names=['column', 'description']
        )

    # Mathur 2017
    elif table == 'm17':
        tablefn = os.path.join(DATADIR, 'kepler_stellar17.csv.gz')
        df = pd.read_csv(tablefn, sep='|', dtype={'st_quarters':'str'})
        namemap = {
            'kepid':'id_kic', 'kepmag':'m17_kepmag', 'teff': 'm17_steff',
            'st_quarters':'st_quarters', 'mass':'m17_smass',
            'st_radius':'m17_srad', 'jmag':'m17_jmag',
            'jmag_err':'m17_jmag_err', 'hmag':'m17_hmag',
            'hmag_err':'m17_hmag_err', 'kmag':'m17_kmag',
            'kmag_err':'m17_kmag_err',
            'degree_ra':'m17_ra', 'degree_dec':'m17_dec'
        }
        df = df.rename(columns=namemap)[namemap.values()]

    # Gaia DR2
    elif table=='gaia2':
        #fn = os.path.join(DATADIR,'xmatch_m17_gaiadr2-result.csv')
        #fn = os.path.join(DATADIR,'xmatch_gaia2_m17_ruwe_tmass-result.vot.gz')
        fn = os.path.join(DATADIR,'xmatch_gaia2_m17_ruwe_tmass-result.csv')
        df = pd.read_csv(fn)
        # Systematic offset from Zinn et al. (2018)
        df['gaia2_sparallax'] += 0.053

    # CDPP table
    elif table == 'cdpp':
        # Pulled from Kepler data characteristics
        lc_per_quarter = {
            0:476,
            1:1639,
            2:4354,
            3:4370,
            4:4397,
            5:4634,
            6:4398,
            7:4375,
            8:3279,
            9:4768,
            10:4573,
            11:4754,
            12:4044,
            13:4421,
            14:4757,
            15:4780,
            16:4203,
            17:1556,
        }
        fn = os.path.join(DATADIR, 'kic_q0_q17.dat')
        df = idl.readsav(fn)
        df = df['kic']
        df = LittleEndian(df) # Deals with the byte order
        df = pd.DataFrame(df)
        df = df.rename(columns={
            'KEPMAG':'kepmag', 'KICID':'id_kic',
            'CDPP3':'cdpp3', 'CDPP6':'cdpp6', 'CDPP12':'cdpp12'}
        )
        df = df['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        for col in 'cdpp3 cdpp6 cdpp12'.split():
            cdpp = np.vstack(df.ix[:, col])
            cdpp = cdpp[:, 1:18] # cut out q 0 which was't used
            observed = (cdpp > 0.0).astype(float) # True if observed 
            cdpp[cdpp == 0.0] = np.nan
            cdppmed = np.nanmedian(cdpp, axis=1)
            df[col] = cdppmed
            df['log'+col] = np.log10(cdppmed)

        # calculate days that target was observed
        lc = np.array([lc_per_quarter[i] for i in range(1, 18)])
        lc = lc.reshape(1, 17) 
        long_cadence_day = 29.7 / 60.0 / 24.0 # long cadence length in days
        days = lc * long_cadence_day
        df['tobs'] = (observed * days).sum(axis=1)

    elif table == 'm17+cdpp+gaia2+ber19':
        m17 = load_table('m17', cache=1)
        # Store to separate cache to prevent long reload times
        cachefn = os.path.join(DATADIR,'cdpp.hdf') 
        cdpp = load_table('cdpp',cache=1, cachefn=cachefn)        
        ber = load_table('berger19', cache=1)
        gaia = load_table('gaia2', cache=1)
        df = pd.merge(m17,cdpp)
        df = pd.merge(df,gaia)
        df = pd.merge(df,ber,on=['id_kic'])

    elif table == 'field-cuts':
        df = load_table('m17+cdpp+gaia2+ber19',cache=1)
        cuttypes = ['none','faint','giant','rizzuto']
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'field')

    elif table == 'planets-cuts1':
        star = load_table('m17+cdpp+gaia2+ber19',cache=1)
        plnt = load_table('koi-thompson18-dr25')
        df = pd.merge(star,plnt)
        cuttypes = ['none','faint','giant','rizzuto','notreliable','lowsnr']
        df.sample = 'koi-thompson18'
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'koi-thompson18')

    # Results from isoclassify table
    elif table == 'iso':
        source = 'cks1'
        star0 = ckscool._isoclassify.load_stellar_parameters(source)
        iso = pd.read_csv('data/isoclassify_{}.csv'.format(source),index_col=0)
        iso['id_koi'] = iso.id_starname.str.slice(start=-5).astype(int)
        cks1 = pd.merge(star0,iso)
        cks1.index = cks1.id_koi
        cks1['cks_sprov'] = 'cks1' # note that obs didn't transfer over here.

        source = 'smemp'
        star0 = ckscool._isoclassify.load_stellar_parameters(source)
        iso = pd.read_csv('data/isoclassify_{}.csv'.format(source),index_col=0)
        iso['id_koi'] = iso.id_starname.str.slice(start=-5).astype(int)
        smemp = pd.merge(star0,iso)
        smemp.index = smemp.id_koi

        source = 'smsyn'
        star0 = ckscool._isoclassify.load_stellar_parameters(source)
        iso = pd.read_csv('data/isoclassify_{}.csv'.format(source),index_col=0)
        iso['id_koi'] = iso.id_starname.str.slice(start=-5).astype(int)
        smsyn = pd.merge(star0,iso)
        smsyn.index = smsyn.id_koi
        
        df = []
        smemplimit = 4800
        for id_koi in smemp.id_koi.drop_duplicates():
            try:
                smemp_steff = smemp.loc[id_koi,'cks_steff']
                if smemp_steff < smemplimit:
                    df.append(smemp.loc[id_koi])
                    continue
            except KeyError:
                continue

            try:
                df.append(cks1.loc[id_koi])
                continue
            except KeyError:
                pass

            try:
                df.append(smsyn.loc[id_koi])
                continue
            except KeyError:
                pass

        df = pd.DataFrame(df)
        #df['cks_svsini'] = df.cks_svsini.fillna(smsyn.cks_svsini)
        df = df.reset_index()

    # Stellar sample
    elif table=='planets+iso':
        #star = load_table('m17+cdpp+gaia2+ber18')
        plnt = load_table('koi-thompson18-dr25')
        iso = load_table('iso',cache=1)
        df = pd.merge(plnt, iso, how='left',on=['id_koi','id_kic'])

        # Add in ReaMatch parameters. If the provenence of the
        # parameter is CKS-I then in means that the star passed the
        # reamatch pipeline
        rm = load_table('reamatch')
        df = pd.merge(df, rm, how='left', on='id_koi')
        idx = df.query('cks_sprov == "cks1"').index
        df.loc[idx,'rm_sb2'] = 1

        fn = os.path.join(DATADIR,'kolbl15/table9.tex')
        k15 = pd.read_table(fn,sep='&',skiprows=9,header=None,nrows=64)
        k15 = k15[[0]]
        k15 = k15.rename(columns={0:'id_koi'})
        k15['id_koi'] = k15.id_koi.str.replace('\t','').str.strip().replace('',None).astype(int)
        k15 = k15.drop_duplicates()
        temp = pd.merge(k15,df[['id_koi']]).drop_duplicates()
        df = df.set_index('id_koi')
        for id_koi in temp.id_koi:
            df.loc[id_koi,'rm_sb2'] = 5

        # Add in Furlan parameters
        df = df.reset_index()
        f17 = load_table('fur17')
        df = pd.merge(df, f17, how='left',on=['id_kic','id_koi'])

        # Add in Kraus parameters
        k16 = load_table('kraus16')
        df = pd.merge(df, k16, how='left', on='id_koi')

        # Add in false positive prob
        fpp = load_table('fpp')
        df = pd.merge(df, fpp)
        
        # Add in expected transit duration
        df['giso_tau0'] = 2.036 * df.koi_period**(1/3.) * df.giso_srho**(-1/3.0)
        
        print "updating parameters"
        df = ckscool.calc.update_planet_parameters(df)
        #df = order_columns(df, verbose=True)

        print "number of planets {}".format(len(df))
        df = df.dropna(subset=['cks_steff'])
        print "number of planets after after removing missing CKS {}".format(len(df))
        
    elif table == 'planets-cuts1+iso':
        df = load_table('planets+iso',cache=1)
        planets1 = load_table('planets-cuts1',cache=2)
        df = ckscool.cuts.occur.add_cuts(df, planets1.cuttypes, 'koi-thompson18')
        df = df[~df.isany]

    elif table == 'planets-cuts2+iso':
        df = load_table('planets-cuts1+iso',cache=2)
        cuttypes = ['none','badvsini','sb2','badspecparallax','badprad','badpradprec','badimpacttau']
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'koi-thompson18')

    elif table == 'reamatch':
        fn = os.path.join(DATADIR, 'reamatch.csv')
        df = pd.read_csv(fn,index_col=None, usecols=range(6))
        df['id_koi'] = df.name.str.slice(start=-5).astype(int)
        namemap = {'id_koi':'id_koi','is_sb2':'rm_sb2'}
        df = df.rename(columns=namemap)[namemap.values()]

    # Spectroscopic parameters
    elif table == 'smemp':
        df = pd.read_csv('data/specmatch-emp_results.csv')
        df = df.dropna(subset=['name'])
        namemap = {
            'obs':'id_obs',
            'name':'id_name',
            'teff':'sm_steff',
            'teff_err':'sm_steff_err',
            'radius':'sm_srad',
            'radius_err':'sm_srad_err',
            'fe':'sm_smet',
            'fe_err':'sm_smet_err',
        }
        df = df.rename(columns=namemap)[namemap.values()]

    elif table == 'nrm-previous':
        # File is from an email that Adam sent me in 2017-11-02
        fn = os.path.join(DATADIR,'kraus/KOIours.txt')
        df = pd.read_table(fn,header=None,names=['name'])
        df['id_koi'] = df.name.str.slice(start=4).astype(int)
        df = df[['id_koi']]

    ############################
    # Tables from other papers #
    ############################
    elif table == 'fpp':
        df = pd.read_csv('data/q1_q17_dr25_koifpp.csv',comment='#')
        namemap = {'kepoi_name':'id_koicand','fpp_prob':'fpp_prob'}
        df = df.rename(columns=namemap)[namemap.values()]

    # KOI tables
    elif table=='koi-thompson18':
        df = load_table_koi(table)

    # KOI tables
    elif table == 'koi-thompson18-dr25-chains':
        df = load_table_koi(table)
        t18.index = t18.id_koicand
        df = pd.read_csv('data/dr25-mcmc-chains_missing.txt',sep=' ',names=['s','x'])
        df['id_koicand'] = df.s.apply(lambda x : x.split('-')[-1]).str.upper()
        t18 = t18.drop(df.id_koicand)
        df = t18
        
    elif table == 'koi-thompson18-dr25':
        df = load_table_koi('koi-thompson18')
        dr25 = pd.read_hdf('data/dr25-mcmc-summary.hdf','dr25',)
        namemap = {
            'dr25_RD1_cum':'dr25_ror',
            'dr25_RD1_cum_err1':'dr25_ror_err1',
            'dr25_RD1_cum_err2':'dr25_ror_err2',
            'dr25_BB1_cum':'dr25_impact',
            'dr25_BB1_cum_err1':'dr25_impact_err1',
            'dr25_BB1_cum_err2':'dr25_impact_err2',
            'dr25_TAU1_cum':'dr25_tau',
            'dr25_TAU1_cum_err1':'dr25_tau_err1',
            'dr25_TAU1_cum_err2':'dr25_tau_err2',
            'fgraz':'dr25_fgraz',
            'autocorr_over_length':'dr25_autocorr_over_length',
        }

        cols = ['id_koicand'] +  namemap.values()
        dr25 = dr25.rename(columns=namemap)[cols]
        df = pd.merge(df,dr25)

    elif table == 'koi-coughlin16':
        df = load_table_koi(table)

    elif table == 'koi-mullally15':
        df = load_table_koi(table)

    # Furlan 2017
    elif table == 'furlan17-tab2':
        fn = os.path.join(DATADIR,'furlan17/Table2.txt')
        df = read_furlan17_table2(fn)

    elif table=='furlan17-tab9':
        fn = os.path.join(DATADIR,'furlan17/Table9.txt')
        df = read_furlan17_table9(fn)

    elif table == 'fur17':
        tab2 = load_table('furlan17-tab2')
        tab9 = load_table('furlan17-tab9')
        cols = 'id_koi f17_rcf_avg f17_rcf_avg_err'.split()
        df = pd.merge(tab2,tab9[cols],how='left',on='id_koi')

    elif table == 'kraus16':
        fn = os.path.join(DATADIR,'kraus16/table7.dat')
        readme = os.path.join(DATADIR,'kraus16/ReadMe')
        df = ascii.read(fn,readme=readme)
        df = df.to_pandas()
        namemap = {
            "KOI":"id_koi",
            'M2/M1':'massratio',
            'Sep':'sep',
            'l_M2/M1':'massratio_ul'
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df['massratio_ul'] = (df.massratio_ul=="<").astype(int)
        df = df[df.massratio_ul.astype(bool)]
        g = df.sort_values(by=['massratio']).groupby('id_koi',as_index=False)
        df = g.nth(-1)
        df = df.rename(columns={'massratio':'max_massratio'})
        df = add_prefix(df,'k16_')

    # Mann et al. (2013)
    elif table == 'mann13':
        fn = os.path.join(DATADIR,'mann13/table3.dat') 
        readme = os.path.join(DATADIR,'mann13/ReadMe') 
        df = ascii.read(fn,readme=readme)
        df = df.to_pandas()
        namemap = {
            'KOI':'id_koi',
            'Teff':'steff',
            'R*':'srad',
            'e_R*':'srad_err',
            'M*':'smass',
            'e_M*':'smass_err',
            'L*':'slum',
            'e_L*':'slum_err',
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df['id_koi'] = df['id_koi'].astype(int)
        df['steff_err'] = 57
        df = add_prefix(df,'m13_')

    elif table=='ckscool-mann13':
        cks = load_table('planets-cuts2+iso').groupby('id_koi',as_index=False).nth(0)
        m13 = load_table('mann13')
        df = pd.merge(cks,m13,on='id_koi')

    # Dressing et al (2013)
    elif table == 'dressing13':
        fn = os.path.join(DATADIR,'dressing13/table1.dat') 
        readme = os.path.join(DATADIR,'dressing13/ReadMe') 
        df = ascii.read(fn,readme=readme)
        df = df.to_pandas()
        namemap = {
            'KIC':'id_kic',
            'Teff':'steff',
            'E_Teff':'steff_err1',
            'e_Teff':'steff_err2',
            'R*':'srad',
            'E_R*':'srad_err1',
            'e_R*':'srad_err2',
        }
        df = df.rename(columns=namemap)[namemap.values()]

        for c in df.columns:
            if c.count('err2')==1:
                df[c] *= -1 
        df = add_prefix(df,'d13_')

    elif table == 'ckscool-dressing13':
        cks = load_table('planets-cuts2+iso').groupby('id_koi',as_index=False).nth(0)
        d13 = load_table('dressing13')
        df = pd.merge(cks,d13,on='id_kic')

    # Berger al. (2018)
    elif table == 'berger18':
        tab = ascii.read('data/berger18/apjaada83t1_mrt.txt')
        df = tab.to_pandas()
        namemap = {
            'KIC':'id_kic',
            'Gaia':'id_gaia2',
            'Teff':'steff',
            'e_Teff':'steff_err',
            'R*':'srad',
            'e_R*':'srad_err',
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df = add_prefix(df,'ber18_')

    # Berger al. (2019)
    elif table == 'berger19':
        df = pd.read_csv('data/GKSPC_InOut_V1.csv')
        df = df.rename(columns={'KIC':'id_kic'})
        namemap = {
            'KIC':'id_kic',
            'iso_teff':'steff',
            'iso_teff_err1':'steff_err1',
            'iso_teff_err2':'steff_err2',
            'iso_rad':'srad',
            'iso_rad_err1':'srad_err1',
            'iso_rad_err2':'srad_err2',
            'iso_mass':'smass',
            'iso_mass_err1':'smass_err1',
            'iso_mass_err2':'smass_err2',
        }

        df = df.rename(columns=namemap)[namemap.values()]
        df = add_prefix(df,'ber19_')
    elif table == 'ckscool-brewer18':
        cks = load_table('planets-cuts2+iso').groupby('id_koi',as_index=False).nth(0)
        m13 = load_table('brewer18')
        df = pd.merge(cks,m13,on='id_koi')


    # Brewer et al. (2018)
    elif table == 'brewer18':
        tab = ascii.read('data/brewer18/apjsaad501t3_mrt.txt')
        df = tab.to_pandas()
        df = df[df.Name.str.contains(r'KOI-.*\d$')]
        df['id_koi']  = df.Name.apply(lambda x : x.split('-')[1]).astype(int)
        df['steff'] = df['Teff']
        df = df['id_koi steff'.split()]
        df = add_prefix(df,'b18_')

    else:
        assert False, "table {} not valid table name".format(table)
    return df

def read_furlan17_table2(fn):
    df = pd.read_csv(fn,sep=r'\s+')
    namemap = {'KOI':'id_koi','KICID':'id_kic','Observatories':'ao_obs'}
    df = df.rename(columns=namemap)[namemap.values()]
    df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] 
    df = add_prefix(df,'f17_')
    return df

def read_furlan17_table9(fn):
    names = """
    id_koi hst hst_err i i_err 692 692_err lp600 lp600_err jmag jmag_err 
    kmag kmag_err jkdwarf jkdwarf_err jkgiant jkgiant_err
    avg avg_err
    """.split()

    df = pd.read_csv(fn,sep=r'\s+',skiprows=2,names=names)
    df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] 
    df = add_prefix(df,'f17_rcf_')
    return df 

def load_table_koi(table):
    """
    Helper function to load KOI tables
    """
    if table=='koi-thompson18':
        csvfn = os.path.join(DATADIR,'q1_q17_dr25_koi.csv')
        df = pd.read_csv(csvfn,comment='#',index_col=0)
        names = [
            'koi_disposition',
            'koi_score',
            'koi_period',
            'koi_period_err1',
            'koi_period_err2',
            'koi_ror',
            'koi_ror_err1',
            'koi_ror_err2',
            'koi_prad',
            'koi_prad_err1',
            'koi_prad_err2',
            'koi_impact',
            'koi_impact_err1',
            'koi_impact_err2',
            'kepid',
            'kepoi_name',
            'kepler_name',
            'koi_max_sngle_ev',
            'koi_max_mult_ev',
            'koi_model_snr',
            'koi_count',
            'koi_num_transits',
            'koi_tce_plnt_num',
            'koi_tce_delivname',
            'koi_quarters',
            'koi_bin_oedp_sig',

        ]
        df = df[names]

    elif table=='koi-coughlin16':
        csvfn = os.path.join(DATADIR,'q1_q17_dr24_koi.csv')
        df = pd.read_csv(csvfn,comment='#',index_col=0)

    elif table=='koi-mullally15':
        csvfn = os.path.join(DATADIR,'q1_q16_koi.csv')
        df = pd.read_csv(csvfn, comment='#', skipinitialspace=True)
    else:
        assert False, "{} not valid table".format(table)

    namemap = {
        'kepid':'id_kic',
        'kepoi_name':'id_koicand',
        'kepler_name':'id_kepler_name',
    }
    
    df = df.rename(columns=namemap)
    df['id_koi'] = df.id_koicand.str.slice(start=1,stop=6).astype(int)
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

def order_columns(df, verbose=False, drop=False):
    df0 = df.copy()
    columns = list(df.columns)
    coldefs = load_table('coldefs',cache=0)
    cols = []
    for col in coldefs.column:
        if columns.count(col) == 1:
            cols.append(col)

    df = df[cols]
    if verbose and (len(cols) < len(columns)):
        mcols = list(df0.drop(cols,axis=1).columns)
        print "following columns are not defined"
        print mcols

    if not drop:
        df = pd.concat([df,df0[mcols]],axis=1)
    return df

def load_object(key,cache=0,verbose=1):
    pklfn = os.path.join(CACHEDIR,key+'.pkl')
    if cache == 1:
        try:
            with open(pklfn,'r') as f:
                obj = pickle.load(f)
                if verbose:
                    print "read {} from {}".format(obj,pklfn)
                return obj

        except IOError:
            print "Could not find cache file: %s" % pklfn
            print "Building cache..."
            cache = 2

    if cache == 2:
        obj = load_object(key, cache=0)
        print "writing {} to {}".format(obj,pklfn)
        with open(pklfn,'w') as f:
            pickle.dump(obj,f)
        return obj
        
    if key.count('cp') == 1:
        occurkey = key.replace('cp','occur')
        occ = load_object(occurkey, cache=1)
        obj = ckscool.plot.occur.load_contour_plotter(occ)

    elif key.count('occur') == 1:
        bits = key.split('_')
        limits = {}
        for bit in bits:
            if bit.count('smass'):
                smass1, smass2 = bit.replace('smass=','').split('-')
                limits['smass1'] = float(smass1)
                limits['smass2'] = float(smass2)

            if bit.count('bmr'):
                bmr1, bmr2 = bit.replace('bmr=','').split('-')
                limits['bmr1'] = float(bmr1)
                limits['bmr2'] = float(bmr2)

        obj = ckscool.occur.load_occur(limits)
        obj.comp.__delattr__('stars')

    return obj
