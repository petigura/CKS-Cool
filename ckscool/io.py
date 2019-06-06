import os
import cPickle as pickle

import pandas as pd
import numpy as np
from astropy.io import ascii
from scipy.io import idl

import cpsutils.io
import cksspec.io

import ckscool.cuts.occur
import ckscool.calc
import ckscool.pdplus
import ckscool.cuts.ckscool
import ckscool.comp
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
long_cadence_day = 29.7 / 60.0 / 24.0 # long cadence measurement in days


DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
CKSGAIA_CACHEFN = os.path.join(DATADIR,'cksgaia_cache.hdf')

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
            df = pd.read_hdf(cachefn,table, mode='a')
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

    if table=='coldefs':
        tablefn = os.path.join(DATADIR,'column-definitions.txt')
        colspecs = [(0,1),(3,4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[20,100],
            names=['column','description']
        )


    # Mathur 2017
    elif table=='m17':
        tablefn = os.path.join(DATADIR, 'kepler_stellar17.csv.gz')
        df = pd.read_csv(tablefn,sep='|',dtype={'st_quarters':'str'})
        namemap = {
            'kepid':'id_kic','kepmag':'kic_kepmag', 'teff': 'kic_steff',
            'st_quarters':'st_quarters','mass':'kic_smass',
            'st_radius':'kic_srad', 'jmag':'kic_jmag',
            'jmag_err':'kic_jmag_err','hmag':'kic_hmag',
            'hmag_err':'kic_hmag_err','kmag':'kic_kmag',
            'kmag_err':'kic_kmag_err',
            'degree_ra':'kic_ra', 'degree_dec':'kic_dec'
        }
        df = df.rename(columns=namemap)[namemap.values()]
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)

    # Gaia DR2
    elif table=='gaia2':
        fn = os.path.join('../CKS-Gaia/data/xmatch_m17_gaiadr2-result.csv')
        df = read_xmatch_gaia2(fn)
        # Systematic offset from Zinn et al. (2018)
        df['gaia2_sparallax'] += 0.053 

    # CDPP table 
    elif table=='cdpp':
        fn = os.path.join(DATADIR,'kic_q0_q17.dat')
        df = idl.readsav(fn) 
        df = df['kic']
        df = ckscool.pdplus.LittleEndian(df) # Deals with the byte order
        df = pd.DataFrame(df)
        df = df.rename(columns={
            'KEPMAG':'kepmag','KICID':'id_kic','CDPP3':'cdpp3','CDPP6':'cdpp6',
            'CDPP12':'cdpp12'}
        )
        df = df['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        for col in 'cdpp3 cdpp6 cdpp12'.split():
            cdpp = np.vstack(df.ix[:,col])
            cdpp = cdpp[:,1:18] # cut out q 0 which was't used
            observed = (cdpp > 0.0).astype(float) # True if observed 
            cdpp[cdpp==0.0] = np.nan
            
            cdppmed = np.nanmedian(cdpp,axis=1)
            df[col] = cdppmed
            df['log'+col] = np.log10(cdppmed)

        # calculate days that target was observed
        lc = np.array([ckscool.io.lc_per_quarter[i] for i in range(1,18)])
        lc = lc.reshape(1,17) 
        days = lc * long_cadence_day
        df['tobs'] = (observed * days).sum(axis=1)
        
    elif table=='m17+ber18+gaia2+cdpp':
        # Needed fro kmag
        #m17 = pd.read_hdf(CKSGAIA_CACHEFN,'m17')
        #import pdb;pdb.set_trace()
        #m17 = m17.rename(columns={'m17_kepmag':'kic_kepmag'})
        #m17 = load_table('ber18+gaia2',cache=1)
        m17 = load_table('m17',cache=1)
        cdpp = load_table('cdpp',cache=1)        
        ber18 = load_table('berger18',cache=1)
        gaia = load_table('gaia2',cache=1)
        df = pd.merge(m17,ber18,on=['id_kic'])
        df = pd.merge(df,gaia,on=['id_gaia2','id_kic'])
        df = pd.merge(df,cdpp)


    elif table=='field-cuts':
        df = load_table('m17+ber18+gaia2+cdpp')
        cuttypes = ['faint','giant','rizzuto',]
        star = ckscool.cuts.occur.add_cuts(df, cuttypes, 'field')
        starc = star.query('isany==False')
        df = starc

        #cuttypes = ckscool.cuts.field_cuttypes
        #df = ckscool.cuts.add_cuts(df, cuttypes, 'field')


    elif table=='planets-cuts1':
        star = load_table('m17+ber18+gaia2+cdpp')
        plnt = load_table('koi-thompson18-dr25')
        df = pd.merge(star,plnt)

        cuttypes = ['faint','giant','rizzuto','notreliable','lowsnr']
        df.sample = 'koi-thompson18'
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'koi-thompson18')

    # Stellar sample
    elif table=='planets-cuts1+iso':
        df = ckscool.io.load_iso_batch_table()

        # Add in isoclassify parameters
        df['id_starname'] = df.id_koi.apply(lambda x : "K{:05d}".format(x))
        fn = os.path.join(DATADIR,'isoclassify_gaia2.csv')
        iso = pd.read_csv(fn)
        df = pd.merge(df, iso, on='id_starname', how='left')

        # Add in ReaMatch parameters
        rm = load_table('reamatch')
        df = pd.merge(df, rm, how='left', on='id_koi')

        # Add in Furlan parameters
        f17 = load_table('fur17')
        df = pd.merge(df, f17, how='left',on=['id_kic','id_koi'])

        # Add in Kraus parameters
        k16 = load_table('kraus16')
        df = pd.merge(df, k16, how='left', on='id_koi')

        # Add in false positive prob
        fpp = ckscool.io.load_table('fpp')
        df = pd.merge(df, fpp)
        
        # Add in expected transit duration
        df['giso_tau0'] = 2.036 * df.koi_period**(1/3.) * df.giso_srho**(-1/3.0)
        df = ckscool.calc.update_planet_parameters(df)
        #df = order_columns(df, verbose=True)

        print "number of plabets {}".format(len(df))
        df = df.dropna(subset=['cks_steff'])
        print "number of plabets after after removing missking CKS {}".format(len(df))
        
    elif table=='planets-cuts2+iso':
        df = load_table('planets-cuts1+iso',cache=1)
        cuttypes = ['none','sb2','badprad','badimpacttau']
        table = 'temp'
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, table)
        df.cuttypes = cuttypes

    elif table=='ckscool-targets-cuts':
        # table = 'koi-mullally15'
        # table = 'koi-thompson18'
        table = 'koi-thompson18-dr25'
        sample = 'koi-thompson18'
        cuttypes = [
            'none','notreliable','badteffphot','faint','giant','badparallax',
            'diluted'
        ]
        df = load_table(table)
        df = ckscool.cuts.ckscool.add_cuts(df,cuttypes,sample)
        df.cuttypes = cuttypes

    elif table=='reamatch':
        fn = os.path.join(DATADIR, 'reamatch.csv')
        df = pd.read_csv(fn,index_col=None, usecols=range(6))
        df['id_koi'] = df.name.str.slice(start=-5).astype(int)
        namemap = {
            'id_koi':'id_koi',
            'is_sb2':'rm_sb2',
        }
        df = df.rename(columns=namemap)[namemap.values()]

    # Spectroscopic parameters
    elif table=='smemp':
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

    elif table=='kbc':
        df = cpsutils.kbc.loadkbc()
        b = ( df.type.str.contains('t') 
              & df.name.str.contains('^K\d{5}$|^CK\d{5}$') )
        df = df[b]
        df['id_koi'] = df.name.str.slice(start=-5).astype(int)
        df = df['obs name id_koi'.split()]
        namemap = {'obs':'id_obs','name':'id_name'}
        df = df.rename(columns=namemap)

    elif table=='ckscool+smemp':
        df = ckscool.io.load_table('ckscool-targets-cuts')
        df = df[~df.isany]

        # Grab observation using kbc
        kbc = load_table('kbc')
        kbc = kbc.groupby('id_koi', as_index=False).nth(-1)
        df = pd.merge(df, kbc)

        # Add in specmatch parameters
        sm = load_table('smemp')
        df = pd.merge(df, sm, on=['id_obs','id_name'], how='left')


    elif table=='nrm-previous':
        # File is from an email that Adam sent me in 2017-11-02
        fn = os.path.join(DATADIR,'kraus/KOIours.txt')
        df = pd.read_table(fn,header=None,names=['name'])
        df['id_koi'] = df.name.str.slice(start=4).astype(int)
        df = df[['id_koi']]

    ############################
    # Tables from other papers #
    ############################
    elif table=='fpp':
        df = pd.read_csv('data/q1_q17_dr25_koifpp.csv',comment='#')
        namemap = {'kepoi_name':'id_koicand','fpp_prob':'fpp_prob'}
        df = df.rename(columns=namemap)[namemap.values()]

    # KOI tables
    elif table=='koi-thompson18':
        df = load_table_koi(table)

    # KOI tables
    elif table=='koi-thompson18-dr25-chains':
        df = load_table_koi(table)
        t18.index = t18.id_koicand
        df = pd.read_csv('data/dr25-mcmc-chains_missing.txt',sep=' ',names=['s','x'])
        df['id_koicand'] = df.s.apply(lambda x : x.split('-')[-1]).str.upper()
        t18 = t18.drop(df.id_koicand)
        df = t18
        
    elif table=='koi-thompson18-dr25':
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

    elif table=='koi-coughlin16':
        df = load_table_koi(table)

    elif table=='koi-mullally15':
        df = load_table_koi(table)

    # Furlan 2017
    elif table=='furlan17-tab2':
        fn = os.path.join(DATADIR,'furlan17/Table2.txt')
        df = read_furlan17_table2(fn)

    elif table=='furlan17-tab9':
        fn = os.path.join(DATADIR,'furlan17/Table9.txt')
        df = read_furlan17_table9(fn)

    elif table=='fur17':
        tab2 = load_table('furlan17-tab2')
        tab9 = load_table('furlan17-tab9')
        cols = 'id_koi f17_rcf_avg f17_rcf_avg_err'.split()
        df = pd.merge(tab2,tab9[cols],how='left',on='id_koi')

    elif table=='kraus16':
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
    elif table=='mann13':
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
        cks = load_table('ckscool-stars-cuts')
        m13 = load_table('mann13')
        df = pd.merge(cks,m13,on='id_koi')

    # Dressing et al (2013)
    elif table=='dressing13':
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

    elif table=='ckscool-dressing13':
        cks = load_table('ckscool-stars-cuts')
        d13 = load_table('dressing13')
        df = pd.merge(cks,d13,on='id_kic')


    # Berger al. (2018)
    elif table=='berger18':
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

    # Brewer et al. (2018)
    elif table=='brewer18':
        tab = ascii.read('data/brewer18/apjsaad501t3_mrt.txt')
        df = tab.to_pandas()
        df = df[df.Name.str.contains('KOI-.*\d$')]
        df['id_koi']  = df.Name.apply(lambda x : x.split('-')[1]).astype(int)
        df['steff'] = df['Teff']
        df = df['id_koi steff'.split()]
        df = add_prefix(df,'b18_')

    elif table=='ckscool-brewer18':
        cks = load_table('ckscool-stars-cuts')
        m13 = load_table('brewer18')
        df = pd.merge(cks,m13,on='id_koi')

    # CKS-VII
    elif table=='cksgaia-planets-filtered':
        df = pd.read_hdf(CKSGAIA_CACHEFN,'cksgaia-planets-filtered')
    elif table=='cksgaia-planets':
        df = pd.read_hdf(CKSGAIA_CACHEFN,'cksgaia-planets')

    else:
        assert False, "table {} not valid table name".format(table)
    return df

def read_furlan17_table2(fn):
    df = pd.read_csv(fn,sep='\s+')
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

    df = pd.read_csv(fn,sep='\s+',skiprows=2,names=names)
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

def load_iso_batch_table(mode='isoclassify'):
    """
    This function defines the values that are sent to isoclassify .
    """
    plnt = ckscool.io.load_table('planets-cuts1')
    plntc = plnt.query('isany==False')
    star = plntc.groupby('id_koi', as_index=False).nth(0)
    star = star.copy()
    star['cks_steff'] = np.nan
    star['cks_steff_err'] = np.nan
    star['cks_smet'] = np.nan
    star['cks_smet_err'] = np.nan
    star['cks_svsini'] = np.nan
    star['cks_svsini_err'] = np.nan
    star['cks_sprov'] = None
    star.index = star.id_koi

    # Load SpecMatch-Emp and set uncertainties
    kbc = ckscool.io.load_table('kbc')
    kbc = kbc.groupby('id_koi', as_index=False).nth(-1)
    df = pd.read_csv('data/specmatch-emp_results.csv')
    df = df.dropna(subset=['name'])
    namemap = {
       'obs':'id_obs',
       'name':'id_name',
       'teff':'cks_steff',
       'teff_err':'cks_steff_err',
       'fe':'cks_smet',
       'fe_err':'cks_smet_err',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df['cks_sprov'] = 'emp'
    df['cks_steff_err'] = 60
    df['cks_smet_err'] = 0.12
    sme = pd.merge(kbc, df, on=['id_obs','id_name'], how='left')
    sme = pd.merge(star[['id_koi']],sme) # only retain values in planet list
    sme.index = sme.id_koi


    # Load up CKS-I and set uncertainties
    cks1 = pd.read_csv('data/cks_physical_merged.csv')
    cks1['id_koi'] = cks1.id_koicand.str.slice(start=1,stop=6).astype(int)
    cks1 = cks1.groupby('id_koi',as_index=False).nth(0)
    cks1['cks_steff_err'] = 100
    cks1['cks_smet_err'] = 0.06
    cks1['cks_sprov'] = 'cks1'
    cks1 = pd.merge(star[['id_koi']],cks1) # only retain values in planet list
    cks1.index=cks1.id_koi

    # Load up SpecMatch-Syn and set uncertainties 
    df = pd.read_csv('data/specmatch-syn_results.csv')
    df = df.dropna(subset=['name'])
    namemap = {
       'obs':'id_obs',
       'name':'id_name',
       'teff':'cks_steff',
       'teff_err':'cks_steff_err',
       'fe':'cks_smet',
       'fe_err':'cks_smet_err',
       'vsini':'cks_svsini',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df['cks_sprov'] = 'syn'
    df['cks_steff_err'] = 100
    sms = pd.merge(kbc, df, on=['id_obs','id_name'], how='left')
    cks1 = pd.merge(star[['id_koi']],cks1) # only retain values in planet list
    sms.index=sms.id_koi

    # Cut off for using emp
    idxsmemp = sme.query('cks_steff < 4700').index
    idxsmsyn = sme.query('cks_steff > 4700').index


    star.fillna(value=sme.loc[idxsmemp],inplace=True) # Cool plnts fill in with emp
    star.fillna(value=cks1,inplace=True) # Stars that are in CKS-I fill, also take vsini for everythinh
    star.fillna(value=sms.loc[idxsmsyn],inplace=True) #
    star['cks_sprov'].fillna(value='None',inplace=True)
    star['ber18_srad'] = star.ber18_srad.astype(float)
    star = star.sort_values(by='id_koi')

    return plntc



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
    df = add_prefix(df, 'gaia2_')
    return df


def xmatch_gaia2(df, gaia, key, gaiadr):
    """
    Crossmatch the sources in Gaia 2

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR2 table
        key (str): key to join on
        gaiadr (str): {'gaiadr1','gaiadr2'}
    """

    id_gaia = 'id_{}'.format(gaiadr)
    assert len(df.id_kic)==len(df.id_kic.drop_duplicates()), "No duplicate stars"
    gaia[id_gaia]= gaia[id_gaia].astype(str)  # handle nans

    # just want the stars
    m = pd.merge(df,gaia,on=key,how='left')
    m[id_gaia]= m[id_gaia].fillna('-99')
    m[id_gaia]= m[id_gaia].astype(np.int64)
    ndf = len(df)
    print "max(gaia_angdist) = {} (arcsec)".format(m[gaiadr+'_angdist'].max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    #m.index = m.id_kic
    g = m.groupby('id_kic')
    m[gaiadr+'_gflux_sum'] = g[gaiadr+'_gflux'].sum()
    m['absdiff_gmag_kepmag'] = np.abs(m['gaia2_gmag'] - m['kic_kepmag'])
    m['gaia2_n_8arcsec'] = g.size()

    # Match candidate
    mbest = m.query('gaia2_angdist < 1 and abs(kic_kepmag - gaia2_gmag) < 0.5') # within 1 arcsec
    mbest['gaia2_n_1arcsec'] = mbest.groupby('id_kic').size()

    print "{} gaia sources within 1 arcsec of {} target sources".format(
        len(mbest),mbest.id_kic.drop_duplicates().count()
    )

    mbest = mbest.sort_values(by=['id_kic','absdiff_gmag_kepmag'])
    g = mbest.groupby('id_kic',as_index=False)
    mbest['gaia2_n_1arcsec'] = g.size()
    mbest = g.nth(0) 
    mbest['gaia2_gflux_ratio'] = mbest.eval('gaia2_gflux_sum / gaia2_gflux')
    return mbest, m 


