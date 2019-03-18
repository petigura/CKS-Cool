import os
import cPickle as pickle

import pandas as pd
import numpy as np
import cksspec.io
import ckscool.cuts
import cpsutils.io
import ckscool.calc
from astropy.io import ascii

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

    if table=='coldefs':
        tablefn = os.path.join(DATADIR,'column-definitions.txt')
        colspecs = [(0,1),(3,4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[20,100],
            names=['column','description']
        )

    elif table=='mathur17':
        df = cksspec.io.load_table('stellar17')
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)

    elif table=='johnson17':
        df = pd.read_csv('data/cks_physical_merged.csv',index_col=0)

    # Mathur 2017
    elif table=='stellar17':
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

    elif table=='m17':
        df = load_table('stellar17')
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

    # Merged tables
    elif table=='m17+gaia2':
        print "performing crossmatch on gaia2"

        df = pd.read_hdf(CKSGAIA_CACHEFN,'m17')
        df = df.rename(columns={'m17_kepmag':'kic_kepmag'})
        gaia = load_table('gaia2')
        stars = df['id_kic kic_kepmag'.split()].drop_duplicates()
        mbest,mfull = xmatch_gaia2(stars,gaia,'id_kic','gaia2')
        df = pd.merge(df,mbest.drop(['kic_kepmag'],axis=1),on='id_kic')

    elif table=='j17+m17':
        df = load_table('johnson17')
        m17 = load_table('mathur17')
        df = pd.merge(df,m17,on='id_kic')

    elif table=='cks-physical-merged':
        df = pd.read_csv('data/cks_physical_merged.csv',index_col=0)

    elif table=='cks-physical-merged+mathur17':
        df = load_table('cks-physical-merged')
        m17 = load_table('mathur17')
        df = pd.merge(df,m17, on='id_kic')

    elif table=='ckscool-targets-cuts':
        #table = 'koi-mullally15'
        #table = 'koi-thompson18'
        table = 'koi-thompson18-dr25'
        sample = 'koi-thompson18'
        cuttypes = [
            'none','notreliable','badteffphot','faint','giant','badparallax',
            'diluted'
        ]
        df = load_table(table)
        df = ckscool.cuts.add_cuts(df,cuttypes,sample)
        df.cuttypes = cuttypes

    elif table=='reamatch':
        fn = os.path.join(DATADIR, 'reamatch.csv')
        df = pd.read_csv(fn,index_col=None, usecols=range(6))
        namemap = {
            'obs':'id_obs',
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
        
    # Stellar sample
    elif table=='ckscool-stars':
        df = load_table('ckscool+smemp')

        # Add in ReaMatch parameters
        rm = load_table('reamatch')
        df = pd.merge(df, rm, how='left', on='id_obs')

        # Add in Furlan parameters
        f17 = load_table('fur17')
        df = pd.merge(df, f17, how='left',on=['id_kic','id_koi'])

        # Add in Kraus parameters
        k16 = load_table('kraus16')
        df = pd.merge(df, k16, how='left', on='id_koi')

        # Add in isoclassify parameters
        df['id_starname'] = df.id_name
        fn = os.path.join(DATADIR,'isoclassify_gaia2.csv')
        iso = pd.read_csv(fn)
        df = pd.merge(df, iso, on='id_starname', how='left')
        
    elif table=='ckscool-stars-cuts':
        df = load_table('ckscool-stars')
        #cuttypes = ['none','sb2','badspecparallax','dilutedao']
        cuttypes = ['none','sb2','dilutedao','badspecparallax',]
        table='temp'
        df = ckscool.cuts.add_cuts(df, cuttypes, table)
        df.cuttypes = cuttypes

    elif table == "ckscool-planets":
        df = load_table('ckscool-stars-cuts')
        df = df[~df.isany]
        df = ckscool.calc.update_planet_parameters(df)
        fpp = ckscool.io.load_table('fpp')
        df = pd.merge(df, fpp)

    elif table=='ckscool-planets-cuts':
        df = load_table('ckscool-planets')
        cuttypes = ['none','badprad','badpradprecision','badimpact','largefpp',]
        table = 'temp'
        df = ckscool.cuts.add_cuts(df, cuttypes, table)
        df.cuttypes = cuttypes
        
    elif table=='nrm-previous':
        # File is from an email that Adam sent me in 2017-11-02
        fn = os.path.join(DATADIR,'kraus/KOIours.txt')
        df = pd.read_table(fn,header=None,names=['name'])
        df['id_koi'] = df.name.str.slice(start=4).astype(int)
        df = df[['id_koi']]

    # ###########################
    # Tables from other papers #
    ############################
    elif table=='fpp':
        df = pd.read_csv('data/q1_q17_dr25_koifpp.csv',comment='#')
        namemap = {'kepoi_name':'id_koicand','fpp_prob':'fpp_prob'}
        df = df.rename(columns=namemap)[namemap.values()]


    # KOI tables
    elif table=='koi-thompson18':
        df = load_table_koi(table)

    elif table=='koi-thompson18-dr25':
        df = load_table_koi('koi-thompson18')
        dr25 = pd.read_hdf('../Kepler-Radius-Ratio/data/kepler_project_chains.hdf','dr25',)
        namemap = {
            'dr25_RD1_cum':'koi_ror',
            'dr25_RD1_cum_err1':'koi_ror_err1',
            'dr25_RD1_cum_err2':'koi_ror_err2'
        }
        cols = ['id_koicand'] +  namemap.values()
        df25 = dr25.rename(columns=namemap)[cols]
        print "swapping in dr25 radius ratios"
        print "       old"  
        print df.head(3)['id_koicand koi_ror koi_ror_err1 koi_ror_err2'.split()]
        df = df.drop(columns=namemap.values())
        df = pd.merge(df,df25)
        print "       new"
        print df.head(3)['id_koicand koi_ror koi_ror_err1 koi_ror_err2'.split()]


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
    star = load_table('m17+gaia2',cache=1)
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

def order_columns(df, verbose=False, drop=True):
    columns = list(df.columns)
    coldefs = load_table('coldefs',cache=0)
    cols = []
    for col in coldefs.column:
        if columns.count(col) == 1:
            cols.append(col)

    df = df[cols]
    if verbose and (len(cols) < len(columns)):
        print "table contains columns not defined in coldef"

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
    m.index = m.id_kic
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
