from __future__ import print_function

"""
Module for CKS-Cool I/O
"""
import os
import warnings
from collections import OrderedDict
import re
import cPickle as pickle

import pandas as pd
import numpy as np
from astropy.io import ascii
from astropy.time import Time
from scipy.io import idl
import tables
from astropy.table import Table
import astropy.io.ascii
import lmfit
import scipy

import ckscool.cuts.occur
import ckscool.calc
import ckscool.comp
from ckscool.pdplus import LittleEndian
import ckscool.gaia
import ckscool.occur
import ckscool._isoclassify
import ckscool.fit
import ckscool.fitdetected
import ckscool.gradient

# Ignore the Natural name warning
warnings.simplefilter('ignore', tables.NaturalNameWarning)
warnings.simplefilter('ignore', pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# Define paths to various cache files. DATADIR stores data tables that
# should not change with different code runs. CACHEDIR stores
# temporary files that should be cleared out when the analysis
# changes. The CACHEDIR can be set for different git branches
FILE = os.path.dirname(__file__)
DATADIR = os.path.join(FILE, '../data/')
CACHEDIR = os.path.join(FILE, '../cache/')
ANALYSISDIR = os.path.join(FILE, '../analysis/')
CACHEFN = os.path.join(CACHEDIR, 'load_table_cache.hdf')
os.system('mkdir -p {}'.format(CACHEDIR)) # creates CACHEDIR if doesn't exist
KBCFN = os.path.join(DATADIR,'kbcvel.csv')
SMEMPLIMIT = 4800

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
            print("read table {} from {}".format(table, cachefn))
            return df
        except IOError:
            print("Could not find cache file: %s" % cachefn)
            print("Building cache...")
            cache = 2
        except KeyError:
            print("Cache not built for table: %s" % table)
            print("Building cache...")
            cache = 2

    if cache == 2:
        df = load_table(table, cache=False)
        print("writing table {} to cache".format(table))
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
        # compute extra flux within circle arcsec due to neighboring stars.
        distarcsec_neighbor = 4 # size of circle
        contrast_bright = 2.5 # threshold contrast for bright companion (mag) 

        # Load up gaia-kepler crossmatch
        fn = os.path.join(DATADIR,'xmatch_gaia2_m17_ruwe_tmass-result.csv')
        df = pd.read_csv(fn)
        df['gaia2_sparallax'] += 0.053 # Apply systematic offset Zinn+2018

        # Neighbor table, all stars within 10 arcsec of gaia-kepler crossmatch
        fn = 'xmatch_gaia2_m17_ruwe_tmass_10arcsec-result.csv'
        fn = os.path.join(DATADIR,fn)
        nei = pd.read_csv(fn)
        nei['distarcsec'] = nei.dist * 3600 # convert degrees to arcsec

        # remove stars outside circle
        nei = nei[nei.distarcsec < distarcsec_neighbor] 

        # identify primary star in neighbor table
        nei = nei.sort_values(by=['id_gaia2_primary','dist'])
        g = nei.groupby('id_gaia2_primary')
        pri = g.first()
        nei = nei.set_index('id_gaia2_primary')

        # compute contrast between all stars and primary (mag)
        nei['gaia2_gmag_contrast'] = nei.gaia2_gmag - pri.gaia2_gmag

        # compute sum of all flux in circle
        pri['gaia2_gflux_dilution'] = (  g['gaia2_gflux'].sum()
                                         / g['gaia2_gflux'].first() )
        pri['gaia2_n_neighbors'] = ( nei[nei.dist > 0]
                                     .groupby('id_gaia2_primary').size() )

        # compute number of birght neighbors 
        b = (nei.dist > 0) & (nei.gaia2_gmag_contrast < contrast_bright)
        pri['gaia2_n_neighbors_bright'] = ( nei[b]
                                            .groupby('id_gaia2_primary')
                                            .size() )
        # some stars have no neighbors within circle so fill in nans
        pri = pri.fillna( dict(gaia2_n_neighbor=0, gaia2_n_neighbor_bright=0) )

        # merge in neighbor columns into gaia-kepler crossmatch
        cols = ['gaia2_gflux_dilution',
                'gaia2_n_neighbors',
                'gaia2_n_neighbors_bright']
        df = pd.merge(df,pri[cols],left_on='id_gaia2',right_index=True)

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

    elif table == 'DR1+DR2+CXM-overlap':
        # Stars that pass the photometric cuts
        cxm = load_table('planets-cuts1',cache=2)
        cxm = cxm[~cxm.isany]
        cxm = cxm[['id_koi']].drop_duplicates()
        cxm['in_cxm'] = True
        dr12 = load_table('DR1+DR2')['id_koi id_name in_dr1 in_dr2 obs_bjd obs_counts'.split()]
        dr12['id_koi'] = dr12.id_koi.replace(False,np.nan)
        cxm['id_koi'] = cxm['id_koi'].astype(float)
        m = pd.merge(cxm, dr12, how='outer', on='id_koi')
        for k in 'in_cxm in_dr1 in_dr2'.split():
            m[k] = m[k].fillna(False)

        df = m.sort_values(by='id_koi')
        
    elif table == 'DR1+DR2' :
        # Stars that pass the photometric cuts
        cxm = load_table('planets-cuts1',cache=2)
        cxm = cxm[~cxm.isany]
        cxm = cxm[['id_koi']].drop_duplicates()
        cxm['in_cxm'] = True

        # Superset of spectroscopic parameters 
        drall = pd.read_csv('data/CKS_Spectroscopic_Parameters.csv')
        drall['id_koi'] = drall[
            drall.id_name.str.contains('^K\d{5}$|^CK\d{5}$')
        ].id_name.str.slice(start=-5).astype(int)

        # Identify the spectra that are in DR1 
        dr1 = pd.read_csv('data/cks1-obs.csv',index_col=None)
        dr1 = (drall.set_index('id_obs')
               .loc[dr1.obs.str.replace('rj','j')]
               .reset_index())
        
        dr1['in_dr1'] = True
        print("{} spectra in dr1".format(len(dr1)))

        # Identify the spectra that are in DR2
        dr2 = drall.set_index('id_name').drop(dr1.id_name).reset_index()
        dr2 = (dr2.set_index('id_koi')
               .drop(dr1.id_koi,errors='ignore')
               .reset_index())
        
        print("{} spectra of stars not in dr1".format(len(dr2)))

        date = '2018-01-01'
        dr2 = dr2[
            (dr2.obs_bjd > Time(date,format='iso').jd)
            | dr2.id_koi.isin(cxm.id_koi)]

        print("{} spectra post {}".format(len(dr2),date))
        dr2 = (dr2.sort_values(by='obs_counts')
               .groupby('id_koi',as_index=False)
               .last())

        print("{} spectra highest snr".format(len(dr2)))
        dr2 = dr2.query('obs_counts > 500')
        print("{} spectra with counts > 500".format(len(dr2)))
        dr2['in_dr2'] = True
        df = pd.concat([dr1,dr2],sort=True)

        # merge in id_kic values
        fn = DATADIR+'cumulative_2020.08.05_15.51.48.csv'
        cm = pd.read_csv(fn,skiprows=53)
        cm = cm.rename(columns={'kepid':'id_kic','kepoi_name':'id_koi'})
        cm = cm['id_kic id_koi'.split()]
        cm['id_koi'] = cm['id_koi'].str.slice(start=1,stop=6).astype(int)
        cm = cm.drop_duplicates()
        df = pd.merge(df,cm,how='left')
        idx = df[df.id_name.str.contains('KIC')].index
        df.loc[idx,'id_kic'] = (df.loc[idx,'id_name']
                                .str.slice(start=3)
                                .astype(int))

        df = df.fillna(False).sort_values(by='id_koi')
        df['id_kic'] = df['id_kic'].astype(int)
       
        
    elif table == 'm17+cdpp+gaia2+ber20':
        m17 = load_table('m17', cache=1)
        # Store to separate cache to prevent long reload times
        cachefn = os.path.join(DATADIR,'cdpp.hdf') 
        cdpp = load_table('cdpp',cache=1, cachefn=cachefn)        
        ber = load_table('berger20', cache=1)
        gaia = load_table('gaia2', cache=1)
        df = pd.merge(m17,cdpp)
        df = pd.merge(df,gaia)
        df = pd.merge(df,ber,on=['id_kic'])

    elif table == 'field-cuts':
        df = load_table('m17+cdpp+gaia2+ber20',cache=1)
        cuttypes = ['none','faint','giantcmd','diluted','ruwe']
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'field')

    elif table == 'planets-cuts1':
        star = load_table('m17+cdpp+gaia2+ber20',cache=1)
        plnt = load_table('koi-thompson18')
        df = pd.merge(star,plnt)
        cuttypes = ['none','faint','giantcmd','diluted','ruwe','notreliable','lowsnr','nomcmc']
        df.sample = 'koi-thompson18'
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'koi-thompson18')

    # All columns that appear in star table
    elif table == 'star':
        # DR1+DR2 
        dr12 = load_table('DR1+DR2')
        nstars1 = len(dr12)
        df2 = []
        for i, row in dr12.iterrows():
            d = dict(row['id_kic id_koi id_name id_obs in_dr1 in_dr2'.split()])
            if row.smemp_teff < SMEMPLIMIT:
                d['cks_steff'] = row.smemp_teff
                d['cks_steff_err'] = row.smemp_teff_err
                d['cks_smet'] = row.smemp_fe
                d['cks_smet_err'] = row.smemp_fe_err
                d['cks_sprov'] = 'emp'
            else:
                d['cks_steff'] = row.smsyn_teff
                d['cks_steff_err'] = row.smsyn_teff_err
                d['cks_slogg'] = row.smsyn_logg
                d['cks_slogg_err'] = row.smsyn_logg_err
                d['cks_smet'] = row.smsyn_fe
                d['cks_smet_err'] = row.smsyn_fe_err
                d['cks_svsini'] = row.smsyn_vsini
                d['cks_sprov'] = 'syn'

            df2.append(d)
            
        df = pd.DataFrame(df2)

        # Add in m17 parameters
        cols = ['id_kic m17_kepmag m17_kmag m17_kmag_err'.split()]
        m17 = load_table('m17', cache=1)
        df = pd.merge(df,m17,how='left')

        # Add in gaia2 parameters
        cols = ['id_kic gaia2_sparallax gaia2_sparallax_err'.split()]
        gaia = load_table('gaia2', cache=1)
        df = pd.merge(df,gaia,how='left')

        # Add in isoclassify parameters
        iso2 = []
        for prov in 'syn emp'.split():
            source = 'sm'+prov
            fn = os.path.join(DATADIR,'isoclassify_{}.csv'.format(source))
            iso = pd.read_csv(fn,index_col=0)
            iso['id_kic'] = iso.id_starname.str.slice(start=3).astype(int)
            iso['cks_sprov'] = prov
            cols = """
            id_kic cks_sprov 
            gdir_srad gdir_srad_err1 gdir_srad_err2
            gdir_avs gdir_avs_err1 gdir_avs_err2
            giso_smass giso_smass_err1 giso_smass_err2
            giso_srad giso_srad_err1 giso_srad_err2
            giso_srho giso_srho_err1 giso_srho_err2
            giso_sage giso_sage_err1 giso_sage_err2
            giso2_sparallax giso2_sparallax_err1 giso2_sparallax_err2
            """.split()
            iso2.append(iso[cols])

        iso2 = pd.concat(iso2,sort=True).reset_index(drop=True)
        iso2['giso_slogage'] = iso2.eval('log10(1e9 * giso_sage)')
        iso2['giso_slogage_err1'] = (
            iso2.eval('log10(giso_sage + giso_sage_err1) - log10(giso_sage)')
        )
        iso2['giso_slogage_err2'] = (
            iso2.eval('log10(giso_sage + giso_sage_err2) - log10(giso_sage)')
        )
        df = pd.merge(df,iso2,how='left',on=['id_kic','cks_sprov'])

        
        # Add in ReaMatch parameters
        rm = load_table('reamatch')
        link = dr12['id_koi id_kic'.split()]
        link['id_koi'] = link.id_koi.astype(int)
        rm = pd.merge(link,rm)['id_kic rm_sb2'.split()] # merge on kic
        df = pd.merge(df,rm,how='left')
        df.loc[df.in_dr1,'rm_sb2'] = 1 # stars in DR1 automatically pass
        nstars2 = len(df)

        assert nstars1==nstars2, "error, table={} should preserve nstars".format(table)
        
    # Stellar sample
    elif table=='planets-cuts1+iso+dr25':
        plnt = load_table('planets-cuts1')
        plnt = plnt[~plnt.isany]
        star = load_table('star')
        
        # select columns needed for updating parameters and making cuts
        cols = """
        cks_steff cks_steff_err 
        cks_smet cks_smet_err 
        cks_svsini rm_sb2
        id_kic cks_sprov 
        gdir_srad gdir_srad_err1 gdir_srad_err2
        gdir_avs gdir_avs_err1 gdir_avs_err2
        giso_smass giso_smass_err1 giso_smass_err2
        giso_srad giso_srad_err1 giso_srad_err2
        giso_srho giso_srho_err1 giso_srho_err2
        giso_sage giso_sage_err1 giso_sage_err2
        giso_slogage giso_slogage_err1 giso_slogage_err2
        giso2_sparallax giso2_sparallax_err1 giso2_sparallax_err2
        """.split()
        star = star[cols]
        df = pd.merge(plnt, star, on=['id_kic'])
        
        # Add in expected transit duration
        df, samp = ckscool.calc.update_planet_parameters(df)

    elif table=='planets-cuts2':
        df = load_table('planets-cuts1+iso+dr25',cache=1)
        cuttypes = [
            'none','badvsini','badspecparallax','sb2',
            'badimpact','badimpacttau','badprad','badpradprec'
        ]
        df = ckscool.cuts.occur.add_cuts(df, cuttypes, 'koi-thompson18')
        
    ############################
    # Tables from other papers #
    ############################
    elif table == 'reamatch':
        fn = os.path.join(DATADIR, 'hires/reamatch.csv')
        df = pd.read_csv(fn,index_col=None, usecols=range(6))
        df['id_koi'] = df.name.str.slice(start=-5).astype(int)
        namemap = {'id_koi':'id_koi','is_sb2':'rm_sb2'}
        df = df.rename(columns=namemap)[namemap.values()]
        assert len(df.id_koi)==len(df.id_koi.drop_duplicates()), "duplicates"
        
    elif table == 'nrm-previous':
        # File is from an email that Adam sent me in 2017-11-02
        fn = os.path.join(DATADIR,'kraus/KOIours.txt')
        df = pd.read_table(fn,header=None,names=['name'])
        df['id_koi'] = df.name.str.slice(start=4).astype(int)
        df = df[['id_koi']]

    elif table == 'fpp':
        fn = os.path.join(DATADIR,'q1_q17_dr25_koifpp.csv')
        df = pd.read_csv(fn,comment='#')
        namemap = {'kepoi_name':'id_koicand','fpp_prob':'fpp_prob'}
        df = df.rename(columns=namemap)[namemap.values()]

    # KOI tables
    elif table=='koi-thompson18':
        df = load_table_koi(table)

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
        cks = load_table('planets-cuts2').groupby('id_koi',as_index=False).nth(0)
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
        cks = load_table('planets-cuts2').groupby('id_koi',as_index=False).nth(0)
        d13 = load_table('dressing13')
        df = pd.merge(cks,d13,on='id_kic')

    # Berger al. (2018)
    elif table == 'berger18':
        fn = os.path.join(DATADIR,'berger18/apjaada83t1_mrt.txt')
        tab = ascii.read(fn)
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

    # Berger al. (2020)
    elif table == 'berger20':
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
        df = add_prefix(df,'ber20_')

    elif table == 'ckscool-brewer18':
        cks = load_table('planets-cuts2').groupby('id_koi',as_index=False).nth(0)
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

    # Silva 2015
    elif table=='silva15':
        fn = os.path.join(DATADIR,'silva15/silva-aguirre15.tex')
        df = read_silva15(fn)

    # Huber 2013
    elif table=='huber13':
        fn = os.path.join(DATADIR,'huber13/J_ApJ_767_127/table2.dat')
        readme = os.path.join(DATADIR,'huber13/J_ApJ_767_127/ReadMe')
        df = read_huber13(fn,readme)

    elif table=='christiansen20':
        df = astropy.io.ascii.read('kplr_dr25_inj1_plti.txt')
        df = df.to_pandas().rename(columns={'KIC_ID':'id_kic'})
        

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

def read_silva15(fn):
    with open(fn,'r') as f:
        lines = f.readlines()
    header = lines[7]

    lines = lines[9:]
    _lines = []
    for line in lines:
        if line.count('&') > 0:
            _lines.append(line)

    lines = _lines

    _lines = []
    i=0
    df = []

    for line in lines:
        d = {}
        line =  line.split('&')
        d = OrderedDict()
        d['id_koi'] = line[0]
        d['id_kic'] = line[1]
        d['teff'] = line[2].split('$\\pm$')[0]
        d['teff_err1'] = line[2].split('$\\pm$')[1]

        d['fe'] = line[3].split('$\\pm$')[0]
        d['fe_err1'] = line[3].split('$\\pm$')[1]

        mass = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[4]).split()
        d['smass'] = mass[0]
        d['smass_err1'] = mass[1]
        d['smass_err2'] = mass[2]

        radius = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[5]).split()
        d['srad'] = radius[0]
        d['srad_err1'] = radius[1]
        d['srad_err2'] = radius[2]

        logg = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[7]).split()
        d['slogg'] = logg[0]
        d['slogg_err1'] = logg[1]
        d['slogg_err2'] = logg[2]

        age = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[9]).split()
        d['sage'] = age[0]
        d['sage_err1'] = age[1]
        d['sage_err2'] = age[2]

        df.append(d)
        i+=1

    df = pd.DataFrame(df).convert_objects(convert_numeric=True)
    df['teff_err2'] = -1.0 * df['teff_err1']
    df['slogage'] = np.log10(df['sage']) + 9 
    df['slogage_err1'] = np.log10(df.sage+df.sage_err1)+9 - df.slogage
    df['slogage_err2'] = np.log10(df.sage+df.sage_err2)+9 - df.slogage
    df['fe_err2'] = -1.0 * df['fe_err1']
    df = add_prefix(df,'s15_')
    return df

def read_huber13(fn, readme):
    df = astropy.io.ascii.read(fn, readme=readme)
    df = df.to_pandas()
    namemap = {
        'KIC':'id_kic',
        'KOI':'id_koi',
        'Mass':'h13_smass',
        'Rad':'h13_srad',
        'e_Mass':'h13_smass_err',
        'e_Rad':'h13_srad_err',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df = df.query('h13_srad > 0.5')
    df.loc[:,'id_kic'] = df.id_kic.astype(int)
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
            'koi_pdisposition',
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
            'koi_fittype'
        ]
        df = df[names]
        csvfn = 'q1_q17_dr25_sup_koi_2020.08.14_12.09.28.csv'
        csvfn = os.path.join(DATADIR,csvfn)
        df2 = pd.read_csv(csvfn,comment='#',index_col=0)
        df = pd.merge(df,df2[['kepoi_name','koi_disposition']],
                      on=['kepoi_name'],how='left',suffixes=['','_sup'])

        
        
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
        print("following columns are not defined")
        print(mcols)

    if not drop:
        df = pd.concat([df,df0[mcols]],axis=1)
    return df


# Code for loading in objects

def load_object(key,cache=0, verbose=1, N_cores=None):
    """
    Objects must adopt the following naming convention
    {objkey}_{params}
    """

    assert key.count('_')==1, "must follow {objkey:}_{params:} convention" 
    objkey, params = key.split('_')
    
    pklfn = os.path.join(CACHEDIR,key+'.pkl')
    if cache == 1:
        try:
            with open(pklfn,'r') as f:
                obj = pickle.load(f)
                if verbose:
                    print("read {} from {}".format(obj,pklfn))
                return obj

        except IOError:
            print("Could not find cache file: %s" % pklfn)
            print("Building cache...")
            cache = 2

    if cache == 2:
        obj = load_object(key, cache=0)
        print("writing {} to {}".format(obj,pklfn))
        with open(pklfn,'w') as f:
            pickle.dump(obj,f)
        return obj
        
    elif objkey=='occur-per-prad' or objkey=='occur-sinc-prad':
        limits = {}
        if params.count('smass'):
            smass1, smass2 = params.replace('smass=','').split('-')
            limits['smass1'] = float(smass1)
            limits['smass2'] = float(smass2)
        obj = load_occur(objkey,limits)

    elif objkey=='comp-per-prad' or objkey=='comp-sinc-prad':
        limits = {}
        if params.count('smass'):
            smass1, smass2 = params.replace('smass=','').split('-')
            limits['smass1'] = float(smass1)
            limits['smass2'] = float(smass2)
        obj = load_comp(objkey,limits)


    elif objkey=='fitper':
        params.split()
        smass1, smass2, prad1, prad2 = (
            params.replace('smass=','').replace('prad=','').split('-')
        ) 
        prad1 = float(prad1)
        prad2 = float(prad2)
        if params.count('prad=1.0-1.7'):
            per1 = 1
            per2 = 30
        elif params.count('prad=1.7-4.0'):
            per1 = 1
            per2 = 300
        elif params.count('prad=1.0-1.5'):
            per1 = 1
            per2 = 30
        elif params.count('prad=2.0-3.0'):
            per1 = 1
            per2 = 300
        else:
            assert False,'key not recognized'
            
        occkey = 'occur-per-prad_smass={}-{}'.format(smass1,smass2)
        occ = load_object(occkey,cache=1)
        fit = ckscool.fit.FitSmoothBrokenPowerLaw(occ, 'per','prad', per1, per2, prad1, prad2)
        fit.create_completeness_spline()

        # Use a non-parameteric approach to get the mean number of
        # planets per star, over the specified range of period. We assume this
        # number is right to a factor of two.

        dlogper = 0.25
        logper = np.arange(np.log10(per1),np.log10(per2),dlogper)
        per = 10**logper
        df = dict(per1=per[:-1], per2=per[1:], prad1=fit.y1, prad2=fit.y2)
        df = pd.DataFrame(df)
        df['perc'] = np.sqrt(df.per1 * df.per2)
        rates = []
        for i, row in df.iterrows():
            rates += [fit.occ.occurrence_box(row)]

        rates = pd.DataFrame(rates)    
        rates = pd.concat([df,rates],ignore_index=False,axis=1)
        f = rates.rate.sum()
        logf = np.log10(f)
        p = lmfit.Parameters()
        p.add('logf',value=logf, min=logf-0.3, max=logf+0.3)
        p.add('k1',value=2, min=0, max=4)
        p.add('k2',value=-1.010, min=-3, max=0,vary=False)
        p.add('logx0',value=1,min=0,max=2)

        # confirm that rate function is normalized
        f = lambda per: fit.rate_lambda(per, p)
        _int, err = scipy.integrate.quad(f, fit.x1, fit.x2)
        assert abs(_int - 1) < 1e-3, "rate function must be normalized"
        fit.fit_max_likelihood(p, method='lbfgsb', nan_policy='omit')
        fit.run_mcmc(short=False)
        obj = fit

    elif objkey=='fitsinc':
        params.split()
        smass1, smass2, prad1, prad2 = (
            params.replace('smass=','').replace('prad=','').split('-')
        ) 
        prad1 = float(prad1)
        prad2 = float(prad2)
        if params.count('prad=1.0-1.7'):
            sinc1 = 10.0
            sinc2 = 3000.0
        elif params.count('prad=1.7-4.0'):
            sinc1 = 1.0
            sinc2 = 3000.0
        else:
            assert False,'key not recognized'

        occkey = 'occur-sinc-prad_smass={}-{}'.format(smass1,smass2)
        occ = load_object(occkey,cache=1)
        fit = ckscool.fit.FitSmoothBrokenPowerLaw(
            occ, 'sinc', 'prad', sinc1, sinc2, prad1, prad2
        )
        fit.create_completeness_spline()

        # Use a non-parameteric approach to get the mean number of
        # planets per star, over the specified range of period. We assume this
        # number is right to a factor of two.

        dlogsinc = 0.25
        eps = 1e-10
        logsinc = np.arange(np.log10(sinc1),np.log10(sinc2)+eps,dlogsinc)
        sinc = 10**logsinc
        df = dict(sinc1=sinc[:-1], sinc2=sinc[1:], prad1=fit.y1, prad2=fit.y2)
        df = pd.DataFrame(df)
        df['sincc'] = np.sqrt(df.sinc1 * df.sinc2)
        rates = []
        for i, row in df.iterrows():
            rates += [fit.occ.occurrence_box(row)]

        rates = pd.DataFrame(rates)    
        rates = pd.concat([df,rates],ignore_index=False,axis=1)
        f = rates.rate.sum()
        logf = np.log10(f)
        p = lmfit.Parameters()
        p.add('logf', value=logf, min=logf-0.3, max=logf+0.3)
        p.add('k1', value=-1.0001, min=-2, max=0, vary=False)
        p.add('k2', value=-3.0, min=-5, max=-2,vary=True)
        p.add('logx0', value=2.0, min=0, max=3)

        # confirm that rate function is normalized
        f = lambda x: fit.rate_lambda(x, p)
        _int, err = scipy.integrate.quad(f, fit.x1, fit.x2)
        assert abs(_int - 1) < 1e-3, "rate function must be normalized"
        fit.fit_max_likelihood(p, method='lbfgsb', nan_policy='omit')
        fit.run_mcmc(short=False)
        obj = fit

    elif objkey.count('grad'):
        obj = ckscool.gradient.load_gradient_chain(key)
    
    elif objkey=='mps':
        if params.count('size-se'):
            mps = ckscool.occur.MeanPlanetSize('per < 30 and 1.0 < prad < 1.7')
        elif params.count('size-sn'):
            mps = ckscool.occur.MeanPlanetSize('per < 100 and 1.7 < prad < 4.0')
        else:
            assert False, "Failed"

        mps.sample_logprad()
        mps.fit_powerlaw()
        obj = mps

    elif key.count('fitdetected_')==1:
        mode = key.split('_')[-1]
        fitter = ckscool.fitdetected.Fitter(mode)
        fitter.compute_samples()
        obj = fitter

    else:
        assert False, "could not resolve key"

    return obj

def load_comp(objkey, limits):

    # Derive completeness object
    method = 'christiansen20-gamma-clip' # treatment for planet detectability
    impact = 0.8 # maximum impact parameter considered.

    field = ckscool.io.load_table('field-cuts',cache=1)
    field = field[~field.isany]
    field = field.rename(columns={'ber20_srad':'srad','ber20_smass':'smass'})

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']
        field = field[field.smass.between(smass1,smass2)]

    n1 = len(field)
    field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
    nstars = len(field)
    print("{}/{} stars remain after droping nulls".format(nstars,n1))
    if objkey=='comp-per-prad':
        xbins = np.round(np.logspace(np.log10(0.1),np.log10(1000),65),4)
        ybins = np.round(np.logspace(np.log10(0.25),np.log10(64),51 ),2)
        xk = 'per'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessPerPrad

    elif objkey=='comp-sinc-prad':
        xbins = np.round(np.logspace(np.log10(0.1),np.log10(100000),65),4)
        ybins = np.round(np.logspace(np.log10(0.25),np.log10(64),51 ),2)
        xk = 'sinc'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessSincPrad

    else:
        assert False, "{} not supported objkey".format(objkey)

    comp_bins_dict = {xk: xbins,yk: ybins}
    spacing_dict = {xk:xscale, yk:yscale}
    grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)
    comp = Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=True)
    comp.compute_grid_prob_tr(verbose=True)
    comp.create_splines()

    return comp

def load_occur(objkey, limits, debug=False):
    """
    Constructs occurrence object
    """

    # Derive completeness object
    method = 'christiansen20-gamma-clip' # treatment for planet detectability
    impact = 0.8 # maximum impact parameter considered.

    field = ckscool.io.load_table('field-cuts',cache=1)
    field = field[~field.isany]
    field = field.rename(columns={'ber20_srad':'srad','ber20_smass':'smass'})
    plnt = ckscool.io.load_table('planets-cuts2')
    plnt = plnt[~plnt.isany]

    namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass',
               'giso_sinc':'sinc'}
    plnt = plnt.rename(columns=namemap)

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']
        field = field[field.smass.between(smass1,smass2)]
        plnt = plnt[plnt.smass.between(smass1,smass2)]

    n1 = len(field)
    field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
    nstars = len(field)

    print("{}/{} stars remain after droping nulls".format(nstars,n1))

    if objkey=='occur-per-prad':
        xbins = np.round(np.logspace(np.log10(0.1),np.log10(1000),65),4)
        ybins = np.round(np.logspace(np.log10(0.25),np.log10(64),51 ),2)
        xk = 'per'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessPerPrad
        Occurrence = ckscool.occur.OccurrencePerPrad

    elif objkey=='occur-sinc-prad':
        xbins = np.round(np.logspace(np.log10(0.1),np.log10(100000),65),4)
        ybins = np.round(np.logspace(np.log10(0.25),np.log10(64),51 ),2)
        xk = 'sinc'
        yk = 'prad'
        xscale = 'log'
        yscale = 'log'
        Completeness = ckscool.comp.CompletenessSincPrad
        Occurrence = ckscool.occur.OccurrenceSincPrad

    else:
        assert False, "{} not supported objkey".format(objkey)

    if debug:
        xbins = xbins[::2]
        ybins = ybins[::2]

    comp_bins_dict = {xk: xbins,yk: ybins}
    spacing_dict = {xk:xscale, yk:yscale}
    grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)
    comp = Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=True)
    comp.compute_grid_prob_tr(verbose=True)
    comp.create_splines()
    occ = Occurrence(plnt, comp, nstars)
    return occ
