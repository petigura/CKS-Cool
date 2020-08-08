import os

import numpy as np
import pandas as pd
import glob 
import ckscool.io

def load_stellar_parameters(source):
    """
    Load up stellar properties
    """
    if source=='cks1':
        df = pd.read_csv('data/cks_physical_merged.csv')
        df['id_koi'] = df.id_koicand.str.slice(start=1,stop=6).astype(int)
        df = df.groupby('id_koi',as_index=False).nth(0)
        df['cks_steff_err'] = 100 # stat+sys uncert
        df['cks_slogg_err'] = 0.1 # stat+sys uncert
        df['cks_smet_err'] = 0.06 # stat+sys uncert
        df['cks_svsini_err'] = df['cks_svsini_err1']
        df['cks_sprov'] = 'cks1'
        cols = [
            'id_kic','id_koi','cks_steff','cks_steff_err','cks_slogg',
            'cks_slogg_err','cks_smet',
            'cks_smet_err','cks_svsini','cks_svsini_err'
        ]
        df = df[cols]

    # Load up SpecMatch-Syn and set uncertainties 
    elif source=='smsyn':
        df = ckscool.io.load_table('DR1+DR2')
        namemap = {
            'id_kic':'id_kic',
            'id_obs':'id_obs',
            'id_name':'id_name',
            'smsyn_teff':'cks_steff',
            'smsyn_logg':'cks_slogg',
            'smsyn_fe':'cks_smet',
            'smsyn_vsini':'cks_svsini',
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df['cks_sprov'] = 'smsyn'
        df['cks_steff_err'] = 100 # stat+sys uncert
        df['cks_slogg_err'] = 0.1 # stat+sys uncert
        df['cks_smet_err'] = 0.06 # stat+sys 

    # Load SpecMatch-Emp and set uncertainties
    elif source=='smemp':
        df = ckscool.io.load_table('DR1+DR2')
        namemap = {
            'id_kic':'id_kic',
            'id_obs':'id_obs',
            'id_name':'id_name',
            'smemp_teff':'cks_steff',
            'smemp_teff_err':'cks_steff_err',
            'smemp_fe':'cks_smet',
            'smemp_fe_err':'cks_smet_err',
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df['cks_sprov'] = 'smemp'
        df['cks_slogg'] = -99
        df['cks_slogg_err'] = -99
        df['cks_steff_err'] = 60 # valid for cool stars
        df['cks_smet_err'] = 0.12 # valid for cool stars
    else:
        assert False, "invalid mode"

    star = ckscool.io.load_table('m17+cdpp+gaia2+ber19')
    star = pd.merge(star,df,on='id_kic')
    star['m17_kmag_err'] = star['m17_kmag_err'].fillna(0.02)
    star = star.sort_values(by='id_kic')
    star = star.reset_index()
    star0 = star.copy() 
    return star0

def create_iso_batch_frames(source):
    """Create Isoclassify Batch Jobs

    Creates input parameters for two runs
    
       1. The direct method with the following constraints
          - teff, logg, fe, parallax, kmag

       2. The grid method with the following constraints
          - teff, logg, met, kmag [and parallax]

       3. The grid method with the following constraints
          - teff, logg, met, kmag [no parallax]

    We default to the direct method. But if the parallax method from
    the grid based method is significantly different than the gaia
    parallax, there is additional flux in the aperture which indicates
    dilution.

    Args:

        source (string): either
            cks1
            smemp
            smsyn

    """
    star = load_stellar_parameters(source)
    star['id_starname'] = star.id_kic.apply(lambda x : "KIC{}".format(x) )
    star['band'] = 'kmag'
    star['dust'] = 'green18'
    star['cks_sprov'] = source
    namemap = {
        'gaia2_ra':'ra',
        'gaia2_dec':'dec',
        'gaia2_sparallax':'parallax',
        'gaia2_sparallax_err':'parallax_err',
        'm17_kmag':'kmag',
        'm17_kmag_err':'kmag_err',
        'cks_steff':'teff',
        'cks_steff_err':'teff_err',
        'cks_slogg':'logg',
        'cks_slogg_err':'logg_err',
        'cks_smet':'feh',
        'cks_smet_err':'feh_err',
    }
    star = star.rename(columns=namemap)
    star['parallax'] /= 1e3 # Convert microarcsec to arcsec
    star['parallax_err'] /= 1e3 
    star0 = star.copy()

    cols = [
        'id_starname','teff','teff_err','logg','logg_err','feh','feh_err',
        'parallax','parallax_err','kmag','kmag_err',
        'ra','dec','band','dust','cks_sprov'
    ]

    # Direct method. Don't use spectroscopic logg values so as to not
    # pollute the parallax radii
    star = star0.copy()
    star['logg'] = -99 # use large uncertainties
    star['logg_err'] = -99 
    star = star[cols]
    star_direct = star.copy()

    # Grid method with parallax. This will return model-dependent
    # values of Mstar, Rstar, age, density, luminosity
    star = star0.copy()
    star['logg'] = -99 # use large uncertainties
    star['logg_err'] = -99 
    star = star[cols]
    star_grid_yes = star.copy()

    # Grid method. Don't set parallax so we can compare later
    star = star0.copy()
    star = star[cols]
    star['parallax'] = -99
    star['parallax_err'] = -99
    star_grid_no = star.copy()
    return star_direct, star_grid_yes, star_grid_no 

def create_iso_table(inpdir,outcsv):
    """
    Read in isochrones csvfiles 

    Args:
        outdir (str): where to look for isochrones.csv files
    """

    fn = os.path.join(inpdir,'direct/*/*.csv')
    dfd = scrape_direct(fn)

    fn = os.path.join(inpdir,'grid-parallax-yes/*/*.csv')
    dfg = scrape_grid_parallax_yes(fn)

    fn = os.path.join(inpdir,'grid-parallax-no/*/*.csv')
    dfg2 = scrape_grid_parallax_no(fn)

    dfm = pd.merge(dfd,dfg,on='id_starname', how='left')
    dfm = pd.merge(dfm,dfg2,on='id_starname', how='left')
    temp = dfm['id_starname'].copy()
    dfm = dfm.drop(['id_starname'],axis=1)
    dfm = dfm.convert_objects(convert_numeric=True)
    dfm['id_starname'] = temp.astype(str)
    #dfm = ckscool.io.order_columns(dfm)
    
    dfm.to_csv(outcsv)
    print "created {}".format(outcsv)

def func(x):
    return x.split('.')[0]

def scrape_direct(fn):
    # isoclassify/direct/*/*.csv
    df = scrape_csv(fn)
    df['id_starname'] = df.id_starname.astype(str).apply(func)
    namemap = {
        'id_starname':'id_starname',
        'dir_rad':'gdir_srad',
        'dir_rad_err1':'gdir_srad_err1',
        'dir_rad_err2':'gdir_srad_err2',
        'dir_avs':'gdir_avs',
        'dir_avs_err1':'gdir_avs_err1',
        'dir_avs_err2':'gdir_avs_err2',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    return df

def scrape_grid_parallax_yes(fn):
    # Grid mode with parallax constraints
    df = scrape_csv(fn)
    namemap = {
        'id_starname':'id_starname',
        'iso_mass':'giso_smass',
        'iso_mass_err1':'giso_smass_err1',
        'iso_mass_err2':'giso_smass_err2',
        'iso_rad':'giso_srad',
        'iso_rad_err1':'giso_srad_err1',
        'iso_rad_err2':'giso_srad_err2',
        'iso_rho':'giso_srho',
        'iso_rho_err1':'giso_srho_err1',
        'iso_rho_err2':'giso_srho_err2',
        'iso_age':'giso_sage',
        'iso_age_err1':'giso_sage_err1',
        'iso_age_err2':'giso_sage_err2',
    }
    df = _rename(df, namemap)
    df['id_starname'] = df.id_starname.astype(str).apply(func)
    return df

def scrape_grid_parallax_no(fn):
    # Grid mode without parallax constraints
    df = scrape_csv(fn)
    temp = df['id_starname'].copy()
    df = df.drop(['id_starname'],axis=1)
    df = df.convert_objects(convert_numeric=True)
    df['id_starname'] = temp.astype(str)
    df['giso2_sparallax'] = 1 / df.iso_dis * 1e3
    df['giso2_sparallax_err1'] = - df['giso2_sparallax'] * df['iso_dis_err2'] / df['iso_dis']
    df['giso2_sparallax_err2'] = - df['giso2_sparallax'] * df['iso_dis_err1'] / df['iso_dis']
    columns = ['id_starname','giso2_sparallax','giso2_sparallax_err1', 'giso2_sparallax_err2',]
    df = df[columns]
    return df

def _rename(df, namemap):
    """
    Perform the rename, if no rows exist, return empty dataframe
    """
    if len(df)==0:
        return pd.DataFrame(columns=namemap.values())

    df = df.rename(columns=namemap)[namemap.values()]
    return df

# Lifted from isoclassify module to remove dependence
def _csv_reader(f):
    row = pd.read_csv(f,header=None,squeeze=True, index_col=0)
    return row

def scrape_csv(path):
    """
    Read in isochrones csvfiles 
    Args:
        outdir (str): where to look for isochrones.csv files
    """
    fL = glob.glob(path)
    df = []

    for i, f in enumerate(fL):
        if i%100==0:
            print(i)
        try:
            df.append(_csv_reader(f))
        except ValueError:
            print("{} failed".format(f))


    df = pd.DataFrame(df)
    df = df.reset_index()
    return df
