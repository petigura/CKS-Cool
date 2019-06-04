import ckscool.io
import pandas as pd
import isoclassify.pipeline
import ckscool.cuts.occur
import numpy as np


def create_iso_batch(args):
    """Create Isoclassify Batch Jobs

    Creates input parameters for two runs
    
       1. The direct method with the following constraints
          - teff, logg, fe, parallax, kmag

       2. The grid method with the following constraints
          - teff, logg, met, kmag [no parallax]

    We default to the direct method. But if the parallax method from
    the grid based method is significantly different than the gaia
    parallax, there is additional flux in the aperture which indicates
    dilution.

    """
    df = load_iso_batch_table()
    df = df.reset_index(drop=True)
    star =  df.groupby('id_koi',as_index=False).nth(0)

    # Direct method with parallax constraints
    star = star.rename(
        columns={
            'gaia2_ra':'ra',
            'gaia2_dec':'dec',
            'gaia2_sparallax':'parallax',
            'gaia2_sparallax_err':'parallax_err',
            'm17_kmag':'kmag',
            'm17_kmag_err':'kmag_err',
            'cks_steff':'teff',
            'cks_steff_err':'teff_err',
            'cks_smet':'feh',
            'cks_smet_err':'feh_err',
        }
    )
    star['id_starname'] = star.id_koi.apply(lambda x : "K{:05d}".format(x))
    star['kmag_err'] = star['kmag_err'].fillna(0.02)
    star['band'] = 'kmag'
    star['dust'] = 'green18'
    star['parallax'] /= 1e3 # Convert microarcsec to arcsec
    star['parallax_err'] /= 1e3 
    star['feh_err'] = 0.12 # Ditto

    star0 = star.copy() 

    cols = [
        'id_starname','teff','teff_err','logg','logg_err','feh','feh_err',
        'parallax','parallax_err','kmag','kmag_err',
        'ra','dec','band','dust'
    ]

    # Direct method. Don't use spectroscopic logg values so as to not
    # pollute the parallax radii
    star = star0.copy()
    star['logg_err'] = 1 # use large uncertainties
    star['logg'] = 4.7 
    star = star[cols]
    fn = 'data/isoclassify-direct.csv'
    star.to_csv(fn)
    print "created {}".format(fn)

    # Grid method with parallax. This will return model-dependent
    # values of Mstar, Rstar, age, density, luminosity
    star = star0.copy()
    star['logg_err'] = 1 # use large uncertainties
    star['logg'] = 4.7 
    star = star[cols]
    fn = 'data/isoclassify-grid-parallax-yes.csv'
    star.to_csv(fn)
    print "created {}".format(fn)

    # Grid method. Don't set parallax so we can compare later
    star = star0.copy()
    star['logg_err'] = 1 # use large uncertainties
    star['logg'] = 4.7
    star = star[cols]
    star['parallax'] = -99
    star['parallax_err'] = 0
    fn = 'data/isoclassify-grid-parallax-no.csv'
    star.to_csv(fn)
    print "created {}".format(fn)


def create_iso_table(args):
    """
    Read in isochrones csvfiles 
    Args:
        outdir (str): where to look for isochrones.csv files
    """
    dfd = scrape_direct()
    dfg = scrape_grid_parallax_yes()
    dfg2 = scrape_grid_parallax_no()

    dfm = pd.merge(dfd,dfg,on='id_starname', how='left')
    dfm = pd.merge(dfm,dfg2,on='id_starname', how='left')
    temp = dfm['id_starname'].copy()
    dfm = dfm.drop(['id_starname'],axis=1)
    dfm = dfm.convert_objects(convert_numeric=True)
    dfm['id_starname'] = temp.astype(str)
    #dfm = ckscool.io.order_columns(dfm)
    
    fn = 'data/isoclassify_gaia2.csv'
    dfm.to_csv(fn)
    print "created {}".format(fn)

def func(x):
    return x.split('.')[0]

def scrape_direct():
    df = isoclassify.pipeline.scrape_csv('isoclassify/direct/*/*.csv')
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

def scrape_grid_parallax_yes():
    # Grid mode with parallax constraints
    fn = 'isoclassify/grid-parallax-yes/*/*.csv'
    df = isoclassify.pipeline.scrape_csv(fn)
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

def scrape_grid_parallax_no():
    # Grid mode without parallax constraints
    fn = 'isoclassify/grid-parallax-no/*/*.csv'
    df = isoclassify.pipeline.scrape_csv(fn)
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
