import ckscool.io
import pandas as pd
import isoclassify.pipeline
import ckscool.cuts.occur
import numpy as np

def load_iso_batch_table():
    # Spectra
    plnt = ckscool.io.load_table('koi-thompson18-dr25')
    cuttypes = ['faint','giant','rizzuto','notreliable','lowsnr']
    plnt.sample = 'koi-thompson18'
    plnt = ckscool.cuts.occur.add_cuts(plnt, cuttypes, 'koi-thompson18')
    plntc = plnt.query('isany==False')
    star = plntc.groupby('id_koi', as_index=False).nth(0)
    star = star['id_koi id_koicand id_kic kepmag ber18_srad ber18_steff m17_kmag m17_kmag_err gaia2_sparallax gaia2_sparallax_err gaia2_ra gaia2_dec'.split()]
    star['steff'] = np.nan
    star['steff_err'] = np.nan
    star['smet'] = np.nan
    star['smet_err'] = np.nan
    star['svsini'] = np.nan
    star['svsini_err'] = np.nan
    star['sprov'] = None
    star.index=star.id_koi

    # Load SpecMatch-Emp and set uncertainties
    kbc = ckscool.io.load_table('kbc')
    kbc = kbc.groupby('id_koi', as_index=False).nth(-1)
    df = pd.read_csv('data/specmatch-emp_results.csv')
    df = df.dropna(subset=['name'])
    namemap = {
       'obs':'id_obs',
       'name':'id_name',
       'teff':'steff',
       'teff_err':'steff_err',
       'fe':'smet',
       'fe_err':'smet_err',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df['sprov'] = 'emp'
    df['steff_err'] = 60
    df['smet_err'] = 0.12
    sme = pd.merge(kbc, df, on=['id_obs','id_name'], how='left')
    sme.index=sme.id_koi

    # Load up CKS-I and set uncertainties
    cks1 = pd.read_csv('data/cks_physical_merged.csv')
    cks1['id_koi'] = cks1.id_koicand.str.slice(start=1,stop=6).astype(int)
    cks1 = cks1.groupby('id_koi').nth(0)
    namemap = {
       'cks_steff':'steff',
       'cks_smet':'smet',
       'cks_svsini':'svsini',
    }
    cks1 = cks1.rename(columns=namemap)
    cks1['steff_err'] = 100
    cks1['smet_err'] = 0.06
    cks1['sprov'] = 'cks1'

    # Load up SpecMatch-Syn and set uncertainties 
    df = pd.read_csv('data/specmatch-syn_results.csv')
    df = df.dropna(subset=['name'])
    namemap = {
       'obs':'id_obs',
       'name':'id_name',
       'teff':'steff',
       'teff_err':'steff_err',
       'fe':'smet',
       'fe_err':'smet_err',
       'vsini':'svsini',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df['sprov'] = 'syn'
    df['steff_err'] = 100
    sms = pd.merge(kbc, df, on=['id_obs','id_name'], how='left')
    sms.index=sms.id_koi

    # Cut off for using emp
    idxsmemp = sme.query('steff < 4700').index
    idxsmsyn = sme.query('steff > 4700').index

    star.fillna(value=sme.loc[idxsmemp],inplace=True) # Cool stars fill in with emp
    star.fillna(value=cks1,inplace=True) # Stars that are in CKS-I fill, also take vsini for everythinh
    star.fillna(value=sms.loc[idxsmsyn],inplace=True) #
    star['sprov'].fillna(value='None',inplace=True)
    star['ber18_srad'] = star.ber18_srad.astype(float)
    return star

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
    star = load_iso_batch_table()

    # Direct method with parallax constraints
    star = star.rename(
        columns={
            'gaia2_ra':'ra',
            'gaia2_dec':'dec',
            'gaia2_sparallax':'parallax',
            'gaia2_sparallax_err':'parallax_err',
            'm17_kmag':'kmag',
            'm17_kmag_err':'kmag_err',
            'steff':'teff',
            'steff_err':'teff_err',
            'smet':'feh',
            'smet_err':'feh_err',
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
