import ckscool.io
import numpy as np
import pandas as pd

class Table(object):
    """
    Helper to create ascii tables
    """
    def to_latex(self,cols,header=True):
        # create formatted string
        ncols = len(cols)
        lines = []
        if header:
            s = ' & '.join(cols)
            s = '%' + s 
            lines.append(s)

        for irow,row in self.df.iterrows():
            s = r""
            for icol,c in enumerate(cols):
                s+=('{:' +self.formats[c] + '}').format(row[c])
                if icol==(ncols-1):
                    s+=' \\\\ '
                else:
                    s+=' & '

            s = s.replace('nan','\\nodata')
            lines.append(s)

        return lines

    def to_csv(self,cols,header=True):
        # create formatted string
        ncols = len(cols)
        lines = []
        if header:
            s = ','.join(cols)
            lines.append(s)

        for irow,row in self.df.iterrows():
            s = r""
            for icol,c in enumerate(cols):
                s+=('{:' +self.formats[c] + '}').format(row[c])
                if icol==(ncols-1):
                    s+=''
                else:
                    s+=','

            lines.append(s)

        return lines

class TableStar(Table):
    def __init__(self):
        df = ckscool.io.load_table('star',cache=1)
        m = ckscool.io.load_table('DR1+DR2+CXM-overlap')
        df = pd.merge(df,m[['id_name','in_cxm']], on='id_name')
        df['id_koi'] = df['id_koi'].replace(False,np.nan)
        df = df.sort_values(by='id_koi')
        self.df = df
        self.formats = {
            'id_koi':'0.0f',
            'id_gaia2':'f',
            'm17_kmag':'0.2f',
            'm17_kmag_err':'0.2f',
            'gaia2_sparallax':'0.2f',
            'gaia2_sparallax_err':'0.2f',
            'cks_steff':'0.0f',
            'cks_steff_err':'0.0f',
            'cks_smet':'0.2f',
            'cks_smet_err':'0.2f',
            'cks_svsini':'0.1f',
            'cks_svsini_err':'0.1f',
            'cks_sprov':'s',
            'gdir_srad':'0.2f',
            'gdir_srad_err1':'0.2f',
            'gdir_srad_err2':'0.2f',
            'giso_smass':'0.2f',
            'giso_smass_err1':'0.2f',
            'giso_smass_err2':'0.2f',
            'giso_srad':'0.2f',
            'giso_srad_err1':'0.2f',
            'giso_srad_err2':'0.2f',
            'giso_srho':'0.2f',
            'giso_srho_err1':'0.2f',
            'giso_srho_err2':'0.2f',
            'giso_sage':'0.1f',
            'giso_sage_err1':'0.1f',
            'giso_sage_err2':'0.1f',
            'giso2_sparallax':'0.2f',
            'giso2_sparallax_err1':'0.2f',
            'giso2_sparallax_err2':'0.2f',
            'rm_sb2':'0.0f',
            'in_cxm':'0.0f'
        }
     

def tab_star():
    t = TableStar()
    cols = """id_koi m17_kmag gaia2_sparallax cks_steff cks_smet cks_svsini cks_sprov gdir_srad giso_smass giso_srad giso_srho giso_sage giso2_sparallax rm_sb2 in_cxm""".split()
    return t.to_latex(cols)

def tab_star_csv():
    t = TableStar()
    cols = """id_koi id_gaia2 m17_kmag m17_kmag_err gaia2_sparallax gaia2_sparallax_err
cks_steff cks_steff_err cks_smet cks_smet_err cks_svsini gdir_srad
gdir_srad_err1 gdir_srad_err2 giso_smass giso_smass_err1
giso_smass_err2 giso_srad giso_srad_err1 giso_srad_err2 giso_srho
giso_srho_err1 giso_srho_err2 giso2_sparallax giso2_sparallax_err1
giso2_sparallax_err2 cks_sprov rm_sb2 in_cxm"""
    cols = cols.split()
    return t.to_csv(cols)

    
class TablePlanet(Table):
    def __init__(self):
        df = ckscool.io.load_table('planets-cuts2',cache=1)
        df = df.sort_values(by='id_koicand')
        df['in_curated'] = ~df.isany
        df['dr25_ror']*=100
        self.df = df
        self.formats = {
            'id_koicand':'s',
            'koi_period':'0.6f',
            'koi_period_err1':'0.6f',
            'koi_period_err2':'0.6f',
            'dr25_ror':'0.2f',
            'dr25_ror_err1':'0.2f',
            'dr25_ror_err2':'0.2f',
            'dr25_tau':'0.2f',
            'dr25_tau_err1':'0.2f',
            'dr25_tau_err2':'0.2f',
            'gdir_prad':'0.2f',
            'gdir_prad_err1':'0.2f',
            'gdir_prad_err2':'0.2f',
            'giso_sma':'0.3f',
            'giso_sma_err1':'0.3f',
            'giso_sma_err2':'0.3f',
            'giso_tau0':'0.2f',
            'giso_tau0_err1':'0.2f',
            'giso_tau0_err2':'0.2f',
            'giso_sinc':'0.2f',
            'giso_sinc_err1':'0.2f',
            'giso_sinc_err2':'0.2f',
            'in_curated':'0.0f',
        }

def tab_planet():
    t = TablePlanet()
    t.formats['dr25_ror'] = '0.2f'
    t.formats['koi_period'] = '0.1f'
    cols = 'id_koicand koi_period dr25_ror dr25_tau gdir_prad giso_tau0 giso_sma giso_sinc in_curated'.split()
    return t.to_latex(cols)

def tab_planet_csv():
    t = TablePlanet()
    cols = 'id_koicand koi_period koi_period_err1 koi_period_err2 dr25_ror dr25_ror_err1 dr25_ror_err2 dr25_tau dr25_tau_err1 dr25_tau_err2 gdir_prad gdir_prad_err1 gdir_prad_err2 giso_tau0 giso_tau0_err1 giso_tau0_err2 giso_sma giso_sma_err1 giso_sma_err2 giso_sinc giso_sinc_err1 giso_sinc_err2 in_curated'.split()
    return t.to_csv(cols)

