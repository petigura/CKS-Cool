"""
High-level analysis routines for CKS-Metallicity 
"""
import re

import numpy as np
import pandas as pd 

import cksmet.io
import cksmet.grid
import cksmet.comp
import cksmet.occur
import cksmet.fit
from matplotlib import pylab as plt

# Some nice, convenient grids
per_bins_dict = {
    # sqrt2 
    'xfine': np.round(np.logspace(np.log10(0.1),np.log10(1000),33),4),
    # sqrt2 
    'fine': [ 
        1.00, 1.41, 2.00,  2.83,  4.00, 5.66, 8.00,  
        11.3, 16., 22.6, 32.0, 45.3, 64.0, 90.5, 128., 
        181,  256 ],
    'coarse': [ 
        1.00, 2.00,  4.00, 8.00,  
        16., 32.0, 64.0, 128., 256 ],
    'four-per-decade': 10**np.arange(-0.5,3.001,0.25)    
}    

# Neceassary to overshoot edge oc completeness zone so that spline
# interpolation works so eventhough we only go up to 32 Re, we compute
# completeness out to 64
prad_bins_dict = {
    'xfine': np.round(np.logspace(np.log10(0.25),np.log10(64),49 ),2),
    'fine': np.round(np.logspace(np.log10(0.25),np.log10(64),25 ),2),
    'two-per-octave': 2**np.arange(-2,6.001,0.5),
    'physical': [1.0,1.7,4.0,8.0,24.0],
    'coarse': [ 
         0.5,  0.71, 1.00, 1.41, 2.0, 2.83, 4.00, 5.66, 8.0, 11.31, 16.0
    ]
}

# Commonly used qualitative indecies
dist = pd.DataFrame(
    index=[1.0,10,300], data=['hot','warm','cool'], columns=['dist']
)
size = pd.DataFrame(
    index=[1.0,1.7,4.0,8.0],data=['se','sn','ss','jup',],columns=['size']
)
smet = pd.DataFrame(index=[-1,0],data=['sub','sup'],columns=['smet'])

def load_completeness():
    method = 'fulton-gamma' # treatment for planet detectability
    impact = 0.9 # maximum impact parameter considered.

    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'m17_srad':'srad','m17_smass':'smass'})

    # Some photometric properties have null values, drop
    n1 = len(field)
    field = field.dropna(subset=cksmet.comp.__STARS_REQUIRED_COLUMNS__)
    n2 = len(field)
    print "{}/{} stars remain after droping nulls ".format(n2,n1)

    comp_per_bins = per_bins_dict['xfine']
    comp_prad_bins = prad_bins_dict['xfine']
    comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.grid.Grid(comp_bins_dict,spacing_dict)

    comp = cksmet.comp.Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=True)
    comp.compute_grid_prob_tr(verbose=True)
    comp.create_splines()
    return comp

def compute_occur_bins(per_bins, prad_bins, smet_bins, verbose=1):
    """
    Compute occurrence over speficied grid of period, radius, and metallicity
    """

    nstars = cksmet.cuts.n_stars_field_pass_eff 

    print "Initializing occurrence object"
    comp = cksmet.io.load_object('comp',cache=1)
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks[~cks.isany]
    cks = cks.dropna(subset=['iso_prad'])
    namemap = {
        'iso_srad':'srad','koi_period':'per','iso_prad':'prad',
        'iso_sma':'smax','cks_smet':'smet','koi_impact':'impact',
        'koi_max_mult_ev':'mes'
    }
    plnt = cks.rename(columns=namemap)

    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    smet_field = lamo.lamo_smet
    occur = cksmet.occur.Occurrence(plnt, comp, nstars, smet_field=smet_field)

    print "Define occurrence grid"

    bins_dict = {'per':per_bins, 'prad': prad_bins, 'smet':smet_bins}
    spacing_dict = {'per':'log', 'prad':'log', 'smet':'linear'}
    grid = cksmet.grid.Grid(bins_dict,spacing_dict)

    print "Loop over occurrence grid"
    out = []
    i_count = 0
    def _print_row(d):
        s = ""
        s+="per={per1:.1f}-{per2:.1f}, "
        s+="prad={prad1:.1f}-{prad2:.1f}, "
        s+="smet={smet1:.1f}-{smet2:.1f}, "
        s+="rate={rate_str:s} "
        s = s.format(**_out)
        print s 

    for i, row in grid.ds.to_dataframe().iterrows():
        row = dict(row)
        _out = occur.occurence_box(row)
        _out = dict(row, **_out)
        out.append(_out)
        i_count+=1
        if (i_count%10==0) and (verbose==1):
            _print_row(out)

    out = pd.DataFrame(out)
    for k in 'per1 perc per2 prad1 pradc prad2 smet1 smetc smet2'.split():
        out[k+'r'] = out[k].round(decimals=1)

    occur.df = out
    return occur



def set_index(occ, mode):
    df = occ.df
    if mode=='dist-size':
        df = pd.merge(df, dist,left_on='per1',right_index=True)
        df = pd.merge(df, size,left_on='prad1',right_index=True)
        df.set_index(['dist', 'size'], inplace=True)

    if mode=='size-smet':
        df = pd.merge(df, size,left_on='prad1',right_index=True)
        df = pd.merge(df, smet,left_on='smet1',right_index=True)
        df.set_index(['size', 'smet'], inplace=True)

    if mode=='size':
        df = pd.merge(df, size,left_on='prad1',right_index=True)
        df.set_index(['size'], inplace=True)

    occ.df = df
    return occ
    
def load_occur(key):
    print "Loading occurrence object {}".format(key)
    eps = 1e-3

    def _per_bins(per):
        per = per.replace('per=','')
        per = float(per)
        per_bins = 10**(np.arange(np.log10(1),np.log10(350)+eps,per))
        return per_bins

    def _smet_bins(smet):
        smet = smet.replace('smet=','')
        if smet=='all':
            smet_bins=[-1,0.5]
        elif smet=='sub':
            smet_bins=[-1,0]
        elif smet=='sup':
            smet_bins=[0,0.5]
        else:
            smet = float(smet)
            smet_bins = np.arange(-0.4,0.4+eps,smet)
        return smet_bins

    # Same as above but with indexes for the period grid (hot, warm)
    if re.match('occur-per=hotwarm-prad=physical-smet=(.*)', key):
        _, per, prad, smet = key.split('-')
        per_bins = [1, 10, 100]
        prad_bins = prad_bins_dict['physical']
        smet_bins = _smet_bins(smet)
        occ = compute_occur_bins(per_bins, prad_bins, smet_bins)
        occ = set_index(occ, 'dist-size')

    # Occurence log uniform per, physical prad, uniform smet, 
    elif re.match('occur-per=(.*)-prad=physical-smet=(.*)', key):
        _, per, prad, smet = key.split('-')
        per_bins = _per_bins(per)
        prad_bins = prad_bins_dict['physical']
        smet_bins = _smet_bins(smet)
        occ = compute_occur_bins(per_bins, prad_bins, smet_bins)
        occ = set_index(occ,'size')

    # Occurence log uniform per, physical prad, uniform smet, 
    elif re.match('occur-per=(.*)-prad=twoperoctave', key):
        _, per, prad = key.split('-')
        per_bins = _per_bins(per)
        prad_bins = prad_bins_dict['two-per-octave']
        smet_bins = [-1,0.5]
        occ = compute_occur_bins(per_bins, prad_bins, smet_bins)
    else:
        assert False, "cannot parse {}".format(key)
    return occ

def load_fit(key):
    print key
    if key.count('fit_smet-')==1:
        _, dist, size = key.split('-')
        dsmet = 0.05 
        dper = 1
        dx = [dsmet]
        key = 'occur-per=hotwarm-prad=physical-smet={}'.format(dsmet)
        occ = cksmet.io.load_object(key,cache=1)

        cut = occ.df.ix[dist,size]
        cut = cut[cut.smetc.between(-0.4,0.4)]
        fit = cksmet.fit.Exponential(cut.smetc, dx, cut.nplnt, cut.ntrial)
        fit.fit()
        fit.mcmc(burn=300, steps=600, thin=1, nwalkers=100)

    elif key.count('fit_per-')==1:
        dper = 0.05 # Size of the bins used in the fitting
        dx = [dper]
        _, smet, size = key.split('-')
        occkey = 'occur-per={}-prad=physical-smet={}'.format(dper, smet)
        occ = cksmet.io.load_object(occkey,cache=1)
        cut = occ.df.ix[size,:]

        if (key.count('se')==1):
            cut = cut[cut.perc.between(1,100) & (cut.prob_det_mean > 0.25)]
        if (key.count('sn')==1):
            cut = cut[cut.perc.between(1,350) & (cut.prob_det_mean > 0.25)]
        if (key.count('ss')==1):
            cut = cut[cut.perc.between(1,350) & (cut.prob_det_mean > 0.25)]
        if (key.count('jup')==1):
            cut = cut[cut.perc.between(1,350) & (cut.prob_det_mean > 0.25)]

        fit = cksmet.fit.PowerLawCutoff(cut.perc, dx, cut.nplnt, cut.ntrial) 
        #fit.p0['binwper'].value = binwper

        if (key.count('sn')==1):
            fit.p0['logkp'].vary = True
            fit.p0['logkp'].value = -0.1

        fit.fit()
        fit.mcmc(burn=300, steps=600, thin=1, nwalkers=100)

    elif key.count('fit_persmet-')==1:
        _, dist, size = key.split('-')
        dper = 0.05 # Size of the bins used in the fitting
        dsmet = 0.05 # Size of the bins used in the fitting
        dx = dper
        
        occkey = 'occur-per={}-prad=physical-smet={}'.format(dper, dsmet)
        occ = cksmet.io.load_object(occkey)
        occ = set_index(occ,'size')
        cut = occ.df.ix[size]
        if key.count('hot')==1:
            cut = cut[cut.perc.between(1,10) & cut.smetc.between(-0.4,0.4)]
        elif key.count('warm'):
            cut = cut[cut.perc.between(10,100)  & cut.smetc.between(-0.4,0.4)]


        x = cut['perc smetc'.split()]
        nplnt = cut[['nplnt']]
        ntrial = cut[['ntrial']]
        fit = cksmet.fit.PerPowerLawExpSmetExp(x, dx, nplnt, ntrial)

        if (key.count('hot')==1):
            if (key.count('se')==1) or (key.count('sn')==1):
                fit.p0['alpha'].vary = True
                fit.p0['alpha'].value = 1.7
            elif (key.count('ss')==1) or (key.count('jup')==1):
                fit.p0['alpha'].vary = True
                fit.p0['alpha'].value = 0
                fit.p0['beta'].max = 10
            else:
                assert False, "error"
        if (key.count('warm')==1):
            if (key.count('se')==1):
                fit.p0['alpha'].vary = True
                fit.p0['alpha'].value = -0.6
            else:
                fit.p0['alpha'].vary = True
                fit.p0['alpha'].value = 0


        fit.fit()
        
        fit.mcmc(burn=300, steps=600, thin=1, nwalkers=100)

    else:
        assert False, "{} not supported".format(key)

    fit.corner()
    plt.gcf().savefig(key+'.pdf')
    fit.print_parameters()
    fit.occur = cut
    return fit


