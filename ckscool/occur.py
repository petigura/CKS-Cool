"""
Defines occurrence object
"""

import cksmet.stats
from statsmodels.distributions.empirical_distribution import ECDF
import ckscool.io
import numpy as np
from numpy import log10, logspace, arange, round

class Occurrence(object):
    def __init__(self, plnt, comp, nstars, smet_field=None):
        self.plnt = plnt
        self.comp = comp
        self.nstars = nstars
        self.plnt = plnt

    def occurence_box(self, limits):
        """Compute occurrence in a little box

        We make the assumption that the dN/dlogP and dN/lopRp is constant
        within a box.

        Args:
            limits (dict): must contain, per1, per2, prad1, prad2, can
                optionally contain, smet1, smet2
        """
        out = dict()
        prad1 = limits['prad1']
        prad2 = limits['prad2']
        per1 = limits['per1']
        per2 = limits['per2']

        # Get planet sample and number of stars
        cut = self.plnt.copy()
        cut = cut[cut.prad.between(prad1,prad2)]
        cut = cut[cut.per.between(per1,per2)]
        nstars = self.nstars
        nplnt = len(cut)
        prob_trdet_mean, prob_det_mean = self.comp.mean_prob_trdet(
            per1, per2, prad1, prad2
        )
        ntrial = nstars * prob_trdet_mean
        rate = nplnt / ntrial

        nsample = int(1e4)
        binom = cksmet.stats.Binomial(ntrial, nplnt)
        samples = binom.sample(nsample) 

        uplim = nplnt==0
        rate = cksmet.stats.samples_to_rate(samples,uplim=uplim)
        out['ntrial'] = ntrial
        out['nplnt'] = nplnt
        out['prob_trdet_mean'] = prob_trdet_mean
        out['prob_det_mean'] = prob_det_mean
        out = dict(out,**rate)
        return out

def load_occur(limits):
    """
    Constructs occurrence object
    """
    
    # Derive completeness object
    method = 'fulton-gamma-clip' # treatment for planet detectability
    impact = 0.7 # maximum impact parameter considered.

    field = ckscool.io.load_table('field-cuts',cache=1)
    field = field[~field.isany]
    field = field.rename(columns={'ber18_srad':'srad','m17_smass':'smass'})

    plnt = ckscool.io.load_table('planets-cuts2+iso')
    plnt = plnt[~plnt.isany]
    namemap = {'gdir_prad':'prad','koi_period':'per','giso_smass':'smass'}
    plnt = plnt.rename(columns=namemap)

    if limits.has_key('smass1'):
        smass1 = limits['smass1']  
        smass2 = limits['smass2']
        field = field[field.smass.between(smass1,smass2)]
        plnt = plnt[plnt.smass.between(smass1,smass2)]

    elif limits.has_key('bmr1'):
        bmr1 = limits['bmr1']  
        bmr2 = limits['bmr2']
        xs = 'gaia2_bpmag - gaia2_rpmag'
        field = field[field.eval(xs).between(bmr1,bmr2)]
        plnt = plnt[plnt.eval(xs).between(bmr1,bmr2)]



    n1 = len(field)
    field = field.dropna(subset=ckscool.comp.__STARS_REQUIRED_COLUMNS__)
    n2 = len(field)
    print "{}/{} stars remain after droping nulls ".format(n2,n1)

    # Define grid of period and radius to compute completeness
    comp_per_bins = np.round(logspace(log10(0.1),log10(1000),65),4)
    comp_prad_bins = np.round(logspace(log10(0.25),log10(64),51 ),2)
    
    # debugging
    if False:
        comp_per_bins = comp_per_bins[:6]
        comp_prad_bins = comp_prad_bins[:6]
        
    
    comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = ckscool.grid.Grid(comp_bins_dict,spacing_dict)

    comp = ckscool.comp.Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det(verbose=True)
    comp.compute_grid_prob_tr(verbose=True)
    comp.create_splines()
    
    nstars = len(field)
    occ = ckscool.occur.Occurrence(plnt, comp, nstars)
    return occ


