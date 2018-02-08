# Code from EAP's K2-24 paper use as a template
'''
def val_stat(return_dict=False):
    d = OrderedDict()
    # load up data from p16
    p16 = ktwo24.io.load_table('petigura16')
    d['p16-ror1'] = p16.ix['pl_ror',1]
    d['p16-ror1_err'] = p16.ix['pl_ror_err',1]
    d['p16-ror2'] = p16.ix['pl_ror',2]
    d['p16-ror2_err'] = p16.ix['pl_ror_err',2]

    vst, meta = cpsutils.io.read_vst('data/rv/vstepic203771098.dat',full_output=True, verbose=0)
    d['rv-n'] = len(vst)
    d['rv-start'] = vst.iloc[0]['day']
    d['rv-stop'] = vst.iloc[-1]['day']
    d['rv-errvel-min'] = "{:.1f}".format(vst.errvel.min())
    d['rv-errvel-max'] = "{:.1f}".format(vst.errvel.max())

    ephem = ktwo24.io.load_table('ephem-lithwick')
    d['ref-ephem-per-1'] = ephem.ix[1,'per']
    d['ref-ephem-tc-1'] = ephem.ix[1,'T']
    d['ref-ephem-per-2'] = ephem.ix[2,'per']
    d['ref-ephem-tc-2'] = ephem.ix[2,'T']

    inpfn = 'analysis/ttv-lithwick/emcee-lithwick-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['lithwick-prior-mu1'] = "{:.0f}".format(mod.mu1*1e6)
    d['lithwick-prior-mu1_err'] = "{:.0f}".format(mod.mu1_err*1e6)
    d['lithwick-prior-mu2'] = "{:.0f}".format(mod.mu2*1e6)
    d['lithwick-prior-mu2_err'] = "{:.0f}".format(mod.mu2_err*1e6)

    inpfn = 'analysis/ttv-ttvfast/mnest-ttvfast-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['ttvfast-tstart'] =  mod.tstart 

    inpfn = 'analysis/ttv-ttvfast/emcee-ttvfast-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['ttvfast-emcee-nwalkers'] =  mod.nwalkers
    d['ttvfast-emcee-nsteps'] =  mod.nsteps
    d['ttvfast-emcee-nburn'] =  mod.nburn

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines
'''
