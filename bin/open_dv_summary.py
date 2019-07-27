#!/usr/bin/env python
from argparse import ArgumentParser

def main():
    psr = ArgumentParser()
    psr2 = subpsr.add_argument('create-iso-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_batch)


    args = psr.parse_args()
    args.func(args)

def create_iso_batch(args):
    ckscool._isoclassify.create_iso_batch(args)

def create_iso_table(args):
    ckscool._isoclassify.create_iso_table(args)

def create_spectra_files(args):
    cut = ckscool.plot.spectra.get_spectra()
    fn = 'data/fig_spectra/fig_spectra-files.txt'
    cut.file.to_csv(fn,header=None,index=None) 
    print "created {}".format(fn)

def create_csv(args):
    df = ckscool.io.load_table(args.name)
    fn = 'data/{}.csv'.format(args.name)
    df.to_csv(fn)
    print "created {}".format(fn)

def create_chain_hdf(args):
    """
    Loop over every candidate in the thompson table and create a chain
    """

    t18 = keprat.io.load_table('t18')
    for i, row in t18.iterrows():
        table = 'chains-dr25-{}'.format(row.id_koicand)
        try:
            cachefn='data/dr25-mcmc-chains.hdf'
            df = keprat.io.load_table(table,cache=1,cachefn=cachefn)
        except IOError:
            with file('data/dr25-mcmc-chains_missing.txt','a') as f:
                s = table+" missing\n"
                print s
                f.write(s)

def create_chain_plots(args):
    import keprat
    from chainconsumer import ChainConsumer
    fmt = {}
    fmt['id_koicand'] = args.id_koicand
    table = "chains-dr25-{id_koicand:}".format(**fmt)
    cachefn='data/dr25-mcmc-chains.hdf'
    dr25 = keprat.io.load_table(table,cachefn=cachefn,cache=1)
    cols = "RHO ZPT EP1 PE1 BB1 RD1".split()
    parameters = [r"$\rho_{\star,\mathrm{circ}}$",r"$F_0$",r"$T_0$", r"$P$",r"$b$",r"$R_P / R_\star$"]
    c = ChainConsumer()
    chain = np.array(dr25[cols].iloc[::10])
    extents = zip(chain.min(axis=0),chain.max(axis=0))
    c.add_chain(chain, parameters=parameters)
    fn = 'mcmc/{id_koicand:}_distributions.png'.format(**fmt)
    c.plotter.plot_walks(parameters=parameters,extents=extents,filename=fn)
    print "created {}".format(fn)

def create_chain_plots_batch(args):
    cachefn='data/dr25-mcmc-chains.hdf'
    t18 = keprat.io.load_table('t18')
    t18.index = t18.id_koicand
    df = pd.read_csv('data/dr25-mcmc-chains_missing.txt',sep=' ',names=['s','x'])
    df['id_koicand'] = df.s.apply(lambda x : x.split('-')[-1]).str.upper()
    t18 = t18.drop(df.id_koicand)
    for i, row in t18.iterrows():
        print "run_ckscool.py create-chain-plots {id_koicand:}".format(**row)
