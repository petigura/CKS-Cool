import pandas as pd
from matplotlib.transforms import blended_transform_factory as btf
from matplotlib.pylab import *
import specmatchemp.specmatch
import specmatchemp.library

import ckscool.io

def get_spectra():


    df = ckscool.io.load_table('planets-cuts2+iso',cache=2)
    rea = pd.read_csv('data/reamatch.csv',usecols=[0,1,2,3,4,5])
    df = pd.merge(df, rea[['id_koi','is_sb2']],how='left')

    cut = df[~df.is_sb2.isin([4,5])]
    bins = arange(3700,5100,100)
    cut = cut.groupby(pd.cut(cut.cks_steff,bins)).nth(0)
    basedir = "/data/user/petigura/public_html/smemp/results/"
    cut['file'] = cut.apply(lambda x : basedir + "{id_name:}_{id_obs:}/{id_name:}_{id_obs:}_sm.hdf".format(**x),axis=1)
    return cut 
    
def fig_spectra():
    import seaborn as sns
    sns.set_context('paper',font_scale=1.3)
    sns.set_style('ticks')

    cut = get_spectra()
    fig = plt.figure(figsize=(8.5,8))
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3,rowspan=2)
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3)
    axL = [ax1,ax2]
    trans = btf(ax1.transAxes,ax1.transData)
    lib = specmatchemp.library.read_hdf()

    cut = cut.sort_values(by='cks_steff')
    lw = 0.7
    shift = 0 
    for i, row in cut.iloc[::2].iterrows():
        fn = 'data/fig_spectra'+row.file
        sm = specmatchemp.specmatch.SpecMatch.read_hdf(fn,lib)
        mt = sm.lincomb_matches[1]
        sca(ax1)
        plot(mt.w,mt.target.s+ shift,'-', color='black',lw=lw)
        plot(mt.w,mt.modified.s +shift, '-', color='red',lw=lw)
        s = "{} {}".format(row.id_name,row.sm_steff)
        text(0.03,shift + 0.8,s,transform=trans)

        sca(ax2)
        resid = mt.modified.s-mt.target.s
        plot(mt.w, resid + shift * 0.5, '-', color='red',lw=lw)
        shift+=0.5
        print "reading in {}".format(fn)
        print "{}, {}".format(mean(mt.target.serr),np.std(resid))

    setp(axL,xlim=(5160,5190))
    #setp([ax1,ax2],xlim=(5160,5190))

    fig.set_tight_layout(True)

    for ax in axL:
        ax.grid()

    setp(ax1,ylabel='Normalized Intensity + Offset',xlabel='Wavelength (Ang)')
    setp(ax2,ylabel='Residuals + Offset',xlabel='Wavelength (Ang)')



def fig_spectra_wide():
    import seaborn as sns
    sns.set_context('paper',font_scale=1.3)
    sns.set_style('ticks')

    df = ckscool.io.load_table('planets-cuts2+iso',cache=2)
    rea = pd.read_csv('data/reamatch.csv',usecols=[0,1,2,3,4,5])
    df = pd.merge(df, rea[['id_koi','is_sb2']],how='left')

    cut = df[~df.is_sb2.isin([4,5])]
    bins = arange(4000,5000,250)
    cut = cut.groupby(pd.cut(cut.cks_steff,bins)).nth(0)
    basedir = "/data/user/petigura/public_html/smemp/results/"

    fig,ax1 = plt.subplots(figsize=(8.5,3))
    trans = btf(ax1.transAxes,ax1.transData)
    lib = specmatchemp.library.read_hdf()

    lw = 0.7
    shift = 0 

    files = ['K00854_rj158.280','CK00747_rj296.91','CK01101_rj293.491']
    steffs = ['3695','4287','4892',]

    for fn, steff in zip(files,steffs):

        _fn = 'data/fig_spectra/data/user/petigura/public_html/smemp/results/{}/{}_sm.hdf'.format(fn,fn)
        sm = specmatchemp.specmatch.SpecMatch.read_hdf(_fn,lib)
        mt = sm.lincomb_matches[1]
        sca(ax1)
        plot(mt.w,mt.target.s+ shift,'-', color='black',lw=lw)
        plot(mt.w,mt.modified.s +shift, '-', color='red',lw=lw)
        name = fn.split('_')[0]
        name = "KOI-{}".format(int(name[-5:]))
        s = r"{} Teff = {}".format(name,steff)
        text(0.01,shift + 1.1,s,transform=trans,size='small')
        shift+=1.1
    setp(ax1,xlim=(5160,5190),ylim=(0,4))

    #setp([ax1,ax2],xlim=(5160,5190))

    fig.set_tight_layout(True)


    setp(ax1,ylabel='Normalized Intensity + Offset',xlabel='Wavelength (Ang)')
    fig.savefig('fig_spectra-wide.pdf')
