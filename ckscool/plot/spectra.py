import pandas as pd
from matplotlib.transforms import blended_transform_factory as btf
from matplotlib.pylab import *
import specmatchemp.specmatch
import specmatchemp.library
import matplotlib.patheffects as PathEffects
import ckscool.io

def fig_spectra_wide():
    import seaborn as sns
    sns.set_context('paper',font_scale=1.1)
    sns.set_style('ticks')

    sns.set_color_codes()
    fig,ax1 = plt.subplots(figsize=(8.5,3.5))
    trans = btf(ax1.transAxes,ax1.transData)
    lib = specmatchemp.library.read_hdf()

    step = 0.5
    lw = 0.7
    shift = 0

    files = ['K00854_rj158.280',
             'K01078_rj159.690',
             'CK00739_rj300.277',
             'CK00747_rj296.91',
             'K00784_rj274.108',
             'CK00740_rj292.812',
             'CK01101_rj293.491']
    steffs = ['3695','3841','4093','4287','4423','4687','4892']

    for fn, steff in zip(files,steffs):

        _fn = 'fig_spectra/data/user/petigura/public_html/smemp/results/{}/{}_sm.hdf'.format(fn,fn)
        sm = specmatchemp.specmatch.SpecMatch.read_hdf(_fn,lib)
        mt = sm.lincomb_matches[1]
        sca(ax1)
        plot(mt.w,mt.target.s+ shift,'-', color='black',lw=lw)
        plot(mt.w,mt.modified.s +shift, '-', color='red',lw=lw)
        name = fn.split('_')[0]
        name = "KOI-{}".format(int(name[-5:]))
        #s = r"{} $T_\mathregular{{eff}}$ = {}".format(name,steff)
        s = r"{}, {} (K)".format(name,steff)
        txt = text(0.01,shift +0.8,s,transform=trans,size='small')

        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        shift+=step

    #setp([ax1,ax2],xlim=(5160,5190))

    fig.set_tight_layout(True)
    setp(ax1,ylabel='Normalized Intensity + Offset',xlabel='Wavelength (Ang)')
    ax1.grid()
    yt = [0,1,2,3,4,5]
    yticks(yt,yt)
    setp(ax1,xlim=(5160,5190),ylim=(0,4))
