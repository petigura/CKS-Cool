from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
import ckscool.io

def fig_planets_warm_smet_smass():
    sns.set_context('paper')
    fig, axL = subplots(ncols=2, figsize=(4.5,6),sharey=True)
    df = ckscool.io.load_table('planets-cuts2+iso')
    df = df[~df.isany]
    df = df.query('10 < koi_period < 100')
    query = '-0.2 < cks_smet < 0.2 and 0.6 < giso_smass < 1.2'
    df = df.query(query)

    sca(axL[0])
    semilogy()
    xk = 'cks_smet'
    yk = 'gdir_prad'
    bins = [-0.2,-0.1,0.0,0.1,0.2]
    plot_with_bins(df,xk,yk,bins)
    xlim(-0.3,0.3) 
    xlabel('Metallicity (dex)')
    ylabel('Planet size (Earth-radii)')
    tight_layout()

    sca(axL[1])
    #loglog()
    xk = 'giso_smass'
    yk = 'gdir_prad'
    bins = [0.6,0.7,0.8,0.9,1.0,1.1,1.2]
    #bins = logspace(log10(0.6),log10(1.2),7)
    plot_with_bins(df,xk,yk,bins)
    xlabel('Mass')
    tight_layout()

def plot_with_bins(df,xk,yk,bins):
    bins = np.array(bins)
    rbins = [1.0, 1.7, 4.0]
    plot(df[xk], df[yk],'.')
    for i in range(2):
        rlo = rbins[i]
        rhi = rbins[i+1]
        cut = df[df.gdir_prad.between(rlo,rhi)]
        g = cut.groupby(pd.cut(cut[xk],bins=bins))
        
        x = vstack([bins[:-1],bins[1:]])
        xbin = x.mean(axis=0)
        ybin = np.array(g[yk].mean())
        yerrbin = np.array(g[yk].std()/sqrt(g[yk].count()))

        # Median statistics
        #ybin = np.array(g[yk].median())
        #yerrbin = np.array(1.5 * g[yk].std()/sqrt(g[yk].count()))
        y = array(ybin)
        y = array(2*[y])
        errorbar(xbin,ybin,yerr=yerrbin,fmt='or',ms=1,lw=2,zorder=5)
        plot(x,y,'r',lw=2,zorder=5)

    for be in rbins:
        axhline(be,ls='--',color='r')

    minorticks_off()
    yt = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]
    yticks(yt,yt)
    ylim(0.7,16)

def draw_fiducial_line():
    plot([-0.5,0.5],[2.15,2.55])
    bins = [-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]
    xk = 'cks_smet'
    yk = 'gdir_prad'

