from matplotlib.pylab import *
import pandas as pd
import seaborn as sns
import scipy
import xarray as xr

import ckscool.io
from .config import *
from .contour import ContourPlotter
from ckscool.occur import gaussian_2d_kde
figsize=(3.5,2.5)
font_scale = 0.8

class Plotter(object):
    def __init__(self, df, xk, nopoints=False, zoom=False,normalize=False,query=None,
                   yerrfac=1, xerrfac=1, for_gradient=False):

        yk = 'gdir_prad'
        self.x = df[xk]
        self.y = df[yk]

        if xk=='koi_period':
            if zoom:
                cp = ContourPlotter(1,100,1,4, xscale='log',yscale='log')
            else:
                cp = ContourPlotter(0.3,300,0.5,16, xscale='log',yscale='log')
            xticks = [0.3,1,3,10,30,100,300]
            xlabel = 'Orbital Period (days)'
            self.bwx = log10(1 + 1)
            self.bwy = log10(1+0.05)
            
        if xk=='giso_sinc':
            if zoom:
                cp = ContourPlotter(3e0, 3e3, 1, 4,xscale='log',yscale='log')
            else:
                cp = ContourPlotter(1, 1e4, 0.5, 16,xscale='log',yscale='log')

            xticks = [1,3,10,30,100,300,1000,3000,10000]
            xlabel = 'Incident Bolometric Flux (Earth-units)'
            self.bwx = log10(1 + 1)
            self.bwy = log10(1+0.05)
            
        if xk=='giso_smass':
            if zoom:
                cp = ContourPlotter(0.5, 1.4, 1, 4,xscale='log',yscale='log')
            else:
                cp = ContourPlotter(0.5, 1.5, 0.5, 16,xscale='log',yscale='log')

            xticks = [0.5,0.7,1.0,1.4]
            xlabel = 'Stellar Mass (Solar-masses)'
            self.bwx = log10(1 + 0.15)
            self.bwy = log10(1+0.05)

        if xk=='cks_smet':
            xerr = 0.15
            if zoom:
                cp = ContourPlotter(-0.4, 0.4, 1, 4,xscale='lin',yscale='log')
            else:
                cp = ContourPlotter(-0.5, 0.5, 0.5,16,xscale='log',yscale='log')
            xticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
            xlabel = '[Fe/H]'
            self.bwx = 0.1
            self.bwy = log10(1+0.05)

        yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]
        ylabel = 'Planet Size (Earth-radii)'
        self.kx = self.x
        self.ky = self.y
        if cp.xscale=='log':
            self.kx = log10(self.x)
        if cp.yscale=='log':
            self.ky = log10(self.y)

        self.cp = cp
        self.xticks = xticks
        self.yticks = yticks
        self.ylabel = ylabel
        self.xlabel = xlabel
        if for_gradient:
            X, Y, Z = cp.compute_density(for_gradient=True)
            return X, Y, Z

    def plot(self):
        ds = self.cp.meshgrid()
        Z = gaussian_2d_kde(array(ds.kxc), array(ds.kyc), self.kx,
                            self.ky, self.bwx, self.bwy)
        Z = Z.reshape(ds.kxc.shape)
        ds['Z'] = (['kx','ky'],Z)
        qc = self.cp.contour(ds['Z'], normalize=True, cmap=plt.cm.afmhot_r)
        self.qc = qc
        self.cp.errorbar(self.x, self.y, yerr=None, fmt='o',
                         mfc='none', elinewidth=0, ms=2, mew=0.4,
                         mec='w', zorder=9)

        self.cp.errorbar(self.x, self.y, yerr=None, fmt='o',
                         mfc='none', elinewidth=0, ms=2, mew=0.2,
                         mec='k', zorder=10)

        self.cp.yticks(self.yticks)
        self.cp.xticks(self.xticks)
        self.cp.yticks(self.yticks)
        self.cp.set_lim()
        setp(gca(),xlabel=self.xlabel,ylabel=self.ylabel)
    
def fig_sample(**kwargs):
    sns.set_context('paper',font_scale=1.0)
    fig,axL = subplots(ncols=2,nrows=2,figsize=(7,6))

    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]

    # Period
    sca(axL[0,0])
    pl = Plotter(df,'koi_period',**kwargs)
    pl.plot()
    fig_label('a')
    ax = axes([0.89, 0.85, 0.008, 0.1])
    cbar = colorbar(pl.qc,cax=ax,format='%.1f')
    cbar.set_label('relative density',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Sinc
    sca(axL[0,1])
    pl = Plotter(df,'giso_sinc',**kwargs)
    pl.plot()
    xl = xlim()
    xlim(xl[1],xl[0])
    fig_label('b')

    # Mass 
    sca(axL[1,0])
    pl = Plotter(df,'giso_smass',**kwargs)
    pl.plot()
    fig_label('c')
    
    # Metallicity
    sca(axL[1,1])
    pl = Plotter(df,'cks_smet',**kwargs)
    pl.plot()
    fig_label('d')
    tight_layout(True)
