from matplotlib.pylab import *
import pandas as pd
import seaborn as sns
import scipy
import xarray as xr

import ckscool.io
from .config import *

figsize=(3.5,2.5)
font_scale = 0.8

class ContourPlotter(object):
    def __init__(self, x, xerr, y, yerr, xscale='log', yscale='log'):
        self._kde_nx = 200
        self._kde_ny = 400
        #self._kde_nx = 10
        #self._kde_ny = 20 #debugging
        self.xscale = xscale
        self.yscale = yscale

        if xscale=='log':
            xerr = log10(1 + xerr / x)
            x = log10(x)

        if yscale=='log':
            yerr = log10(1 + yerr / y)
            y = log10(y)

        self.x = x
        self.y = y
        self.xerr = xerr
        self.yerr = yerr
        self.n = len(self.x)

    def compute_density(self, for_gradient=False):
        ndim = 2
        kx = np.linspace(self.xmin, self.xmax, self._kde_nx)
        ky = np.linspace(self.ymin, self.ymax, self._kde_ny)
        if self.xscale=='log':
            kx = np.linspace(log10(self.xmin), log10(self.xmax), self._kde_nx)
        if self.yscale=='log':
            ky = np.linspace(log10(self.ymin), log10(self.ymax), self._kde_ny)

        coords = {'kx':kx, 'ky':ky}
        _kx, _ky = meshgrid(kx,ky,indexing='ij')
        data = {
            'kxc': (['kx', 'ky'], _kx),
            'kyc': (['kx', 'ky'], _ky),
        }
        ds = xr.Dataset(data,coords=coords)

        # midpoints
        mu = np.zeros((self.n,ndim))
        mu[:,0] = self.x
        mu[:,1] = self.y

        # construct covarience matrix
        cov = zeros((self.n,ndim,ndim))
        cov[:,0,0] = self.xerr**2
        cov[:,1,1] = self.yerr**2

        pos = np.empty(_kx.shape + (2,))
        pos[:, :, 0] = _kx
        pos[:, :, 1] = _ky
        Z = gaussian(pos, mu, cov)
        ds['Z'] = (['kx','ky'],Z)
        ds['Z'] /= ds['Z'].max() # rescale
        self.ds = ds

        if for_gradient:
            return kx, ky, Z

    def normalize_density(self):
        frac = self.ds['Z'].sum(dim='ky')
        frac = frac / frac.mean()
        self.ds['Z'] /= frac

    def plot_contour(self):
        _vmax = self.ds.Z.max()
        fac = 1.2
        Z = self.ds.Z / fac

        qc = Z.plot.contourf(
            x='kx', cmap=plt.cm.afmhot_r, levels=arange(0,1.001,0.1),zorder=0,vmax=1,
            add_colorbar=False
        )
        #add_anchored('$N_p$ = {}'.format(len(self.x)),2,prop=dict(size='small'))
        self.xlim(self.xmin,self.xmax)
        self.ylim(self.ymin,self.ymax)
        return qc
        

    def xlim(self, *args):
        if self.xscale=='log':
            args2 = (log10(args[0]),log10(args[1]))
        else:
            args2 = args
        xlim(*args2)

    def ylim(self, *args):
        if self.yscale=='log':
            args2 = (log10(args[0]),log10(args[1]))
        else:
            args2 = args
        ylim(*args2)

    def xticks(self, ticks):
        if self.xscale=='log':
            args2 = (log10(ticks),ticks)
        else:
            args2 = (ticks,ticks)
        xticks(*args2)


    def yticks(self, ticks):
        if self.yscale=='log':
            args2 = (log10(ticks),ticks)
        else:
            args2 = (ticks,ticks)
        yticks(*args2)

    def plot_integrated(self):
        pass

    def plot_points(self):
        #errorbar(self.x,self.y,yerr=self.yerr,fmt='o',mfc='none',elinewidth=0,ms=2,mew=0.5,errorevery=0)
        errorbar(self.x,self.y,yerr=None,fmt='o',mfc='none',elinewidth=0,ms=2,mew=0.4,mec='w',zorder=9)
        errorbar(self.x,self.y,yerr=None,fmt='o',mfc='none',elinewidth=0,ms=2,mew=0.2,mec='k',zorder=10)

class Plotter(object):
    def __init__(self,xk,nopoints=False,zoom=False,normalize=False,query=None,
                   yerrfac=1,xerrfac=1,for_gradient=False):

        df = ckscool.io.load_table('planets-cuts2',cache=1)
        df = df[~df.isany]
        yk = 'gdir_prad'
        yerrk = 'gdir_prad_err1'
        y = df[yk]
        x = df[xk]
        yerr = df[yerrk] * yerrfac

        if xk=='koi_period':
            xerr = 0.5 * x # plottting purposes
            p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
            if zoom:
                p1.xmin = 1
                p1.xmax = 100
                p1.ymin = 1
                p1.ymax = 4
            else:
                p1.xmin = 0.3
                p1.xmax = 300
                p1.ymin = 0.5
                p1.ymax = 16

            xlim = (p1.xmin,p1.xmax)
            ylim = (p1.ymin,p1.ymax)
            xticks = [0.3,1,3,10,30,100,300]
            xlabel = 'Orbital Period (days)'

        if xk=='giso_sinc':
            xerr = x # plottting purposes
            p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
            if zoom:
                p1.xmin = 3e0
                p1.xmax = 3e3
                p1.ymin = 1.0
                p1.ymax = 4
            else:
                p1.xmin = 1
                p1.xmax = 1e4
                p1.ymin = 0.5
                p1.ymax = 16

            xlim = (p1.xmax,p1.xmin)
            ylim = (p1.ymin,p1.ymax)
            xticks = [1,3,10,30,100,300,1000,3000,10000]
            xlabel = 'Incident Bolometric Flux (Earth-units)'

        if xk=='giso_smass':
            xerr = 0.2 * x # plottting purposes
            p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
            if zoom:
                p1.xmin = 0.5
                p1.xmax = 1.4
                p1.ymin = 1
                p1.ymax = 4
            else:
                p1.xmin = 0.5
                p1.xmax = 1.5
                p1.ymin = 0.5
                p1.ymax = 16

            xlim = (p1.xmin,p1.xmax)
            ylim = (p1.ymin,p1.ymax)
            xticks = [0.5,0.7,1.0,1.4]
            xlabel = 'Stellar Mass (Solar-masses)'

        if xk=='cks_smet':
            xerr = 0.15
            p1 = ContourPlotter(x, xerr, y, yerr,xscale='lin',yscale='log')
            if zoom:
                p1.xmin = -0.4
                p1.xmax = 0.4
                p1.ymin = 1
                p1.ymax = 4
            else:
                p1.xmin = -0.5
                p1.xmax = 0.5
                p1.ymin = 0.5
                p1.ymax = 16

            xlim = (p1.xmin,p1.xmax)
            ylim = (p1.ymin,p1.ymax)
            xticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
            xlabel = '[Fe/H]'

        if xk=='giso_slogage':
            xerr = 0.15
            p1 = ContourPlotter(x, xerr, y, yerr,xscale='lin',yscale='log')
            if zoom:
                p1.xmin = 8.8
                p1.xmax = 10.2
                p1.ymin = 1
                p1.ymax = 4
            else:
                p1.xmin = 8.5
                p1.xmax = 10.5
                p1.ymin = 0.5
                p1.ymax = 16

            xlim = (p1.xmin,p1.xmax)
            ylim = (p1.ymin,p1.ymax)
            xticks = [8.0,8.5,9.0,9.5,10.0,10.5]
            xlabel = 'logage'

        yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]
        ylabel = 'Planet Size (Earth-radii)'

        self.p1 = p1
        self.xticks = xticks
        self.yticks = yticks
        self.xlim = xlim
        self.ylim = ylim
        self.ylabel = ylabel
        self.xlabel = xlabel
        if for_gradient:
            X, Y, Z = p1.compute_density(for_gradient=True)
            return X, Y, Z

    def setp(self):
        self.p1.xticks(self.xticks)
        self.p1.yticks(self.yticks)
        self.p1.xlim(*self.xlim)
        self.p1.ylim(*self.ylim)
        setp(gca(),xlabel=self.xlabel,ylabel=self.ylabel)
        
def fig_sample(**kwargs):
    sns.set_context('paper',font_scale=1.0)
    fig,axL = subplots(ncols=2,nrows=2,figsize=(7,6))

    # Period, radius
    sca(axL[0,0])
    pl = Plotter('koi_period',**kwargs)
    pl.p1.compute_density()
    pl.p1.plot_points()
    qm = pl.p1.plot_contour()
    pl.setp()
    fig_label('a')
    ax = axes([0.9, 0.85, 0.008, 0.1])
    cbar = colorbar(qm,cax=ax,format='%.1f')
    cbar.set_label('relative density',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Sinc, radius
    sca(axL[0,1])
    pl = Plotter('giso_sinc',**kwargs)
    pl.p1.compute_density()
    pl.p1.plot_points()
    pl.p1.plot_contour()
    pl.setp()
    fig_label('b')

    # Metallicity and age
    sca(axL[1,0])
    pl = Plotter('giso_smass',**kwargs)
    pl.p1.compute_density()
    pl.p1.plot_points()
    pl.p1.plot_contour()
    pl.setp()
    fig_label('d')
    
    # Mass and age
    sca(axL[1,1])
    pl = Plotter('cks_smet',**kwargs)
    pl.p1.compute_density()
    pl.p1.plot_points()
    pl.p1.plot_contour()
    pl.setp()
    fig_label('d')

    tight_layout(True)

def gaussian(pos, mu, cov):

    """
    Compute the sum of 2D gaussians

    last dimension of x must equal len(mu)
    """
    assert mu.shape[0]==cov.shape[0]
    assert mu.shape[1]==cov.shape[1]

    out = []
    for i in range(mu.shape[0]):
        _out = scipy.stats.multivariate_normal(mu[i], cov[i]).pdf(pos)
        out.append(_out)

    out = np.array(out)
    out = np.sum(out,axis=0)
    return out
