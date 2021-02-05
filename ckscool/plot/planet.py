from matplotlib.pylab import *
import seaborn as sns
import xarray as xr

import ckscool.io
from ckscool.occur import gaussian_2d_kde
import ckscool.gradient
from .config import *
    
def fig_sample(plot_gradient=False, **kwargs):
    sns.set_context('paper',font_scale=1.0)
    fig,axL = subplots(ncols=2,nrows=2,figsize=(7,6))

    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]

    # Period
    sca(axL[0,0])
    pl = NDPlotter(df,'koi_period',smass_lims=[0.5,1.4],**kwargs)
    pl.plot()
    if plot_gradient:
        key = 'grad-per-prad-det_smass=0.5-1.4'
        grads = ckscool.io.load_object(key, cache=1) 
        pl.cp.plot_gradients(grads)
        
    fig_label('a')
    ax = axes([0.89, 0.85, 0.008, 0.1])
    cbar = colorbar(pl.qc,cax=ax,format='%.1f')
    cbar.set_label('relative density',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Sinc
    sca(axL[0,1])
    pl = NDPlotter(df,'giso_sinc',smass_lims=[0.5,1.4],**kwargs)
    pl.plot()
    if plot_gradient:
        key = 'grad-sinc-prad-det_smass=0.5-1.4'
        grads = ckscool.io.load_object(key, cache=1) 
        pl.cp.plot_gradients(grads)

    xl = xlim()
    xlim(xl[1],xl[0])
    fig_label('b')

    # Mass 
    sca(axL[1,0])
    pl = NDPlotter(df,'giso_smass',**kwargs)
    pl.plot()
    fig_label('c')
    
    # Metallicity
    sca(axL[1,1])
    pl = NDPlotter(df,'cks_smet',**kwargs)
    pl.plot()
    fig_label('d')
    tight_layout(True)


def fig_gupta_comparison(plot_gradient=False, **kwargs):
    sns.set_context('paper',font_scale=1.0)
    fig,axL = subplots(ncols=3,nrows=2,figsize=(10,6))

    kwline = dict(linestyle='--',color='white')
    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]

    # Mass
    pl = NDPlotter(df,'giso_smass',zoom=True,**kwargs)
    xi = linspace(log10(0.5),log10(1.4))
    yi = 0.33 * xi + log10(1.6)

    sca(axL[0,0])
    pl.plot()
    plot(xi,yi,**kwline)

    sca(axL[1,0])
    pl.plot(column_normalize=True)
    plot(xi,yi,**kwline)
    
    # Metalicity
    pl = NDPlotter(df,'cks_smet',zoom=True)
    xi = linspace(-0.4,0.4)
    yi = 0.1 * xi + log10(2.3)

    sca(axL[0,1])
    pl.plot()
    plot(xi,yi,**kwline)

    sca(axL[1,1])
    pl.cp.xmin=-0.3
    pl.cp.xmax=0.3
    pl.plot(column_normalize=True)
    pl.cp.xlim(-0.4,0.4)
    plot(xi,yi,**kwline)
    
    # Age
    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]
    df = df.query('cks_steff > 5500')
    pl = NDPlotter(df,'giso_sage',zoom=True)

    xi = linspace(log10(1),log10(10))
    yi = -0.1 * xi + log10(2.5)
    sca(axL[0,2])
    pl.plot()
    pl.cp.xlim(1,10)
    plot(xi,yi,**kwline)

    sca(axL[1,2])

    pl.cp.xmin=1.2
    pl.cp.xmax=10
    pl.plot(column_normalize=True)
    pl.cp.xlim(1,10)
    plot(xi,yi,**kwline)
    tight_layout(True)

'''    
    
    

plot(df.giso_slogage,df.gdir_prad,'.')

# Period
#sca(axL[0,0])


pl = NDPlotter(df,'giso_sage',zoom=True)



df = ckscool.io.load_table('planets-cuts2',cache=1)
df = df[~df.isany]
#df = df.query('-0.2 < cks_smet < 0.2 and 0.5 < giso_smass < 1.0')
pl = NDPlotter(df,'cks_smet',zoom=True)
pl.cp.xmin=-0.35
pl.cp.xmax=0.35
pl.plot(column_normalize=True)
pl.cp.xlim(-0.4,0.4)



figure()

pl.plot(column_normalize=False)



    
    # Period
    sca(axL[0,0])
    pl = NDPlotter(df,'koi_period',smass_lims=[0.5,1.4],**kwargs)
    pl.plot()
    if plot_gradient:
        key = 'grad-per-prad-det_smass=0.5-1.4'
        grads = ckscool.io.load_object(key, cache=1) 
        pl.cp.plot_gradients(grads)
        
    fig_label('a')
    ax = axes([0.89, 0.85, 0.008, 0.1])
    cbar = colorbar(pl.qc,cax=ax,format='%.1f')
    cbar.set_label('relative density',size='x-small')
    cbar.ax.tick_params(labelsize='xx-small')

    # Sinc
    sca(axL[0,1])
    pl = NDPlotter(df,'giso_sinc',smass_lims=[0.5,1.4],**kwargs)
    pl.plot()
    if plot_gradient:
        key = 'grad-sinc-prad-det_smass=0.5-1.4'
        grads = ckscool.io.load_object(key, cache=1) 
        pl.cp.plot_gradients(grads)

    xl = xlim()
    xlim(xl[1],xl[0])
    fig_label('b')

    # Mass 
    sca(axL[1,0])
    pl = NDPlotter(df,'giso_smass',**kwargs)
    pl.plot()
    fig_label('c')
    
    # Metallicity
    sca(axL[1,1])
    pl = NDPlotter(df,'cks_smet',**kwargs)
    pl.plot()
    fig_label('d')
    tight_layout(True)
'''
    
class NDPlotter(object):
    """
    Class to facillitate plotting of occurrence number density
    """
    def __init__(self, df, xk, smass_lims=False, zoom=False):
        yk = 'gdir_prad'
        self.x = df[xk]
        self.y = df[yk]
        self.xk = xk
        self.smass_lims = smass_lims
        self.bwy = log10(1+0.07)

        if xk=='koi_period':
            if zoom:
                xmin = 1
                xmax = 300
                ymin = 1
                ymax = 4
            else:
                xmin = 0.3
                xmax = 300
                ymin = 0.5
                ymax = 16

            xscale='log'
            yscale='log'
            xticks = [0.3,1,3,10,30,100,300]
            xlabel = 'Orbital Period (days)'
            self.bwx = log10(1 + 0.5)
            
        if xk=='giso_sinc':
            if zoom:
                xmin = 3e0
                xmax = 3e3
                ymin = 1.0
                ymax = 4
            else:
                xmin = 1
                xmax = 1e4
                ymin = 0.5
                ymax = 16
            xscale='log'
            yscale='log'
            xticks = [1,3,10,30,100,300,1000,3000,10000]
            xlabel = 'Incident Bolometric Flux (Earth-units)'
            self.bwx = log10(1 + 1)
           
        if xk=='giso_smass':
            if zoom:
                xmin = 0.5
                xmax = 1.4
                ymin = 1
                ymax = 4
            else:
                xmin = 0.5
                xmax = 1.5
                ymin = 0.5
                ymax = 16

            xscale='log'
            yscale='log'
            xticks = [0.5,0.7,1.0,1.4]
            xlabel = 'Stellar Mass (Solar-masses)'
            self.bwx = log10(1 + 0.15)

        if xk=='cks_smet':
            xerr = 0.15
            if zoom:
                xmin = -0.4
                xmax = 0.4
                ymin = 1
                ymax = 4
            else:
                xmin = -0.5
                xmax = 0.5
                ymin = 0.5
                ymax = 16

            xticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
            xlabel = '[Fe/H]'
            self.bwx = 0.1
            xscale='lin'
            yscale='log'

        if xk=='giso_sage':
            if zoom:
                xmin = 1
                xmax = 10
                ymin = 1
                ymax = 4
            else:
                xmin = 1
                xmax = 15
                ymin = 0.5
                ymax = 16

            xscale='log'
            yscale='log'
            xticks = [1,3,10]
            xlabel = 'Stellar Age (Gyr)'
            self.bwx = log10(1 + 0.25)
            

        yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]
        ylabel = 'Planet Size (Earth-radii)'
        cp = ContourPlotter(xmin, xmax, ymin, ymax,xscale=xscale,yscale=yscale)
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

    def plot(self, gradient_array=False, column_normalize=False):
        ds = self.cp.meshgrid()
        Z = gaussian_2d_kde(array(ds.kxc), array(ds.kyc), self.kx,
                            self.ky, self.bwx, self.bwy)
        Z = Z.reshape(ds.kxc.shape)

        if gradient_array:
            return ds['kx'].values, ds['ky'].values, Z
    

        ds['Z'] = (['kx','ky'],Z)
        qc = self.cp.contour(
            ds['Z'], normalize=True, cmap=plt.cm.afmhot_r,zorder=0,
            column_normalize=column_normalize
        )
        self.qc = qc
        self.cp.errorbar(
            self.x, self.y, yerr=None, fmt='o', mfc='none', elinewidth=0,
            ms=2, mew=0.4, mec='w', zorder=9
        )

        self.cp.errorbar(self.x, self.y, yerr=None, fmt='o',
                         mfc='none', elinewidth=0, ms=2, mew=0.2,
                         mec='k', zorder=10)

        self.cp.yticks(self.yticks)
        self.cp.xticks(self.xticks)
        self.cp.yticks(self.yticks)
        self.cp.set_lim()
        setp(gca(),xlabel=self.xlabel,ylabel=self.ylabel)

class ContourPlotter(object):

    """This class facilliates plotting contours where some of the
    variables have been tranformed from linear to logspace. Just keeps
    track of the conversion

    """
    def __init__(self, xmin, xmax, ymin, ymax, xscale='log', yscale='log'):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xscale = xscale
        self.yscale = yscale
        self._kde_nx = 200 # default, can change
        self._kde_ny = 400

    def meshgrid(self):

        if self.xscale=='log':
            kx = np.linspace(log10(self.xmin), log10(self.xmax), self._kde_nx)
            x = 10**kx
        else:
            x = kx = np.linspace(self.xmin, self.xmax, self._kde_nx)
            
        if self.yscale=='log':
            ky = np.linspace(log10(self.ymin), log10(self.ymax), self._kde_ny)
            y = 10**ky
        else:
            y = ky = np.linspace(self.ymin, self.ymax, self._kde_ny)

        coords = {'kx':kx, 'ky':ky}
        _kx, _ky = meshgrid(kx,ky,indexing='ij')
        _x, _y = meshgrid(x,y,indexing='ij')
        data = {
            'kxc': (['kx', 'ky'], _kx),
            'kyc': (['kx', 'ky'], _ky),
            'xc': (['kx', 'ky'], _x),
            'yc': (['kx', 'ky'], _y),
        }
        ds = xr.Dataset(data,coords=coords)
        return ds
       
    def contour(self, Z, column_normalize=False, normalize=False, **kwargs):
        if column_normalize:
            Z = Z / Z.sum('ky')

        if normalize:
            fac = 1.2
            Z /= Z.max()
            Z /= fac
            qc = Z.plot.contourf(
                x='kx', levels=arange(0,1.001,0.05),
                vmax=1,add_colorbar=False,**kwargs)
        else:
            qc = Z.plot.contourf(x='kx', add_colorbar=False,**kwargs)
            



        return qc

    def set_lim(self):
        self.xlim(self.xmin,self.xmax)
        self.ylim(self.ymin,self.ymax)

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

    def errorbar(self, *args, **kwargs):
        x = args[0]
        y = args[1]

        if self.xscale=='log':
            x = log10(x)
        if self.yscale=='log':
            y = log10(y)
        
        errorbar(x,y, **kwargs)

    def plot_gradient(self, grad):
        kx = np.linspace(log10(grad.x0), log10(grad.x1), 100)
        ky = grad.func(grad.out.params, kx)
        plot(kx,ky)
        axvline(log10(grad.x0))
        axvline(log10(grad.x1))
        
    def plot_gradients(self, grads):
        grad = grads[0]
        kx = np.linspace(log10(grad.x0), log10(grad.x1), 100)
        ky = []
        for grad in grads:
            ky += [grad.func(grad.out.params, kx)]
        ky = np.vstack(ky)
        lo, hi = np.percentile(ky,[16,84],axis=0)
        fill_between(kx, lo, hi, color='b', alpha=0.4)
        

        
            
