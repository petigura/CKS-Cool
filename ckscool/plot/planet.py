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

    ndplotkey = ['koi_period','giso_sinc','giso_smass','cks_smet']
    labels = 'abcd'

    gradkey = ['grad-per-prad-det_smass=0.5-1.4','grad-sinc-prad-det_smass=0.5-1.4','grad-smass-prad-det_smass=0.5-1.4','grad-smet-prad-det_smass=0.5-1.4']

    for i in range(4):
        ax = axL.flatten()[i]
        sca(ax)
        
        pl = NDPlotter(df,ndplotkey[i],smass_lims=[0.5,1.4],**kwargs)
        pl.plot()
        
        fig_label(labels[i])
        if plot_gradient:
            grad = ckscool.gradient.Gradient(gradkey[i])
            grad.load_csv()
            grad.plot_gradients('band')

        if i==0:
            ax = axes([0.89, 0.85, 0.008, 0.1])
            cbar = colorbar(pl.qc,cax=ax,format='%.2f')
            cbar.set_label('relative density',size='x-small')
            cbar.ax.tick_params(labelsize=0.4 * rcParams['font.size'])

        if i==1:
            xl = xlim()
            xlim(xl[1],xl[0])

    tight_layout(True)


def fig_sample_smass():
    sns.set_context('paper',font_scale=1.1)
    fig,ax = subplots(ncols=1,nrows=1,figsize=(6,4))

    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]

    ndplotkey = 'giso_smass'
    gradkey = 'grad-smass-prad-det_smass=0.5-1.4'

    
    sca(ax)

    pl = NDPlotter(df,ndplotkey,zoom=True)
    pl.plot()

    grad = ckscool.gradient.Gradient(gradkey)
    grad.load_csv()
    #ylim(1,4.0)
    xl = xlim()
    logsmass= linspace(xl[0],xl[1],100)

    grad.plot_gradients('band',xi =logsmass )

    logprad = log10(2.5) + 0.25 * logsmass 
    plot(logsmass, logprad,'r-')
    text(log10(1.45),log10(2.8),r'$\alpha = 0.25$')
    
    logprad = log10(1.4) + 0.25 * logsmass 
    plot(logsmass, logprad,'r-')
    text(log10(1.45),log10(1.5),r'$\alpha = 0.25$')

    mps = ckscool.io.load_object('mps_size-se',cache=1)
    fill_between(mps.logsmassci, log10(mps.q16),log10(mps.q84) ,alpha=0.5)


    mps = ckscool.io.load_object('mps_size-sn',cache=1)
    fill_between(mps.logsmassci, log10(mps.q16),log10(mps.q84) ,alpha=0.5)


def fig_gupta_comparison(plot_gradient=False, **kwargs):
    sns.set_context('paper',font_scale=1.3)
    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,'xtick.major.size':3.0,
            'xtick.direction': u'in','ytick.direction': u'in'
        }
    )
    #sns.set_style('whitegrid')
    fig,axL = subplots(ncols=2,nrows=2,figsize=(8,5))
    #fig,axL = subplots(ncols=2,nrows=2,figsize=(10,6.25))
    kwline = dict(linestyle='--',color='white',zorder=10,lw=1)
    kwline2 = dict(linestyle='-',color='black',zorder=9,lw=2)
    def _plot(*args):
        plot(xi,yi,**kwline)
        plot(xi,yi,**kwline2)
        
    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]

    # Metalicity
    pl = NDPlotter(df,'cks_smet',zoom=True)
    xi = linspace(-0.4,0.4)
    yi = 0.1 * xi + log10(2.3)
    xt = [-0.4,-0.2,0.0,0.2,0.4]
    yt = [1.0,1.5,2.4,3.5]
    
    sca(axL[0,0])
    pl.plot(add_colorbar=True,cbar_title='relative density of planets')
    pl.cp.xlim(-0.5,0.5)
    pl.cp.xticks(xt)
    pl.cp.yticks(yt)
    _plot(xi,yi)

    sca(axL[1,0])
    pl.cp.xmin=-0.4
    pl.cp.xmax=0.4
    pl.plot(add_colorbar=True,column_normalize=True, cbar_title='relative density of planets \n (at const [Fe/H])')
    pl.cp.xlim(-0.5,0.5)
    pl.cp.xticks(xt)
    pl.cp.yticks(yt)
    _plot(xi,yi)
    
    # Age
    xt = [1,3,10]

    df = ckscool.io.load_table('planets-cuts2',cache=1)
    df = df[~df.isany]
    df = df.query('cks_steff > 5500')
    pl = NDPlotter(df,'giso_sage',zoom=True)

    xi = linspace(log10(1),log10(10))
    yi = -0.1 * xi + log10(2.5)
    sca(axL[0,1])
    pl.cp.xmin=1.2
    pl.cp.xmax=12
    pl.plot(add_colorbar=True,cbar_title='relative density of planets')
    pl.cp.xlim(0.5,12)
    pl.cp.xticks(xt)
    pl.cp.yticks(yt)
    _plot(xi,yi)

    sca(axL[1,1])

    pl.cp.xmin=1.2
    pl.cp.xmax=12
    pl.plot(add_colorbar=True,column_normalize=True,cbar_title='relative density of planets \n (at const. age)')
    pl.cp.xlim(0.5,12)
    _plot(xi,yi)
    tight_layout(True)

    labels = 'cdef'
    for i in range(4):
        sca(axL.flatten()[i])
        grid()
        fig_label(labels[i])
    
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

    def plot(self, gradient_array=False, column_normalize=False, **kwargs):
        ds = self.cp.meshgrid()
        Z = gaussian_2d_kde(array(ds.kxc), array(ds.kyc), self.kx,
                            self.ky, self.bwx, self.bwy)
        Z = Z.reshape(ds.kxc.shape)

        if gradient_array:
            return ds['kx'].values, ds['ky'].values, Z
    

        ds['Z'] = (['kx','ky'],Z)
        qc = self.cp.contour(
            ds['Z'], normalize=True, cmap=plt.cm.afmhot_r,zorder=0,
            column_normalize=column_normalize, **kwargs
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
       
    def contour(self, Z, column_normalize=False, normalize=False, add_colorbar=False, cbar_title=None,  **kwargs):
        if column_normalize:
            Z = Z / Z.sum('ky')

        if normalize:
            fac = 1.2 # Fac is a factor that determins how saturated the darkest contour is 
            Z /= Z.max()
            Z /= fac
            qc = Z.plot.contourf(
                x='kx', levels=arange(0,fac + 0.001,0.05),
                vmax=1,add_colorbar=False,**kwargs)


        else:
            qc = Z.plot.contourf(x='kx', add_colorbar=False,**kwargs)
            
        if add_colorbar:
            cbar = colorbar(qc,format='%.2f',shrink=0.5)
            cbar.set_label(cbar_title,size='small')
            cbar.ax.tick_params(labelsize=0.6 * rcParams['font.size'])

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
        

        
            
