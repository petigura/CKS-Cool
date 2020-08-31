"""
Module to assist with contour plotting
"""
from matplotlib.pylab import *
import xarray as xr

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
        return ds
        
    def contour(self, Z, normalize=False, **kwargs):
        if normalize:
            fac = 1.2
            Z /= Z.max()
            Z /= fac
            qc = Z.plot.contourf(
                x='kx', levels=arange(0,1.001,0.05),zorder=0,
                vmax=1,add_colorbar=False,**kwargs)
        else:
            qc = Z.plot.contourf(x='kx', zorder=0,add_colorbar=False,**kwargs)
            
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
