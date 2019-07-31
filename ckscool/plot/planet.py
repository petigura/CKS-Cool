import seaborn as sns
from sklearn.neighbors import KernelDensity
from matplotlib.pylab import *
import ckscool.io
import pandas as pd
import xarray as xr
import scipy
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from astropy import units as u


class ContourPlotter(object):
    def __init__(self, x, xerr, y, yerr, xscale='log', yscale='log'):
        self._kde_nx = 200
        self._kde_ny = 400
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
        
    def compute_density(self):
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
        self.ds = ds

    def normalize_density(self):
        frac = self.ds['Z'].sum(dim='ky')
        frac = frac / frac.mean()
        self.ds['Z'] /= frac

    def plot_contour(self):
        self.ds.Z.plot.contourf(x='kx', cmap=plt.cm.afmhot_r, levels=20)
        add_anchored('$N_p$ = {}'.format(len(self.x)),2,prop=dict(size='small'))
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

    def plot_integrated(self):
        pass

    def plot_points(self):
        errorbar(self.x,self.y,yerr=self.yerr,fmt='.',elinewidth=0.5,ms=4)

def fig_per_prad(**kwargs):
    sns.set_context('paper',font_scale=1.1)
    fig, axL = subplots(figsize=(5,4))
    df = ckscool.io.load_table('planets-cuts2+iso',cache=2)
    df = df[~df.isany]
    _per_prad(df, **kwargs)

def _per_prad(df,nopoints=False,zoom=False,query=None,yerrfac=1,xerrfac=1):
    if query is not None:
        df = df.query(query)
    print len(df)
    xk = 'koi_period'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    x = df[xk]
    y = df[yk]
    xerr = x * xerrfac
    yerr = df[yerrk] * yerrfac

    p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
    p1.xmin = 0.3
    p1.xmax = 300
    p1.ymin = 0.5
    p1.ymax = 16
    xticks = [0.3,1,3,10,30,100,300]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.plot_contour()
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    p1.xticks(xticks)
    p1.yticks(yticks)

    if not nopoints:
        p1.plot_points()
    if zoom:
        p1.xlim(1,100)
        p1.ylim(1,4)

    tight_layout()


def fig_sinc_prad(nopoints=False,zoom=False,query=None,yerrfac=1,xerrfac=1):
    df = ckscool.io.load_table('planets-cuts2+iso',cache=1)
    df = df[~df.isany]
    if query is not None:
        df = df.query(query)

    print len(df)
    sns.set_context('paper')
    fig, axL = subplots(figsize=(5,4))
    xk = 'giso_sinc'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'

    x = df[xk]
    y = df[yk]
    xerr = x * xerrfac
    yerr = df[yerrk] * yerrfac

    p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
    if zoom:
        p1.xmin = 1
        p1.xmax = 1e4
        p1.ymin = 1.0
        p1.ymax = 4
    else:
        p1.xmin = 1
        p1.xmax = 1e4
        p1.ymin = 0.5
        p1.ymax = 16

    xticks = [1,3,10,30,100,300,1000,3000,1e4]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.plot_contour()
    xlabel('Incident Bolometric Flux (Earth-units)')
    ylabel('Planet Size (Earth-radii)')
    p1.xticks(xticks)
    p1.yticks(yticks)

    p1.xlim(p1.xmax,p1.xmin)
    p1.ylim(p1.ymin,p1.ymax)

    if not nopoints:
        p1.plot_points()

    tight_layout()


def fig_intfxuv_prad(nopoints=False,zoom=False):
    sns.set_context('paper')
    fig, axL = subplots(figsize=(5,4))
    xk = 'giso_sintfxuv'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    df2 = ckscool.io.load_table('cksgaia-planets-filtered')
    namemap = {'giso_insol':'giso_sinc','giso_insol_err1':'giso_sinc_err1','giso_insol_err2':'giso_sinc_err2',}
    df2 = df2.rename(columns=namemap)
    df = pd.concat([df,df2])

    df =  df.dropna(subset=['giso_smass'])
    df = df.query('giso_smass < 0.8')

    data = np.load('mcdonald/posterior_padova125_1000_lookup.npy')
    mass = pd.read_csv('mcdonald/posterior_mass_vec.csv',names=['mass'],squeeze=True)
    logage = pd.read_csv('mcdonald/posterior_logage_vec.csv',names=['logage'],squeeze=True)
    coords={'mass':mass,'logage':logage,'sample':arange(1000)}
    dims = ('mass','logage','sample')
    data = xr.DataArray(data,coords=coords,dims=dims)
    med = data.median(dim='sample')

    logage = ones(len(df))*log10(5e9)
    logage = log10(df.giso_sage * 1e9)
    logage = xr.DataArray(logage, dims='z')
    mass = xr.DataArray(df['giso_smass'], dims='z')
    df['giso_sintlxuv'] = med.interp(mass=mass, logage=logage,kwargs={'fill_value': None})
    sma = np.array(df['giso_sma']) * u.AU
    intlxuv = np.array(df['giso_sintlxuv']) * u.erg
    df['giso_sintfxuv'] = intlxuv / 4 / pi / sma**2

    #df = df.query('giso_smass < 0.8')

    x = df[xk]
    y = df[yk]
    xerr = x * 1
    yerr = df[yerrk] * 1.5

    p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
    if zoom:
        p1.xmin = 1e45
        p1.xmax = 1e49
        p1.ymin = 0.7
        p1.ymax = 4
    else:
        p1.xmin = 1e45
        p1.xmax = 1e49
        p1.ymin = 0.5
        p1.ymax = 16

    xticks = [1e45,1e46,1e47,1e48,1e49]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.plot_contour()
    xlabel(r'$\int F_{xuv} dt$ (ergs/cm^2)')
    ylabel('Planet Size (Earth-radii)')
    p1.xticks(xticks)
    p1.yticks(yticks)

    p1.xlim(p1.xmax,p1.xmin)
    p1.ylim(p1.ymin,p1.ymax)

    if not nopoints:
        p1.plot_points()

    tight_layout()

def fig_smass_prad(nopoints=False,zoom=False, boost=True):
    """
    Boost: whether to normalize KDE so that each mass bin gets equal weight
    """
    sns.set_context('paper')
    fig, axL = subplots(figsize=(5,4))
    xk = 'giso_smass'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    df = ckscool.io.load_table('planets-cuts2+iso',cache=1)
    df = df.dropna(subset=['giso_smass'])
    df = df[~df.isany]

    x = df[xk]
    y = df[yk]
    yerr = df[yerrk] * 1.5
    xerr = x * 0.1

    p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')

    if zoom:
        p1.xmin = 0.6
        p1.xmax = 1.4
        p1.ymin = 0.7
        p1.ymax = 4
    else:
        p1.xmin = 0.5
        p1.xmax = 1.5
        p1.ymin = 0.5
        p1.ymax = 16

    xticks = [0.5,0.7,1.0,1.4]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.normalize_density()
    p1.xticks(xticks)
    p1.yticks(yticks)
    p1.plot_contour()
    xlabel('Stellar Mass (Solar-masses)')
    ylabel('Planet Size (Earth-radii)')
    p1.xlim(p1.xmin,p1.xmax)
    p1.ylim(p1.ymin,p1.ymax)

    if not nopoints:
        p1.plot_points()

    tight_layout()


def fig_smet_prad(nopoints=False,zoom=False):
    sns.set_context('paper',font_scale=1.1)
    fig, axL = subplots(figsize=(5,4))
    xk = 'cks_smet'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    df = ckscool.io.load_table('planets-cuts2+iso',cache=1)
    df = df[~df.isany]
    df = df.dropna(subset=[xk,yk,yerrk])
    
    x = df[xk] 
    y = df[yk]
    yerr = df[yerrk] * 1.5
    xerr = 0.1

    p1 = ContourPlotter(x, xerr, y, yerr, xscale='lin', yscale='log')
    if zoom:
        p1.xmin = -0.4
        p1.xmax = 0.5
        p1.ymin = 1
        p1.ymax = 4
    else:
        p1.xmin = -0.5
        p1.xmax = 0.5
        p1.ymin = 0.5
        p1.ymax = 16

    xticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.normalize_density()

    p1.plot_contour()
    xlabel('[Fe/H]')
    ylabel('Planet-size (Earth-radii)')
    p1.xticks(xticks)
    p1.yticks(yticks)

    if not nopoints:
        p1.plot_points()
    if zoom:
        p1.plot_points()
        p1.xlim(-0.5,0.5)
        p1.ylim(1,4)

    #tight_layout()


def add_anchored(*args,**kwargs):
    ax = gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)

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
