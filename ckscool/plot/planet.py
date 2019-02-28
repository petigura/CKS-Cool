import seaborn as sns
from sklearn.neighbors import KernelDensity
from matplotlib.pylab import *
import ckscool.io
import pandas as pd
import xarray as xr
import scipy
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText


class ContourPlotter(object):
    def __init__(self, x, xerr, y, yerr, xscale='log', yscale='log'):
        self._kde_nx = 100
        self._kde_ny = 100
        self.xscale = xscale
        self.yscale = xscale

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
        kx = np.linspace(log10(self.xmin), log10(self.xmax), self._kde_nx)
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

    def plot_contour(self):
        self.ds.Z.plot.contourf(x='kx', cmap=plt.cm.afmhot_r, levels=20)
        add_anchored('$N_p$ = {}'.format(len(self.x)),2,prop=dict(size='small'))


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

def fig_per_prad(nopoints=False,zoom=False):
    sns.set_context('paper')
    fig, axL = subplots(figsize=(5,4))
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    x = df['koi_period']
    xerr = x * 1.1
    y = df['gdir_prad']
    yerr = df['gdir_prad_err1']

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
        p1.plot_points()
        p1.xlim(1,100)
        p1.ylim(1,4)


def fig_smass_prad(nopoints=False,zoom=False, boost=True):
    """
    Boost: whether to normalize KDE so that each mass bin gets equal weight
    """
    sns.set_context('paper')
    fig, axL = subplots(figsize=(5,4))
    xk = 'giso_smass'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    df = df.query('giso_smass < 0.8')
    df2 = ckscool.io.load_table('cksgaia-planets-filtered')
    df2 = df2.query('giso_smass > 0.8')
    df = pd.concat([df,df2])
    df = df2
    df = df.dropna(subset=[xk,yk,yerrk])

    x = df[xk]
    y = df[yk]
    yerr = df[yerrk]
    xerr = x * 0.15

    p1 = ContourPlotter(x, xerr, y, yerr,xscale='log',yscale='log')
    p1.xmin = 0.5
    p1.xmax = 1.5
    p1.ymin = 0.5
    p1.ymax = 16
    xticks = [0.5,0.7,1.0,1.4]
    yticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]

    p1.compute_density()
    p1.plot_contour()
    xlabel('Stellar Mass (Solar-masses)')
    ylabel('Planet Size (Earth-radii)')
    p1.xticks(xticks)
    p1.yticks(yticks)
    if not nopoints:
        p1.plot_points()
    if zoom:
        p1.plot_points()
        p1.xlim(0.5,1.5)
        p1.ylim(1,4)
 




    '''
    if boost:
        #imshow(np.rot90((Z/fac[:,np.newaxis])),vmax=3,aspect='auto')
        #import pdb;pdb.set_trace()
        Z = Z / Z.sum(axis=1)
        imshow(
            np.rot90(Z), cmap=plt.cm.afmhot_r, extent=[xmin, xmax, ymin, ymax],
            aspect='auto',zorder=1,vmax=8
        )
    else:
        imshow(
            np.rot90(Z), cmap=plt.cm.afmhot_r, extent=[xmin, xmax, ymin, ymax],
            aspect='auto',zorder=1,vmax=8
        )

    '''
def fig_smet_prad(nopoints=False,zoom=False):
    fig, axL = subplots(figsize=(5,4))
    xk = 'sm_smet'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    xfac = 0.3
    yfac = 1.25
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    df = df.dropna(subset=[xk,yk,yerrk])
    sns.set_context('talk')
    
    x = df[xk] * xfac
    y = np.log10(df[yk]) * yfac
    yerr = (np.log10(df[yk] + df[yerrk]) - np.log10(df[yk])) * yfac
    xy = np.vstack([x,y])

    d = xy.shape[0]
    n = xy.shape[1]

    bw = 0.03
    print('bw: {}'.format(bw))


    kde = KernelDensity(
        bandwidth=bw, metric='euclidean', kernel='gaussian', algorithm='ball_tree')
    kde.fit(xy.T)

    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    X, Y = np.mgrid[xmin:xmax:300j, ymin:ymax:300j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(np.exp(kde.score_samples(positions.T)), X.shape)

    imshow(
        np.rot90(Z), cmap=plt.cm.afmhot_r, extent=[xmin, xmax, ymin, ymax],
        aspect='auto',zorder=1,vmax=12
    )

    xt = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]
    yt = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.8,8.0,11.3,16]
    xt = np.array(xt)
    yt = np.array(yt)
    xticks(xt*xfac,xt)
    yticks(np.log10(yt)*yfac,yt)

    if not nopoints:
        errorbar(x, y, yerr=yerr,fmt='.',ms=5,zorder=10,elinewidth=1,color='k')

    xlabel('[Fe/H]')
    ylabel('Planet-size (Earth-radii)')
    xlim(xmin, xmax)
    ylim(ymin, ymax)

    if zoom:
        ylim(np.log10(1)*yfac, np.log10(4)*yfac)

    colorbar()
    tight_layout()

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
