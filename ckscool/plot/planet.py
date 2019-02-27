import seaborn as sns
from sklearn.neighbors import KernelDensity
from matplotlib.pylab import *
import ckscool.io
import pandas as pd

def fig_per_prad(nopoints=False,zoom=False):
    fig, axL = subplots(figsize=(5,4))
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    sns.set_context('talk')
    xfac = 0.2
    yfac = 1

    x = np.log10(df.koi_period) * xfac
    y = np.log10(df.gdir_prad) * yfac
    yerr = (np.log10(df.gdir_prad + df.gdir_prad_err1) - np.log10(df.gdir_prad)) * yfac
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
    xt = [0.3,1,3,10,30,100,300]
    yt = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.8,8.0,11.3,16]
    xt = np.array(xt)
    yt = np.array(yt)
    xticks(np.log10(xt)*xfac,xt)
    yticks(np.log10(yt)*yfac,yt)

    if not nopoints:
        errorbar(x, y, yerr=yerr,fmt='.',ms=5,zorder=10,elinewidth=1,color='k',mec='white')
    
    xlabel('Orbital Period (days)')
    ylabel('Planet-size (Earth-radii)')
    ylim(np.log10(0.5)*yfac, np.log10(20)*yfac)
    xlim(np.log10(0.3)*xfac, np.log10(300)*xfac)
    if zoom:
        xlim(np.log10(1)*xfac, np.log10(100)*xfac)
        ylim(np.log10(1)*yfac, np.log10(4)*yfac)
    colorbar()
    tight_layout()



def fig_smass_prad(nopoints=False,zoom=False, boost=True):

    """
    Boost: whether to normalize KDE so that each mass bin gets equal weight
    """
    fig, axL = subplots(figsize=(5,4))
    xk = 'giso_smass'
    yk = 'gdir_prad'
    yerrk = 'gdir_prad_err1'
    xfac = 1
    yfac = 1.5
    df = ckscool.io.load_table('ckscool-planets-cuts',cache=1)
    df = df[~df.isany]
    df = df.query('giso_smass < 0.8')
    df2 = ckscool.io.load_table('cksgaia-planets-filtered')
    df2 = df2.query('giso_smass > 0.8')

    df = pd.concat([df,df2])
    df = df.dropna(subset=[xk,yk,yerrk])

    sns.set_context('talk')
    
    x = np.log10(df[xk]) * xfac
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
    X, Y = np.mgrid[xmin:xmax:300j, ymin:ymax:200j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(np.exp(kde.score_samples(positions.T)), X.shape)

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


    xt = [0.3,0.5,0.7,1,1.4,2]
    yt = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.8,8.0,11.3,16]
    xt = np.array(xt)
    yt = np.array(yt)
    xticks(np.log10(xt)*xfac,xt)
    yticks(np.log10(yt)*yfac,yt)
    
    if not nopoints:
        errorbar(x, y, yerr=yerr,fmt='.',ms=5,zorder=10,elinewidth=1,color='k')


    xlabel('Orbital Period (days)')
    xlabel('Stellar Mass (Solar Masses)')
    ylabel('Planet-size (Earth-radii)')
    xlim(xmin, xmax)
    ylim(np.log10(1)*yfac,np.log10(10)*yfac)

    if zoom:
        #xlim(np.log10(1)*xfac, np.log10(100)*xfac)
        ylim(np.log10(1)*yfac, np.log10(4)*yfac)

    axvline(np.log10(0.8)*xfac,zorder=10,color='m',lw=2,alpha=1)
    colorbar()
    tight_layout()


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


