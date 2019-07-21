import string

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy import constants as c
from ckscool.plot.config import *
import ckscool.io

sns.set_style('ticks')
sns.set_color_codes()

from matplotlib.pylab import *
import seaborn as sns
import pandas as pd
from astropy.io import ascii 
import seaborn as sns
import cksspec.io

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from cksspec.plotting import texdict

errorbar_kw = dict(fmt='.',markersize=6,color='b')

figsize = (3.5,4)

def fig_compare(table):
    if table=='ckscool-mann13':
        df = ckscool.io.load_table('ckscool-mann13')
        df = df.groupby('id_koi',as_index=False).first()
        xk = 'cks_steff'
        yk = 'm13_steff'
        label2 = 'M13'

    elif table=='ckscool-dressing13':
        df = ckscool.io.load_table('ckscool-dressing13')
        df = df.groupby('id_koi',as_index=False).first()
        xk = 'cks_steff'
        yk = 'd13_steff'
        label2 = 'D13'

    elif table=='ckscool-brewer18':
        df = ckscool.io.load_table('ckscool-brewer18')
        df = df.groupby('id_koi',as_index=False).first()
        xk = 'cks_steff'
        yk = 'b18_steff'
        label2 = 'B18'
    else:
        assert False, "failed"

    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,
            'xtick.major.size':3.0,
            'xtick.direction': u'in',
            'ytick.direction': u'in',}
    )
    sns.set_context('paper',font_scale=1.1)

    fig = figure(figsize=(7,3.5))
    shape = (10,2)
    ax1 = subplot2grid(shape, (0,0), rowspan=8)
    ax2 = subplot2grid(shape, (8,0), rowspan=2, sharex=ax1)
    axL = [ax1,ax2]

    df = df[~df.isany]
    steff1 = df[xk]
    steff2 = df[yk]
    comparison('steff', steff1, steff2, axL0=axL, label1='SM', label2=label2)

    '''
    steff1 = df[xk]
    steff2 = df[yk]
    comparison('steff', steff1, steff2, axL0=axL, label1='SM', label2=label2)
    '''



def add_anchored(*args,**kwargs):
    ax = gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)

def comparison(key, x1, x2, label1='CKS',label2='Comp',fig0=None, axL0=None):
    ebarkw = dict(fig0=fig0, axL0=axL0, **errorbar_kw)
    x3 = x2 - x1 

    fig, axL = subplots_compare(x1, x2, x3, **ebarkw)
    if key=='steff':
        _ylabel0 = '{} (K) [{}]'.format(texdict[key],label2)
        _xlim0 = (3500,6500) 
        _ylim1 = (-400,400)
        _ylabel1 = '{} - {}'.format(label2,label1)
        _xlabel1 = '{} (K) [{}]'.format(texdict[key],label1)

    if key=='slogg':
        _ylabel0 = '{} (dex) [{}]'.format(texdict[key],label2)
        _xlim0 = (3.0,5.0)
        _ylim1 = (-0.5,0.5)
        _ylabel1 = '{} - {}'.format(label2,label1)
        _xlabel1 = '{} (dex) [{}]'.format(texdict[key],label1)

    if key=='smet':
        _ylabel0 = '{} (dex) [{}]'.format(texdict[key],label2)
        _xlim0 = (-0.6,0.6)
        _ylim1 = (-0.3,0.3)
        _ylabel1 = '{} - {}'.format(label2,label1)
        _xlabel1 = '{} (dex) [{}]'.format(texdict[key],label1)

    setp(axL[0],ylabel=_ylabel0, xlim=_xlim0, ylim=_xlim0)
    setp(axL[1],ylabel=_ylabel1, xlim=_xlim0, ylim=_ylim1, xlabel=_xlabel1)

    for t in axL[0].get_xticklabels():
        t.set_size(0)
    
    axL[1].yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    axL[1].xaxis.set_major_locator(MaxNLocator(7,prune='both'))
    axL[0].yaxis.set_major_locator(MaxNLocator(7,prune='both'))
    fmt = dict(param=key, mdiff=np.mean(x3),sdiff=np.std(x3),ndiff=x3.count())
    sca(axL[0])
    grid()
    one2one(color='g',linestyle='--')
    sca(axL[1])
    axhline(0,color='g',linestyle='--')
    grid()

    if key=='steff':
        s = """\
Mean($\Delta$) = {mdiff:.0f} K
RMS($\Delta$) = {sdiff:.0f} K""".format(**fmt)

    if key=='slogg':
        s = """\
Mean($\Delta$) = {mdiff:.2f} dex
RMS($\Delta$) = {sdiff:.2f} dex""".format(**fmt)

    if key=='smet':
        s = """\
Mean($\Delta$) = {mdiff:.3f} dex
RMS($\Delta$) = {sdiff:.3f} dex""".format(**fmt)
        
    print " ".join("{}: {}".format(k,fmt[k]) for k in 'param ndiff mdiff sdiff'.split())

    sca(axL[0])
    add_anchored(s,loc=2, prop=dict(size='small'))
    fig.set_tight_layout(True)


def provision_figure():
    fig = figure(figsize=figsize)
    ax1 = subplot2grid((4,1), (0,0), rowspan=3)
    ax2 = subplot2grid((4,1), (3,0), rowspan=1, sharex=ax1)
    axL = [ax1,ax2]
    return fig,axL 

def subplots_compare(x1, x2, x3, xerr3=None, fig0=None, axL0=None, **kwargs):
    if fig0 is None:
        fig, axL = provision_figure()
    else:
        fig = fig0
        axL = axL0

    fig.set_tight_layout(False)

    #fig.subplots_adjust(hspace=0.4,left=0.17,top=0.95,right=0.90)
    fig.subplots_adjust(hspace=0.001,left=0.17,top=0.95,right=0.90)
    sca(axL[0])
    errorbar(x1,x2,**kwargs)
    sca(axL[1])
    errorbar(x1,x3,**kwargs)
    return fig,axL

def one2one(**kwargs):
    xl = xlim()
    plot(xl,xl,**kwargs)

