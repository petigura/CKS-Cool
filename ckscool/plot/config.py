from matplotlib.pylab import * 
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import seaborn as sns


def add_anchored(*args,**kwargs):
    """
    Parameters
    ----------
    s : string
        Text.

    loc : str
        Location code.

    pad : float, optional
        Pad between the text and the frame as fraction of the font
        size.

    borderpad : float, optional
        Pad between the frame and the axes (or *bbox_to_anchor*).

    prop : `matplotlib.font_manager.FontProperties`
        Font properties.
    """

    bbox = {}
    if kwargs.has_key('bbox'):
        bbox = kwargs.pop('bbox')
    at = AnchoredText(*args, **kwargs)
    if len(bbox.keys())>0:
        plt.setp(at.patch,**bbox)

    ax = plt.gca()
    ax.add_artist(at)

def fig_label(text):
    add_anchored(
        text, loc=2, frameon=True, 
        prop=dict(size='large', weight='bold'),
        bbox=dict(ec='none', fc='w', alpha=0.0)
    )

def sns_set_style(key):
    if key=='ticks':
        sns.set(
            style='ticks',
            rc={'ytick.major.size':3.0,'xtick.major.size':3.0,
                'xtick.direction': u'out','ytick.direction': u'out'
            }
        )
        sns.set_context('paper',font_scale=1.1)
    if key=='whitegrid':
        sns.set_style('whitegrid')
        sns.set_context('paper',font_scale=1.1)
