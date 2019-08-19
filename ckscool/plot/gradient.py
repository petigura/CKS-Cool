import ckscool.comp
import ckscool.io
import ckscool.grid
import ckscool.occur
import ckscool.plot.occur
from ckscool.gradient import R

import numpy as np
from numpy import log10, logspace, sqrt
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
import corner


def occurence_gradient_plot(smass_bins=[0.5,0.8,1.1,1.4]):

    for i in range(len(smass_bins)-1):

        smass1, smass2 = smass_bins[i], smass_bins[i+1]
        chain_array = np.loadtxt("./data/chain_occ_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')

        plt.figure(i+1)

        labels = ["m", r"R$_p$(10)"]
        figure = corner.corner(chain_array[0:,:],
                               quantiles=[0.16, 0.5, 0.84],
                               show_titles=True,
                               labels=labels,
                               smooth1d=True)

        m, m_std = np.mean(chain_array[:,0]), np.std(chain_array[:,0])
        Rp10, Rp10_std = np.mean(chain_array[:,1]), np.std(chain_array[:,1])
        plt.suptitle(r'M$_*$ = {0} - {1} M$_\odot$'.format(smass1, smass2), fontsize=16)
        print m, m_std, Rp10, Rp10_std


        plt.figure(0)

        plt.subplot(len(smass_bins)-1,1,i+1)
        occ = ckscool.occur.load_occur({'smass1':smass1,'smass2':smass2})
        cp = ckscool.plot.occur.load_contour_plotter(occ)
        ckscool.plot.occur.contour(cp)

        P_domain = np.log10(np.logspace(0,2,100))
        plt.plot(P_domain, [log10(R(10**i, m, Rp10)) for i in P_domain],color='b', label=r"m = {0:.2f}$\pm${1:.2f}, R$_p$(10)={2:.2f}$\pm${3:.2f}".format(m, m_std, Rp10, Rp10_std))

        # error representation
        N_grid = len(P_domain)
        midpoint = P_domain[N_grid//2-1]
        # match up shaded region at 10 day period
        eps = log10( R(10**midpoint,m-m_std,Rp10+Rp10_std) / R(10**midpoint,m+m_std,Rp10+Rp10_std) )

        plt.fill_between(P_domain[:N_grid//2],
                             [log10(R(10**i,m+m_std,Rp10-Rp10_std)) for i in P_domain[:N_grid//2]],
                             [log10(R(10**i,m-m_std,Rp10+Rp10_std)) for i in P_domain[:N_grid//2]],
                             color='b',
                             linewidth = 0.0,
                             alpha = 0.2)

        plt.fill_between(P_domain[N_grid//2-1:],
                             [log10(R(10**i,m-m_std,Rp10-Rp10_std))-eps for i in P_domain[N_grid//2-1:]],
                             [log10(R(10**i,m+m_std,Rp10+Rp10_std))+eps for i in P_domain[N_grid//2-1:]],
                             color='b',
                             linewidth = 0.0,
                             alpha = 0.2)

        plt.title(r'M$_*$ = {0} - {1} M$_\odot$'.format(smass1, smass2), fontsize=16)
        plt.legend(loc=1)

    plt.show()

# ---------------------------------------------------------------------------- #


def detections_gradient_plot(smass_bins=[0.5,0.8,1.1,1.4]):

    for i in range(len(smass_bins)-1):

        smass1, smass2 = smass_bins[i], smass_bins[i+1]
        chain_array = np.loadtxt("./data/chain_detections_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')

        plt.figure(i+1)

        labels = ["m", r"R$_p$(10)"]
        figure = corner.corner(chain_array[0:,:],
                               quantiles=[0.16, 0.5, 0.84],
                               show_titles=True,
                               labels=labels,
                               smooth1d=True)

        m, m_std = np.mean(chain_array[:,0]), np.std(chain_array[:,0])
        Rp10, Rp10_std = np.mean(chain_array[:,1]), np.std(chain_array[:,1])
        plt.suptitle(r'M$_*$ = {0} - {1} M$_\odot$'.format(smass1, smass2), fontsize=16)
        print m, m_std, Rp10, Rp10_std


        plt.figure(0)

        plt.subplot(len(smass_bins)-1,1,i+1)
        ckscool.plot.planet.fig_per_prad(zoom=True, query="{0}<=giso_smass<={1}".format(smass1,smass2))

        P_domain = np.log10(np.logspace(0,2,100))
        plt.plot(P_domain, [log10(R(10**i, m, Rp10)) for i in P_domain],color='b', label=r"m = {0:.2f}$\pm${1:.2f}, R$_p$(10)={2:.2f}$\pm${3:.2f}".format(m, m_std, Rp10, Rp10_std))

        # error representation
        N_grid = len(P_domain)
        midpoint = P_domain[N_grid//2-1]
        # match up shaded region at 10 day period
        eps = log10( R(10**midpoint,m-m_std,Rp10+Rp10_std) / R(10**midpoint,m+m_std,Rp10+Rp10_std) )

        plt.fill_between(P_domain[:N_grid//2],
                             [log10(R(10**i,m+m_std,Rp10-Rp10_std)) for i in P_domain[:N_grid//2]],
                             [log10(R(10**i,m-m_std,Rp10+Rp10_std)) for i in P_domain[:N_grid//2]],
                             color='b',
                             linewidth = 0.0,
                             alpha = 0.2)

        plt.fill_between(P_domain[N_grid//2-1:],
                             [log10(R(10**i,m-m_std,Rp10-Rp10_std))-eps for i in P_domain[N_grid//2-1:]],
                             [log10(R(10**i,m+m_std,Rp10+Rp10_std))+eps for i in P_domain[N_grid//2-1:]],
                             color='b',
                             linewidth = 0.0,
                             alpha = 0.2)

        plt.title(r'M$_*$ = {0} - {1} M$_\odot$'.format(smass1, smass2), fontsize=16)
        plt.legend(loc=1)

    plt.show()
