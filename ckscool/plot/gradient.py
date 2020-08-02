from __future__ import division

import ckscool.comp
import ckscool.io
import ckscool.grid
import ckscool.occur
import ckscool.plot.occur
from ckscool.gradient import R, R_sinc

import numpy as np
from numpy import log, log10, logspace, sqrt, random, exp
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
import corner


def occurence_gradient_plot(smass_bins, corner_plot=False, sinc=False):

    [smass1, smass2] = smass_bins

    if sinc:
        chain_array = np.loadtxt("./data/chain_occ_SincPrad_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')
    else:
        chain_array = np.loadtxt("./data/chain_occ_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')

    # if corner_plot:
    #     plt.figure(i+1)
    #     labels = ["m", r"R$_p$(10)"]
    #     figure = corner.corner(chain_array[0:,:],
    #                            quantiles=[0.16, 0.5, 0.84],
    #                            show_titles=True,
    #                            labels=labels,
    #                            smooth1d=True)

    m = np.median(chain_array[:,0])
    m_sigma = np.percentile(chain_array[:,0], [16,84])

    Rp10 = np.median(chain_array[:,1])
    Rp10_sigma = np.percentile(chain_array[:,1], [16,84])

    m_err = [m-m_sigma[0], m_sigma[1]-m]
    Rp10_err = [Rp10-Rp10_sigma[0], Rp10_sigma[1]-Rp10]

    print r'm = {0:.3f} + {1:.3f} / - {2:.3f}'.format(m, m_err[0], m_err[1])
    print r'Rp10 = {0:.3f} + {1:.3f} / - {2:.3f}'.format(Rp10, Rp10_err[0], Rp10_err[1])

    if sinc:
        if os.path.exists("./data/cp_sinc_{0}-{1}.pkl".format(smass_bins[0], smass_bins[1])):
            f = open("./data/cp_sinc_{0}-{1}.pkl".format(smass_bins[0], smass_bins[1]),'rb')
            cp = pickle.load(f)
            f.close()
        else:
            occ = ckscool.occur.load_occur({'smass1':smass1,'smass2':smass2}, sinc=True)
            cp = ckscool.plot.occur.load_contour_plotter_sinc(occ)

        ckscool.plot.occur.contour_sinc(cp, ntrials_min=100, normalize=True)
        P_domain = np.log10(np.logspace(0,4,100))
    else:
        if os.path.exists("./data/cp_{0}-{1}.pkl".format(smass_bins[0], smass_bins[1])):
            f = open("./data/cp_{0}-{1}.pkl".format(smass_bins[0], smass_bins[1]),'rb')
            cp = pickle.load(f)
            f.close()
        else:
            occ = ckscool.occur.load_occur({'smass1':smass1,'smass2':smass2}, sinc=True)
            cp = ckscool.plot.occur.load_contour_plotter_sinc(occ)

        ckscool.plot.occur.contour(cp, ntrials_min=100, normalize=True)
        P_domain = np.log10(np.logspace(0,2,100))

    label_str = r"m = {0:.2f}$^{{{1:.2f}}}_{{{2:.2f}}}$, R$_p$(10)={3:.2f}$^{{{4:.2f}}}_{{{5:.2f}}}$".format(m, m_err[1], m_err[0], Rp10, Rp10_err[1], Rp10_err[0])
    plt.plot(P_domain, [log10(R(10**i, m, Rp10)) for i in P_domain], color='b', label=label_str)

    R_upper = []
    R_lower = []

    for i in P_domain:
        R_values_i = [R(10**i, chain_array[j,0], chain_array[j,1]) for j in range(len(chain_array[:,0]))]
        R_lim_i = np.percentile(R_values_i, [16,84])
        R_upper.append(R_lim_i[0])
        R_lower.append(R_lim_i[1])

    plt.fill_between(P_domain, log10(R_lower), log10(R_upper), color='b', alpha=0.2)
    plt.title(r'M$_*$ = {0} - {1} M$_\odot$'.format(smass1, smass2), fontsize=16)
    plt.legend(loc=1)

# ---------------------------------------------------------------------------- #


def detections_gradient_plot(smass_bins=[0.5,0.8,1.1,1.4], corner_plot=False):

    for i in range(len(smass_bins)-1):

        smass1, smass2 = smass_bins[i], smass_bins[i+1]
        chain_array = np.loadtxt("./data/chain_detections_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')

        if corner_plot:
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


def gradient_panels(sinc=False):

    smassbins = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1]
    smass_mid = np.zeros(len(smassbins))

    m_det_mean = np.zeros(len(smassbins))
    m_det_std = np.zeros(len(smassbins))
    Rp10_det_mean = np.zeros(len(smassbins))
    Rp10_det_std = np.zeros(len(smassbins))

    m_occ_mean = np.zeros(len(smassbins))
    m_occ_std = np.zeros(len(smassbins))
    Rp10_occ_mean = np.zeros(len(smassbins))
    Rp10_occ_std = np.zeros(len(smassbins))


    for i in range(len(smassbins)):

        smass1, smass2 = smassbins[i], smassbins[i] + 0.3
        smass_mid[i] = ( smass1 + smass2 ) / 2

        if sinc:
            occ_array = np.loadtxt("./data/chain_occ_SincPradv2_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')
            det_array = np.loadtxt("./data/chain_det_SincPradv2_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')
        else:
            occ_array = np.loadtxt("./data/chain_occv2_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')
            det_array = np.loadtxt("./data/chain_detv2_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')

        # occurrence
        m_occ_i, m_occ_std_i = np.mean(occ_array[:,0]), np.std(occ_array[:,0])
        Rp10_occ_i, Rp10_occ_std_i = np.mean(occ_array[:,1]), np.std(occ_array[:,1])

        m_det_i, m_det_std_i = np.mean(det_array[:,0]), np.std(det_array[:,0])
        Rp10_det_i, Rp10_det_std_i = np.mean(det_array[:,1]), np.std(det_array[:,1])

        m_occ_mean[i] = m_occ_i
        m_occ_std[i] = m_occ_std_i
        Rp10_occ_mean[i] = Rp10_occ_i
        Rp10_occ_std[i] = Rp10_occ_std_i

        m_det_mean[i] = m_det_i
        m_det_std[i] = m_det_std_i
        Rp10_det_mean[i] = Rp10_det_i
        Rp10_det_std[i] = Rp10_det_std_i

    # -------------------------- DETECTIONS + M ------------------------------ #
    plt.subplot(221)
    plt.title('Detections', fontsize=14)
    plt.errorbar(smass_mid, m_det_mean, m_det_std, capsize=2)

    if sinc:
        plt.ylim([-0.05,0.15])
    else:
        plt.ylim([-0.24,0.05])

    plt.ylabel('m', fontsize=14)
    plt.grid()

    # ----------------------- DETECTIONS + Rp10 ------------------------------ #
    plt.subplot(223)
    plt.errorbar(smass_mid, Rp10_det_mean, Rp10_det_std, capsize=2)
    plt.xlabel(r'Stellar Mass [M$_\odot$]', fontsize=14)

    if sinc:
        plt.ylim([1.5,2.25])
        plt.ylabel(r'R$_p$(100)', fontsize=14)
    else:
        plt.ylim([1.5,2.25])
        plt.ylabel(r'R$_p$(10)', fontsize=14)

    plt.grid()

    # -------------------------- OCCURRENCE + M ------------------------------ #
    plt.subplot(222)
    plt.title('Occurence', fontsize=14)
    plt.errorbar(smass_mid, m_occ_mean, m_occ_std, capsize=2)

    if sinc:
        plt.ylim([-0.05,0.15])
    else:
        plt.ylim([-0.24,0.05])

    plt.grid()

    # ----------------------- OCCURRENCE + Rp10 ------------------------------ #
    plt.subplot(224)
    plt.errorbar(smass_mid, Rp10_occ_mean, Rp10_occ_std, capsize=2)

    if sinc:
        plt.ylim([1.5,2.25])
        plt.suptitle('Incident Stellar Flux vs. Planet Size', fontsize=20)
    else:
        plt.ylim([1.5,2.25])
        plt.suptitle('Orbital Period vs. Planet Size', fontsize=20)

    plt.xlabel(r'Stellar Mass [M$_\odot$]', fontsize=14)
    plt.grid()


    plt.tight_layout(rect=[0, 0, 1, 0.92])


def flux_at_fixed_radius(R_choices=[1.4,1.6,1.8,2.0]):

    for j in range(len(R_choices)):

        R_fixed = R_choices[j]

        smassbins = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1]

        F_R_occ = np.zeros(len(smassbins))
        F_R_occ_std = np.zeros(len(smassbins))
        smass_mid = np.zeros(len(smassbins))

        for i in range(len(smassbins)):

            smass1, smass2 = smassbins[i], smassbins[i] + 0.3
            smass_mid[i] = ( smass1 + smass2 ) / 2

            # occurrence
            occ_array = np.loadtxt("./data/chain_occ_SincPradv2_{0}-{1}-smass.csv".format(smass1, smass2), delimiter=',')
            m_i, m_std_i = np.mean(occ_array[:,0]), np.std(occ_array[:,0])
            Rp100_i, Rp100_std_i = np.mean(occ_array[:,1]), np.std(occ_array[:,1])

            u = 1 / m_i

            # invert power law and calculate error
            lnF_i = log(100) + u*log(R_fixed / Rp100_i)
            lnF_err_i = ( u * u * log(Rp100_i / R_fixed) * m_std_i )**2 \
                      + ( u * (1 / Rp100_i) * Rp100_std_i )**2
            lnF_err_i = np.sqrt(lnF_err_i)

            # append
            F_R_occ[i] = lnF_i
            F_R_occ_std[i] = lnF_err_i


        plt.subplot(len(R_choices),1,j+1)
        plt.errorbar(smass_mid, F_R_occ, F_R_occ_std, capsize=2)
        # plt.plot(smass_mid, F_R_occ)
        plt.ylim([-10,10])
        plt.ylabel(r'log(F$_p$({} R$_\oplus$))'.format(R_fixed), fontsize=10)
        if j==len(R_choices)-1:
            plt.xlabel(r'Stellar Mass [M$_\odot$]', fontsize=14)
        plt.grid()


    plt.tight_layout(True)
