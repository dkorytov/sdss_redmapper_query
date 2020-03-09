#!/usr/bin/env python3

from __future__ import print_function, division
import numpy as np
import matplotlib
import os
# checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import dtk
import sys

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)


def get_zmr(param_fname):
    zmr = {}
    param = dtk.Param(param_fname)
    result_fname =param.get_string('result_folder')+"type1_weight1_mag1_clr1_result.hdf5"
    hfile = h5py.File(result_fname, 'r')
    print(hfile.keys())
    zmr['z_bins'] = hfile['z_bins'].value
    zmr['m_bins'] = hfile['m_bins'].value
    zmr['r_bins'] = hfile['r_bins'].value
    zmr['zm_counts'] = hfile['zm_counts'].value
    zmr['zmr_gal_density'] = hfile['zmr_gal_density'].value
    zmr['zmr_gal_density_err'] = hfile['zmr_gal_density_err'].value
    zmr['z_bins_cen'] = dtk.bins_avg(zmr['z_bins'])
    zmr['m_bins_cen'] = dtk.bins_avg(zmr['m_bins'])
    zmr['r_bins_cen'] = dtk.bins_avg(zmr['r_bins'])
    zmr['label'] = get_label(param_fname)
    return zmr

def plot_zmr(zmr):
    plt.figure()
    for m_i in range(0, len(zmr['m_bins_cen'])):
        plt.plot(zmr['r_bins_cen'], zmr['zmr_gal_density'][0,m_i,:])
    plt.yscale('log')

def plot_zmrs(zmr1, zmr2):
    f, axs = plt.subplots(6,2, figsize=(10,15))

    colors = ['tab:red', 'tab:blue', 'r', 'c', 'm', 'y']
    for m_i in range(0, len(zmr1['m_bins_cen'])):
        if m_i > 5:
            break
        xi = m_i//2
        yi = m_i%2
        label = zmr1['label']+"[{:.0f}]".format(zmr1['zm_counts'][0,m_i])
        axs[2*xi, yi].plot(zmr1['r_bins_cen'], zmr1['zmr_gal_density'][0,m_i,:], color=colors[m_i%len(colors)], label=label)
        label = zmr2['label']+"[{:.0f}]".format(zmr2['zm_counts'][0,m_i])
        axs[2*xi, yi].plot(zmr2['r_bins_cen'], zmr2['zmr_gal_density'][0,m_i,:], color=colors[m_i%len(colors)], label=label, ls='--')
        axs[2*xi, yi].set_ylim([1e-1,1e3])
        axs[2*xi, yi].set_yscale('log')
        axs[2*xi, yi].legend(loc='best')
        axs[2*xi, yi].set_title("{:.2f}<M200<{:.2f}".format(np.log10(zmr1['m_bins'][m_i]), np.log10(zmr1['m_bins'][m_i+1])))
        
        relative_err = (zmr1['zmr_gal_density'][0,m_i,:]-zmr2['zmr_gal_density'][0,m_i,:])/((zmr1['zmr_gal_density'][0,m_i,:] + zmr2['zmr_gal_density'][0,m_i,:])/2.0)
        axs[2*xi+1, yi].plot(zmr1['r_bins_cen'], relative_err, color=colors[m_i])
        axs[2*xi+1, yi].axhline(0, ls='--', color='k')
        axs[2*xi+1, yi].set_ylim([-1,+1])
    f.tight_layout()


def plot_zmrs2(zmr1, zmr2):
    nrow = 3
    ncol = 2
    f, axs = plt.subplots(ncol, nrow, figsize=(15, 6), sharex=True, sharey=True)
    colors=['tab:red', 'tab:blue', 'g', 'c', 'm', 'y'] 
    f2, axs2 = plt.subplots(ncol, nrow, figsize=(15, 4), sharex=True, sharey=True)

    
    for m_index in range(0, 6):
        x_index = m_index//nrow
        y_index = m_index%nrow
        print(x_index, y_index)
        ax = axs[x_index, y_index]
        ax2 = axs2[x_index, y_index]
        for zmr, color in zip([zmr1, zmr2], ['tab:red', 'tab:blue']):
            y1 = zmr['zmr_gal_density'][0,m_index,:]
            yerr1 = zmr['zmr_gal_density_err'][0,m_index,:]
            y2 = zmr1['zmr_gal_density'][0,m_index,:]
            yerr2 = zmr1['zmr_gal_density_err'][0,m_index,:]
            yerr1_min = np.clip(y1-yerr1, a_min=1e-3, a_max=np.inf)
            ax.plot(zmr['r_bins_cen'], y1, color=color)
            ax.fill_between(zmr['r_bins_cen'], y1+yerr1, yerr1_min, color=color, alpha=0.3)
            #relative difference
            ax2.plot(zmr['r_bins_cen'], y1/y2, color=color)
            ax2.fill_between(zmr['r_bins_cen'], y1/y2-yerr1/y2, y1/y2+yerr1/y2, color=color, alpha=0.3)
        ax.set_yscale('log')
        mass_label="{:.2f}$<$M$_{{200m}}<${:.2f}".format(np.log10(zmr1['m_bins'][m_index]), np.log10(zmr1['m_bins'][m_index+1]))
        ax.text(0.1, 0.1, mass_label, transform=ax.transAxes)
        ax2.text(0.1, 0.1, mass_label, transform=ax2.transAxes)
        ax.set_ylim(1e-1, 1e3)
        ax2.set_xlim(0,1)
        ax2.set_ylim(0,2)

        if y_index == 0:
            ax.set_ylabel("$\Sigma_{gal}$ [R$_{200m}^{-2}$]")
            ax2.set_ylabel("$\Sigma_{gal}/\Sigma_{gal}^{redMaPPer}$")
        if x_index == 1:
            ax.set_xlabel("r/R$_{200m}$")
            ax2.set_xlabel("r/R$_{200m}$")
        if x_index == 0 and y_index== 0:
            ax.plot([],[], 'tab:red', lw=2, label='RedMaPPer')
            ax.legend(loc='upper right', framealpha=0.0)            

        if x_index == 0 and y_index== 1:
            ax.plot([],[], 'tab:blue', lw=2, label='SPIDERS')
            ax.legend(loc='upper right', framealpha=0.0)
        # 2nd plot comparing ratios
        # y = zmr2['zmr_gal_density'][0,m_index,:]/zmr1['zmr_gal_density'][0,m_index,:]
        # y_err = zmr2['zmr_gal_density_err'][0,m_index,:]/zmr1['zmr_gal_density'][0,m_index,:]
        # y1_err = zmr2['zmr_gal_density_err'][0,m_index,:]/zmr1['zmr_gal_density'][0,m_index,:]
        # ax2.plot(zmr1['r_bins_cen'], y, color=colors[m_index])
        # ax2.fill_between(zmr1['r_bins_cen'], y-y_err, y+y_err, alpha=0.3, color=colors[m_index])
        # ax2.fill_between(zmr1['r_bins_cen'], 1-y1_err, 1-y_err, alpha=0.3, color=colors[m_index])
        # ax2.set_xlabel('r/R$_{200m}$')
        # ax2.set_ylabel('$\Sigma^{SPIDERS}_{gal}/\Sigma^{redMaPPer}_{gal}$')
    
    f.tight_layout()
    f2.tight_layout()


def get_label(param_fname):
    if "spider" in param_fname:
        if "richness" in param_fname:
            return "SPIDER richness"
        elif "xray" in param_fname:
            return "SPIDER x-ray"
        else:
            return "SPIDER X-ray"
    else:
        return "Redmapper"

def compare_zmrs(param_fname1, param_fname2):
    zmr1 = get_zmr(param_fname1)
    zmr2 = get_zmr(param_fname2)
    plot_zmr(zmr1)
    plot_zmr(zmr2)
    plot_zmrs(zmr1, zmr2)
    plot_zmrs2(zmr1, zmr2)
    

if __name__ == "__main__":
    if len(sys.argv) > 2:
        compare_zmrs(sys.argv[1], sys.argv[2])
        dtk.save_figs("figs/"+__file__+"/", extension='.pdf')
    else:
        arg1 = 'params/rad_profile/mstar0_wmap7_simet_mean4.param'
        # arg2 = 'params/rad_profile/mstar0_wmap7_spider_mean_bcg.param'
        arg2 = 'params/rad_profile/spider/mstar0_wmap7_spider_crit_xray_bcg.param'
        compare_zmrs(arg1, arg2)
        dtk.save_figs("figs"+__file__+"/redmapper_vs_spiders/", extension='.pdf')
    plt.show()
