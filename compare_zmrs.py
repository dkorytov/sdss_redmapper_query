#!/usr/bin/env python2.7

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



def get_zmr(param_fname):
    zmr = {}
    param = dtk.Param(param_fname)
    result_fname =param.get_string('result_folder')+"type1_weight1_mag1_clr1_result.hdf5"
    hfile = h5py.File(result_fname, 'r')
    zmr['z_bins'] = hfile['z_bins'].value
    zmr['m_bins'] = hfile['m_bins'].value
    zmr['r_bins'] = hfile['r_bins'].value
    zmr['zm_counts'] = hfile['zm_counts'].value
    zmr['zmr_gal_density'] = hfile['zmr_gal_density'].value
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

    colors = ['b', 'g', 'r', 'c', 'm', 'y']



    for m_i in range(0, len(zmr1['m_bins_cen'])):
        if m_i > 5:
            break
        xi = m_i//2
        yi = m_i%2
        axs[2*xi, yi].plot(zmr1['r_bins_cen'], zmr1['zmr_gal_density'][0,m_i,:], color=colors[m_i%len(colors)], label=zmr1['label']+"[{:.0f}]".format(zmr1['zm_counts'][0,m_i]))
        axs[2*xi, yi].plot(zmr2['r_bins_cen'], zmr2['zmr_gal_density'][0,m_i,:], color=colors[m_i%len(colors)], label=zmr2['label']+"[{:.0f}]".format(zmr2['zm_counts'][0,m_i]), ls='--')
        axs[2*xi, yi].set_ylim([1e-1,1e3])
        axs[2*xi, yi].set_yscale('log')
        axs[2*xi, yi].legend(loc='best')
        axs[2*xi, yi].set_title("{:.2f}<M200<{:.2f}".format(np.log10(zmr1['m_bins'][m_i]), np.log10(zmr1['m_bins'][m_i+1])))
        
        relative_err = (zmr1['zmr_gal_density'][0,m_i,:]-zmr2['zmr_gal_density'][0,m_i,:])/((zmr1['zmr_gal_density'][0,m_i,:] + zmr2['zmr_gal_density'][0,m_i,:])/2.0)
        axs[2*xi+1, yi].plot(zmr1['r_bins_cen'], relative_err, color=colors[m_i])
        axs[2*xi+1, yi].axhline(0, ls='--', color='k')
        axs[2*xi+1, yi].set_ylim([-1,+1])

    f.tight_layout()

def get_label(param_fname):
    if "spider" in param_fname:
        if "richness" in param_fname:
            return "SPIDER richness"
        elif "xray" in param_fname:
            return "SPIDER x-ray"
    else:
        return "Redmapper"

def compare_zmrs(param_fname1, param_fname2):
    zmr1 = get_zmr(param_fname1)
    zmr2 = get_zmr(param_fname2)
    plot_zmr(zmr1)
    plot_zmr(zmr2)
    plot_zmrs(zmr1, zmr2)
    plt.show()

if __name__ == "__main__":
    compare_zmrs(sys.argv[1], sys.argv[2])
