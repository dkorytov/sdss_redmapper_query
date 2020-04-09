#!/usr/bin/env python3

from __future__ import print_function, division


from matplotlib import rc

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'family':'serif','serif':['Palatino']})
# rc('font',**{'size':
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)
# rc('legend', fontsize=13)

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk

from background import *

def load_zmr(param_fname):
    param = dtk.Param(param_fname)
    result_folder     = param.get_string("result_folder")
    
    fname = result_folder+"type1_weight1_mag1_clr1_result.hdf5"
    print(fname)
    zmr = h5py.File(fname, 'r')
    return zmr

def plot_zmr(param_fname):
    zmr = load_zmr(param_fname)
    mass_bins = zmr['m_bins'].value
    r_bins = zmr['r_bins'].value
    r_bins_center = dtk.bins_avg(r_bins)
    plt.figure()
    colors = ['b', 'g', 'r', 'c', 'm']
    for i in range(0, len(mass_bins)-1):
        if zmr['zm_counts'][0][i]>2:
            rad_profile = zmr['zmr_gal_density'].value[0, i]
            rad_profile_err = zmr['zmr_gal_density_err'].value[0, i]
            label = "--"#"{:.2f} < log$_{{10}}$ M$_{{200c}}$ < {:.2f}".format(np.log(mass_bins[i]), np.log(mass_bins[i+1]))
            label = "{:.2f} $<$ log$_{{10}}$ M$_{{200c}}$ $<$ {:.2f}".format(np.log10(mass_bins[i]), np.log10(mass_bins[i+1]))
            plt.errorbar(r_bins_center, rad_profile,
                         yerr=rad_profile_err, fmt='-o', markeredgecolor='none',
                         label=label, markerfacecolor=colors[i], color=colors[i],
                         capsize=0)
            yerr_min = np.clip(rad_profile-rad_profile_err, a_max=np.inf, a_min=1e-3)
            plt.fill_between(r_bins_center,
                             rad_profile+rad_profile_err, yerr_min,
                             color=colors[i], alpha=0.3)
    plt.legend(loc='best', framealpha=0.0)
    plt.yscale('log')
    plt.xlabel('r/R$_{200c}$')
    plt.ylabel('$\Sigma_{galaxy} [(R_{200c})^{-2}$]')
    dtk.save_figs("figs/"+param_fname+"/"+__file__+"/")

def check_error_zmr(param_fname):
    zmr = load_zmr(param_fname)
    mass_bins = zmr['m_bins'].value
    r_bins = zmr['r_bins'].value
    r_bins_center = dtk.bins_avg(r_bins)
    for i in range(0, len(mass_bins)-1):
        if zmr['zm_counts'][0][i] > 1:
            print( "\n{:.2e} < M200 < {:.2e} ".format(mass_bins[i], mass_bins[i+1]))
            print("counts: ", zmr['zmr_counts'][0][i][:])
            print("density:", zmr['zmr_gal_density'][0][i][:])
            print("error:  ", zmr['zmr_gal_density_err'][0][i][:])
    
if __name__ == "__main__":
    param_file = sys.argv[1]
    # plot_zmr(param_file)
    check_error_zmr(param_file)
    plt.show()
