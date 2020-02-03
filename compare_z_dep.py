#!/usr/bin/env python3

import numpy as np
import matplotlib
import os
# checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')




import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.cosmology import WMAP9 as cosmo
from background import *
import dtk
import sys
import h5py

def compare_z_dep(param_fname):
    param = dtk.Param(param_fname)
    result_folder = param.get_string('result_folder')
    hfile = h5py.File(result_folder+'type1_weight1_mag1_clr1_result.hdf5', 'r')
    
    zmr = hfile['zmr_gal_density'][()]
    zmr_err = hfile['zmr_gal_density_err'][()]
    zmr_cnt = hfile['zmr_counts'][()]
    r_bins  = hfile['r_bins'][()]
    r_bins_cen = dtk.bins_avg(r_bins)
    z_bins = hfile['z_bins'][()]
    print(r_bins_cen)

    shape = zmr.shape
    for i in range(shape[1]):
        if not np.isfinite(np.sum(zmr[0,i,:])):
            continue
        plt.figure()
        for j in range(shape[0]):
            plt.plot(r_bins_cen, zmr[j,i,:], label='{:.2f}<z<{:.2f}'.format(z_bins[j], z_bins[j+1]))
            plt.fill_between(r_bins_cen, zmr[j,i,:]+zmr_err[j,i,:], zmr[j,i,:]-zmr_err[j,i,:], alpha=0.3)

        print(np.sum(zmr_cnt[0,i,:]), np.sum(zmr_cnt[1,i,:]), np.nanmean(np.sum(zmr_cnt[0,i,:])/np.sum(zmr_cnt[1,i,:])))
        print(zmr[0,i,:])
        print(zmr[1,i,:])
        print(np.nanmean(zmr[0,i,:]/zmr[1,i,:]))
        print('\n')
        plt.yscale('log')
        plt.legend()
        np.log10(np.nanmean(zmr[0,i,:]/zmr[1,i,:]))
        plt.show()



if __name__ == '__main__':
    compare_z_dep(sys.argv[1])
