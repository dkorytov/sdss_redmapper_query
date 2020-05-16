#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import h5py
import dtk

import background

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

def plot_mstar_background(param_fname):    
    param = dtk.Param(param_fname)
    background_folder=param.get_string('background_folder')
    galaxy_type = param.get_int('galaxy_type')
    galaxy_weight = param.get_int('galaxy_weight')
    background.set_gal_type_weight(galaxy_type, galaxy_weight)
    background.set_gal_clr_type(param.get_int('galaxy_color_type'), param.get_int_list('galaxy_mag_lim1'))
    background.set_gal_mag_type(param.get_int('galaxy_mag_type'), param.get_int_list('galaxy_mag_lim1'))
    plt.figure()
    labels = {-1:'$>2.50 $L$_\star$',
              -0.5:'$>1.58 $L$_\star$',
              0:'$> 1.00 $L$_\star$',
              0.5:'$> 0.63 $L$_\star$',
              1:'$> 0.40 $L$_\star$',
        }
    for mstar_cut in [-1, -0.5, 0, 0.5, 1.0][::-1]:
        background.set_mstar_cut(mstar_cut)
        a, b  = background.get_background_estimate_sqdeg(background_folder, galaxy_type, galaxy_weight)
        z_range = np.linspace(0.1, 0.35, 100)
        plt.plot(z_range, a(z_range), label=labels[mstar_cut])
        plt.yscale('log')
        plt.xlabel('redshift')
        plt.ylabel('background galaxy density [deg$^{-2}$]')
        plt.legend(framealpha=0.0, labelspacing=0.0)
        plt.tight_layout()


if __name__ == "__main__":
    param_fname = sys.argv[1]
    plot_mstar_background(param_fname)
    dtk.save_figs('figs/'+__file__+'/', '.pdf')
    plt.show()
