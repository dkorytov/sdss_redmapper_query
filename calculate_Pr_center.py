#!/usr/bin/env python3

from __future__ import print_function, division 
import numpy as np
import matplotlib
import os
import h5py
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr
from astropy.io import fits as pyfits
from astropy.cosmology import WMAP7 as cosmowmap7
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)



def calculate_Pr_center(fname):
    cat = {}
    hdulist = pyfits.open(fname)
    tdata = hdulist[1].data
    print(hdulist[1].columns)
    # help(hdulist[1])
    print(tdata)
    cat['id'] = tdata.field('id')
    cat['ra'] = tdata.field('ra')
    cat['dec']= tdata.field('dec')
    cat['pcen']  = tdata.field('p')
    cat['pfree'] = tdata.field('p_free')
    cat['r']  = tdata.field('r')
    print(cat['id'])
    for cat_id in np.unique(cat['id']):
        print(cat_id)
        slct = cat['id'] == cat_id
        plt.figure()

        plt.scatter(cat['ra'][slct], cat['dec'][slct], s=cat['pcen'][slct]*10)
        plt.plot(cat['ra'][slct], cat['dec'][slct], 'xr')

        print(cat['r'][slct])
        print(cat['id'][slct])
        print(cat['pcen'][slct])
        print(cat['pfree'][slct])
        plt.show()
        
def deg2rad(deg):
    return (deg/360.0)*(2.0*np.pi)

def rad2deg(rad):
    return (rad/(2.0*np.pi))*360.0

def rad_dist2(ra1, dec1, ra2, dec2):
    ra1=deg2rad(ra1)
    ra2=deg2rad(ra2)
    dec1=deg2rad(dec1)
    dec2=deg2rad(dec2)
    del_dec = dec1-dec2
    del_ra  = ra1-ra2
    radial_dist= 2.0*np.arcsin(np.sqrt(np.sin(del_dec/2.0)**2+np.cos(dec1)*np.cos(dec2)*np.sin(del_ra/2.0)**2))
    return rad2deg(radial_dist)

def select_prob(probs):
    if ~np.isclose(np.sum(probs),  1.0):
        print("Probabilities do not add up to 1", np.sum(probs))
        raise
    rand = np.random.rand()
    probs_cum = np.cumsum(probs)
    i =0
    while(i<len(probs) and probs_cum[i] < rand):
        i+=1
    return i

def select_radius(probs, radii):
    i = select_prob(probs)
    return radii[i]

def save_radii(fname, radii_to_save):
    hfile = h5py.File(fname, 'w')
    hfile['radii'] = np.array(radii_to_save, dtype='float')
    hfile.close()

def calculate_Pr_center2(fname):
    cat = {}
    hdulist = pyfits.open(fname)
    tdata = hdulist[1].data
    print(hdulist[1].columns)
    # help(hdulist[1])
    print(tdata)
    cat['p_cen'] = tdata.field('p_cen')
    cat['z_lambda'] = tdata.field('z_lambda')
    cat['ra_cen'] = tdata.field('ra_cen')
    cat['dec_cen'] = tdata.field('dec_cen')
    print(cat['p_cen'])
    probs = []
    d_kpcs = []
    non_center_avg = []
    non_center_avg_prob = []
    correct_prop = []
    radii_to_save = []
    print("progress")
    for i in range(len(cat['z_lambda'])):
    # for i in range(1000):
        z = cat['z_lambda'][i]
        a = 1.0/(1.0+z)
        prob = cat['p_cen'][i]
        ra = cat['ra_cen'][i]
        dec = cat['dec_cen'][i]
        d_r = rad_dist2(ra, dec, ra[0], dec[0])
        d_kpc = cosmowmap7.kpc_proper_per_arcmin(z).value * (60.0*d_r)
        probs += list(prob)
        d_kpcs += list(d_kpc)
        ## Converting from physical kpc to comoving kpc/h
        ## 1.27323
        non_center_avg.append(np.average(d_kpcs[1:])*0.7/a)
        non_center_avg_prob.append(np.average(d_kpc[1:], weights=prob[1:])*0.7/a)
        correct_prop.append(prob[0])
        radii_to_save.append(select_radius(prob, d_kpc)*0.7/a)
        if i%1000 == 0:
            print('\t{:.3f}\r'.format(i/len(cat['z_lambda'])), end='')
    print()
    save_radii('data/Pr_center/distribution.hdf5', np.array(radii_to_save)/1000)
    save_radii('data/Pr_center/distribution.0mpc.hdf5', np.zeros_like(radii_to_save))
    save_radii('data/Pr_center/distribution.1mpc.hdf5', np.ones_like(radii_to_save))
    probs = np.array(probs)
    d_kpcs = np.array(d_kpcs)
    
    # probs[d_kpcs>0.2] = 0
    plt.figure()
    plt.hist(probs, bins=np.linspace(0,1,100))
    plt.yscale('log')
    
    plt.figure()
    plt.hist(d_kpcs/1000.0, bins=64, histtype='step', label='unweighted')
    plt.hist(d_kpcs/1000.0, bins=64, weights=probs, histtype='step', label='weighted')
    plt.ylabel('Counts')
    plt.xlabel('Clustering Miscentering Distance\n[h$^{-1}$ Mpc Comoving]')
    plt.legend(framealpha=0.0)
    # plt.yscale('log')
    plt.tight_layout()

    num_clusters = len(cat['z_lambda'])
    avg_displaced = np.sum((probs*d_kpcs))/num_clusters
    print("Average Miscentering: {:2f} [comoving kpc/h]".format(avg_displaced))
    print("Average Distance when Miscentered: {:.2f} [comoving kpc/h]".format(np.average(non_center_avg_prob)))
    print("Miscentering Prob: {:.2f}".format(np.average(correct_prop)))

def test_saved_data(fname='data/Pr_center/distribution1.hdf5'):
    hfile = h5py.File(fname,'r')
    data = hfile['radii'][:]
    plt.figure()
    plt.hist(data/1000, bins=100, histtype='step')
    plt.ylabel('Counts')
    plt.xlabel('Distance [Mpc Physical]')
    plt.yscale('log')

    
if __name__ == "__main__":
    # fname = 'redmapper_dr8_public_v6.3_members.fits.gz'
    # calculate_Pr_center(fname)
    fname = 'redmapper_dr8_public_v6.3_catalog.fits.gz'
    calculate_Pr_center2(fname)
    test_saved_data()
    dtk.save_figs('figs/'+__file__+'/', '.pdf')
    dtk.save_figs('figs/'+__file__+'/', '.png')
    plt.show()
