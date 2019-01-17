#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr

from astropy.cosmology import WMAP9 as cosmowmap9
from astropy.cosmology import Planck13 as cosmoplanck15
from astropy.cosmology import WMAP7 as cosmowmap7

from astropy.io import fits as pyfits
import h5py

import dtk
import background 



if __name__ == "__main__":
    richness = np.logspace(np.log10(20), np.log10(300), 100)
    m200m_rykoff = background.lambda_to_m200m_Rykoff( richness, 0.25)
    m200c_rykoff = background.lambda_to_m200c_Rykoff( richness, 0.25)
    m200m_simet  = background.lambda_to_m200m_Simet(  richness, 0.25)
    m200m_baxter = background.lambda_to_m200m_Baxter( richness, 0.25)
    plt.figure()
    plt.plot(richness, m200m_rykoff, label='Rykoff et al, 2012')
    plt.plot(richness, m200m_simet, label='Simet et al, 2016')
    plt.plot(richness, m200m_baxter, label='Baxter et al, 2016')
    plt.legend(loc='best', framealpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'M$_{200m}$ [h$^{-1}$ M$_\odot$]')
    plt.grid()

    plt.figure()
    plt.plot(richness, m200m_rykoff, label='M200_m, Rykoff et al, 2012')
    plt.plot(richness, m200c_rykoff, label='M200_c, Rykoff et al, 2012')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'M$_{200m}$ [h$^{-1}$ M$_\odot$]')
    plt.grid()

    plt.figure()
    plt.plot(richness, m200c_rykoff/m200m_rykoff, '-', label='M200_c/M200_m')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'M$_{200m}$ [h$^{-1}$ M$_\odot$]')
    plt.grid()
    
    #header list can be found in http://arxiv.org/pdf/1303.3562v2.pdf
    hdulist = pyfits.open("redmapper_dr8_public_v6.3_catalog.fits")
    hdurnd  = pyfits.open("redmapper_dr8_public_v6.3_randoms.fits")
    tdata =  hdulist[1].data
    rdata =  hdurnd[1].data
    # red_ra = tdata.field('ra')
    # red_dec= tdata.field('dec')
    red_z  = tdata.field('z_lambda')
    red_lambda=tdata.field('lambda')
    # red_pcen=tdata.field('p_cen')

    m200m_rykoff = background.lambda_to_m200m_Rykoff( red_lambda, 0.25)
    m200m_simet  = background.lambda_to_m200m_Simet(  red_lambda, 0.25)
    m200m_baxter = background.lambda_to_m200m_Baxter( red_lambda, 0.25)
    
    xbins = np.logspace(13.5,15.5,100)
    xbins_avg = dtk.bins_avg(xbins)
    
    plt.figure()
    h, _ = np.histogram(m200m_rykoff, bins = xbins)
    plt.plot(xbins_avg, h, label='rykoff')

    h, _ = np.histogram(m200m_simet, bins = xbins)
    plt.plot(xbins_avg, h, label='simet')

    h, _ = np.histogram(m200m_baxter, bins = xbins)
    plt.plot(xbins_avg, h, label='baxter')

    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Counts')
    plt.xlabel('M$_{200m}$ [h$^{-1}$ M$_\odot$]')
    plt.legend(loc='best', framealpha=1.0)
    plt.grid()
    plt.show()
    

    
