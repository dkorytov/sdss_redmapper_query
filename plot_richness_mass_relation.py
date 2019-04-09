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
from halotools.empirical_models import NFWProfile
from colossus.halo.mass_so import M_to_R, R_to_M

def plot_richness_mass_relation():
    background.set_cosmology("wmap7")
    richness = np.logspace(np.log10(20), np.log10(300), 100)
    m200m_rykoff, r200m_rykoff = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Rykoff_mean")
    m200m_simet,  r200m_simet  = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Simet_mean")
    m200m_baxter, r200m_baxter = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Baxter_mean")
    m200c_rykoff, r200c_rykoff = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Rykoff_crit")
    m200c_simet,  r200c_simet  = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Simet_crit")
    m200c_baxter, r200c_baxter = background.lambda_to_m200_r200( richness, 0.25, richness_mass_author = "Baxter_crit")
    plt.figure()
    plt.plot(richness, m200m_rykoff,'-r')
    plt.plot(richness, m200m_simet, '-b')
    plt.plot(richness, m200m_baxter,'-g')
    plt.plot(richness, m200c_rykoff,'--r')
    plt.plot(richness, m200c_simet, '--b')
    plt.plot(richness, m200c_baxter,'--g')
    plt.plot([],[],'r',label='rykoff')
    plt.plot([],[],'b',label='simet')
    plt.plot([],[],'g',label='baxter')
    plt.plot([],[],'-k',label='M200m')
    plt.plot([],[],'--k',label='M200c')
    plt.legend(loc='best', framealpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'M$_{200}$')
    plt.grid()

    plt.figure()
    plt.plot(richness, r200m_rykoff,'-r')
    plt.plot(richness, r200m_simet, '-b')
    plt.plot(richness, r200m_baxter,'-g')
    plt.plot(richness, r200c_rykoff,'--r')
    plt.plot(richness, r200c_simet, '--b')
    plt.plot(richness, r200c_baxter,'--g')
    plt.plot([],[],'r',label='rykoff')
    plt.plot([],[],'b',label='simet')
    plt.plot([],[],'g',label='baxter')
    plt.plot([],[],'-k',label='R200m')
    plt.plot([],[],'--k',label='R200c')
    plt.legend(loc='best', framealpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'R$_{200}$')
    plt.grid()


    plt.figure()
    plt.plot(richness, m200c_rykoff/m200m_rykoff, '-r', label='rykoff')
    plt.plot(richness, m200c_simet/m200m_simet, '-b', label='simet')
    plt.plot(richness, m200c_baxter/m200m_baxter, '-g', label='baxter')
    plt.legend(loc='best', framealpha=0.3)
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'M200c/M200m')
    plt.grid()

    plt.figure()
    plt.plot(richness, r200c_rykoff/r200m_rykoff, '-r', label='rykoff')
    plt.plot(richness, r200c_simet/r200m_simet, '-b', label='simet')
    plt.plot(richness, r200c_baxter/r200m_baxter, '-g', label='baxter')
    plt.legend(loc='best', framealpha=0.3)
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel(r'R200c/R200m')
    plt.grid()


    plt.figure()
    a = 1.0/(1.0+0.25)
    # plt.plot(m200m_simet, r200m_simet, '-b',  label='my simet mean')
    # plt.plot(m200c_simet, r200c_simet, '--b', label='my simet crit')
    plt.plot(m200m_simet, r200m_simet, '-b',  label='my simet mean')
    plt.plot(m200c_simet, r200c_simet, '--b', label='my simet crit')

    nfw_200m = NFWProfile(mdef='200m', redshift=0.25)
    nfw_200c = NFWProfile(mdef='200c', redshift=0.25)
    plt.plot(m200m_simet, nfw_200m.halo_mass_to_halo_radius(m200m_simet)*1000, '-r',  label='ht simet mean')
    plt.plot(m200c_simet, nfw_200c.halo_mass_to_halo_radius(m200c_simet)*1000, '--r', label='ht simet crit')

    plt.plot(m200m_simet, M_to_R(m200m_simet, 0.25, '200m'), '-g', label='cls simet mean')
    plt.plot(m200c_simet, M_to_R(m200c_simet, 0.25, '200c'), '--g', label='cls simet crit')
    
    plt.legend(loc='best', framealpha=0.3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('M200 [Msun]')
    plt.ylabel('R200 [proper kpc]')
    plt.grid()

    plt.show()
    
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

    plt.figure()
    h, xbins = np.histogram(red_lambda, bins=256)
    plt.plot(dtk.bins_avg(xbins), h)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Richness')
    plt.ylabel('count')

    m200m_rykoff = background.lambda_to_m200m_Rykoff( red_lambda, 0.25)
    m200m_simet  = background.lambda_to_m200m_Simet(  red_lambda, 0.25)
    m200m_baxter = background.lambda_to_m200m_Baxter( red_lambda, 0.25)
    m200c_rykoff = background.lambda_to_m200c_Rykoff( red_lambda, 0.25)
    m200c_simet  = background.lambda_to_m200c_Simet(  red_lambda, 0.25)
    m200c_baxter = background.lambda_to_m200c_Baxter( red_lambda, 0.25)
    
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

    plt.figure()
    h, _ = np.histogram(m200c_rykoff, bins = xbins)
    plt.plot(xbins_avg, h, label='rykoff')
    h, _ = np.histogram(m200c_simet, bins = xbins)
    plt.plot(xbins_avg, h, label='simet')
    h, _ = np.histogram(m200c_baxter, bins = xbins)
    plt.plot(xbins_avg, h, label='baxter')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Counts')
    plt.xlabel('M$_{200c}$ [h$^{-1}$ M$_\odot$]')
    plt.legend(loc='best', framealpha=1.0)
    plt.grid()

    plt.show()
    

    
if __name__ == "__main__":
    plot_richness_mass_relation()
