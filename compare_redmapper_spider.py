#!/usr/bin/env python2.7

from __future__ import print_function, division 
import numpy as np
import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from matplotlib import rc
from colossus.halo import mass_adv
from colossus.cosmology import cosmology
cosmology.setCosmology('WMAP7')
from query import load_spider_fits, load_redmapper_cluster_fits
import background
from astropy import units as u
from astropy.coordinates import SkyCoord


rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

def select_dict(dic, slct):
    new_dic = {}
    for key in dic.keys():
        new_dic[key]= dic[key][slct]
    return new_dic

def mass_estimate_spider(spider_cat):
    r200c_deg = spider_cat['r200c_deg']
    rad = r200c_deg * 60
    r200 = background.arcmin_to_r200(rad, spider_cat['z'])
    mass = background.r200c_to_m200c(r200, spider_cat['z'])
    spider_cat['m200c'] = mass
    m200m_list = []
    for m200c, redshift in zip(mass, spider_cat['z']):
        m200m, r200m, c200m =mass_adv.changeMassDefinitionCModel(m200c, redshift, 
                                                                 "200c", "200m",
                                                                 c_model='child18')
        m200m_list.append(m200m)
    mass = np.array(m200m_list)
    r200m_r200c_ratio = r200m/r200
    rad *= r200m_r200c_ratio
    r200 *= r200m_r200c_ratio
    spider_cat['m200m'] = np.array(m200m_list)
    spider_cat['rad'] = rad

def mass_esitmate_redmapper(redmapper_cat):
    richnesses = redmapper_cat['lambda']
    zs = redmapper_cat['z']
    m200m_list = []
    m200c_list = []
    for z, richness in zip(zs, richnesses):
        m200m, r200 = background.lambda_to_m200_r200(richness, z, richness_mass_author="Simet_mean")
        m200c, r200 = background.lambda_to_m200_r200(richness,z, richness_mass_author="Simet_crit")
        m200m_list.append(m200m)
        m200c_list.append(m200c)
    redmapper_cat['m200c'] = np.array(m200c_list)
    redmapper_cat['m200m'] = np.array(m200m_list)

def get_radec_mat(cat):
    return np.vstack((cat['ra'], cat['dec'], cat['z']*10)).T

def get_radec_separation(ra1, dec1, ra2, dec2):
    distance = np.zeros_like(ra1)
    for i in range(0, len(ra1)):
        c1 = SkyCoord(ra1[i]*u.deg, dec1[i]*u.deg, frame='icrs')
        c2 = SkyCoord(ra2[i]*u.deg, dec2[i]*u.deg, frame='icrs')
        distance[i] = c1.separation(c2).deg
    return distance
def find_match(spider_cat, redmapper_cat):
    spider_radec = get_radec_mat(spider_cat)
    redmapper_radec = get_radec_mat(redmapper_cat)
    ckdtree = cKDTree(redmapper_radec, compact_nodes=False, balanced_tree=False)
    redmapper_distance, redmapper_indexes = ckdtree.query(spider_radec)
    distances = get_radec_separation(spider_cat['ra'], spider_cat['dec'], redmapper_cat['ra'][redmapper_indexes], redmapper_cat['dec'][redmapper_indexes])
    slct = (redmapper_distance<0.24) & (distances < 0.01) & (spider_cat['z'] < 0.35) 

    redmapper_cat = select_dict(redmapper_cat, redmapper_indexes)
    mass_estimate_spider(spider_cat)
    mass_esitmate_redmapper(redmapper_cat)


    plt.figure()
    h, xbins = np.histogram(redmapper_distance[slct], bins=np.logspace(-3, 1, 32))
    plt.semilogx(dtk.bins_avg(xbins), h, '-')


    plt.figure()
    h, xbins = np.histogram(distances[slct], bins=np.logspace(-3, 1, 32))
    plt.semilogx(dtk.bins_avg(xbins), h, '-')
    plt.xlabel('Separation [deg]')
    plt.ylabel('count')
    

    plt.figure()
    plt.plot(spider_cat['z'][slct], redmapper_cat['z'][slct], '.')
    plt.xlabel('SPIDERS redshift')
    plt.ylabel('redMaPPer redshift')

    plt.figure()
    plt.plot(spider_cat['ra'][slct], spider_cat['dec'][slct], 'xb', alpha=0.3, label='SPIDERS')
    plt.plot(redmapper_cat['ra'][slct], redmapper_cat['dec'][slct], '.r', alpha=0.3, label='redmapper')
    plt.legend(loc='best')
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.figure()
    plt.loglog(spider_cat['m200m'][slct], redmapper_cat['m200m'][slct], '.')
    plt.xlabel(r'SPIDERS M$_{200m}$ [Msun/h]')
    plt.ylabel('redMaPPer M$_{200m}$ [Msun/h]')

    plt.figure()
    plt.loglog(spider_cat['m200c'][slct], redmapper_cat['m200c'][slct], '.')
    plt.ylim([1e14, 1e16])
    plt.xlabel(r'SPIDERS M$_{200c}$ [Msun/h]')
    plt.ylabel('redMaPPer M$_{200c}$ [Msun/h]')
    
def compare_redmapper_spider():
    spider_cat = load_spider_fits("/media/luna1/dkorytov/data/spider_xray/catCluster-SPIDERS_RASS_CLUS-v2.0.fits")
    redmapper_cat = load_redmapper_cluster_fits("redmapper_dr8_public_v6.3_catalog.fits")
    print(spider_cat.keys())
    print(redmapper_cat.keys())
    # plt.figure()
    # plt.plot(redmapper_cat['ra'], redmapper_cat['dec'], '.', alpha=0.3)
    # plt.plot(spider_cat['ra'], spider_cat['dec'], '.r')
    
    find_match(spider_cat, redmapper_cat)


if __name__ == "__main__":
    background.set_cosmology('wmap7')
    compare_redmapper_spider()
    dtk.save_figs("figs/"+__file__+"/")
    plt.show()
