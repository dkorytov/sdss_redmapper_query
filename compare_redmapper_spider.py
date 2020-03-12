#!/usr/bin/env python3

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
from query import load_spider_fits, load_redmapper_cluster_fits, load_spiders_bcg_fits, combine_spiders_bcg
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
    m200m_richness_list = []
    m200c_richness_list = []
    for m200c, redshift, richness in zip(mass, spider_cat['z'], spider_cat['lambda']):
        m200m, r200m, c200m =mass_adv.changeMassDefinitionCModel(m200c, redshift, 
                                                                 "200c", "200m",
                                                                 c_model='child18')
        m200m_list.append(m200m)
        m200m, r200m = background.lambda_to_m200_r200(richness, redshift, richness_mass_author="Simet_mean")
        m200c, r200c = background.lambda_to_m200_r200(richness, redshift, richness_mass_author="Simet_crit")
        m200m_richness_list.append(m200m)
        m200c_richness_list.append(m200c)
    mass = np.array(m200m_list)
    r200m_r200c_ratio = r200m/r200
    rad *= r200m_r200c_ratio
    r200 *= r200m_r200c_ratio
    spider_cat['m200m'] = np.array(m200m_list)
    spider_cat['rad'] = rad
    spider_cat['m200m_lambda'] = np.array(m200m_richness_list)
    spider_cat['m200c_lambda'] = np.array(m200c_richness_list)

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
    slct = (redmapper_distance<0.24) & (distances < 0.1) & (spider_cat['z'] < 0.35) 
    print(np.sum(slct))
    print(np.sum(spider_cat['z']<0.35))
    redmapper_old = redmapper_cat
    redmapper_cat = select_dict(redmapper_cat, redmapper_indexes)
    mass_estimate_spider(spider_cat)
    mass_esitmate_redmapper(redmapper_cat)


    # plt.figure()
    # h, xbins = np.histogram(redmapper_distance[slct], bins=np.logspace(-3, 1, 32))
    # plt.semilogx(dtk.bins_avg(xbins), h, '-')

    plt.figure()
    h, xbins = np.histogram(distances[slct], bins=np.logspace(-7, 1, 32))
    plt.semilogx(dtk.bins_avg(xbins), h, '-')
    plt.xlabel('Separation [deg]')
    plt.ylabel('count')
    plt.tight_layout()
    
    plt.figure()
    plt.plot(redmapper_old['ra'], redmapper_old['dec'], 'x')
    plt.plot(spider_cat['ra'], spider_cat['dec'], '.r')
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.tight_layout()
    
    plt.figure()
    plt.plot(spider_cat['z'][slct], redmapper_cat['z'][slct], '.')
    plt.xlabel('SPIDERS redshift')
    plt.ylabel('redMaPPer redshift')
    plt.tight_layout()
    
    plt.figure()
    plt.plot(spider_cat['ra'][slct], spider_cat['dec'][slct], 'xb', alpha=0.3, label='SPIDERS')
    plt.plot(redmapper_cat['ra'][slct], redmapper_cat['dec'][slct], '.r', alpha=0.3, label='redmapper')
    plt.legend(loc='best')
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.tight_layout()
        
    plt.figure()
    plt.ylim([1e14, 3e15])
    plt.xlim([1e14, 3e15])
    plt.loglog(spider_cat['m200m'][slct], redmapper_cat['m200m'][slct], '.')
    plt.xlabel(r'SPIDERS M$_{200m}$ [M$_\odot$/h]')
    plt.ylabel('redMaPPer M$_{200m}$ [M$_\odot$/h]')
    # plt.axis('equal')
    plt.plot([8e13, 3e15], [8e13, 3e15], '--k', label='one-to-one')
    plt.legend(framealpha=0.0)

    plt.tight_layout()

    
    plt.figure()
    plt.loglog(spider_cat['m200c'][slct], redmapper_cat['m200c'][slct], '.')
    plt.ylim([8e13, 3e15])
    plt.xlim([8e13, 3e15])
    plt.xlabel(r'SPIDERS M$_{200c}$ [Msun/h]')
    plt.ylabel('redMaPPer M$_{200c}$ [Msun/h]')
    # plt.axis('equal')
    plt.plot([8e13, 3e15], [8e13, 3e15], '--k', label='one-to-one')
    plt.legend(framealpha=0.0)
    plt.tight_layout()
    

    plt.figure()
    plt.loglog(spider_cat['lambda'][slct], redmapper_cat['lambda'][slct], '.')
    plt.xlabel(r'SPIDERS Richness')
    plt.ylabel('redMaPPer Richness')
    plt.tight_layout()
    
    plt.figure()
    plt.loglog(spider_cat['m200m_lambda'][slct], redmapper_cat['m200m'][slct], '.')
    plt.xlabel(r'SPIDERS M200m from richness')
    plt.ylabel('redMaPPer M200m')
    plt.tight_layout()
    
    plt.figure()
    plt.semilogx(redmapper_cat['m200m'][slct], spider_cat['m200m_lambda'][slct]/redmapper_cat['m200m'][slct], '.')
    plt.ylabel(r'SPIDERS M200m from richness/RedMaPPer M200m')
    plt.xlabel('redMaPPer M200m')
    plt.tight_layout()
    
    plt.figure()
    plt.loglog(spider_cat['m200m'], spider_cat['m200m_lambda'], '.')
    plt.ylabel(r'SPIDERS M200m from richness')
    plt.xlabel('SPIDERS M200m from x-ray')
    plt.tight_layout()
    
    plt.figure()
    plt.plot(spider_cat['z'][slct], spider_cat['m200m'][slct]/redmapper_cat['m200m'][slct],  '.')
    plt.ylabel('SPIDER mass/redMaPPer mass')
    plt.xlabel('Cluster redshift')
    plt.tight_layout()

    plt.figure()
    plt.semilogx(redmapper_cat['m200m'][slct], spider_cat['lambda'][slct]/redmapper_cat['lambda'][slct], '.')
    plt.ylabel(r'SPIDERS richness/Redmapper Richness')
    plt.xlabel('RedMaPPer M200m')
    plt.tight_layout()
    
    plt.figure()
    plt.semilogx(redmapper_cat['lambda'][slct], redmapper_cat['m200m'][slct],  '.')
    plt.xlabel('RedMaPPer richness')
    plt.ylabel('RedMaPPer M200m')
    plt.tight_layout()
    # plt.ylabel(r'SPIDERS richness/Redmapper Richness')

    
    plt.figure()
    plt.semilogx(distances[slct], spider_cat['lambda'][slct]/redmapper_cat['lambda'][slct], '.')
    plt.ylabel(r'SPIDERS richness/Redmapper Richness')
    plt.xlabel('Center separation [deg]')
    plt.tight_layout()
    
    plt.figure()
    plt.semilogx(distances, spider_cat['lambda']/redmapper_cat['lambda'], '.')
    plt.ylabel(r'SPIDERS richness/Redmapper Richness')
    plt.xlabel('Center separation [deg]')
    plt.tight_layout()
    
    plt.figure()
    plt.semilogx(distances[slct], spider_cat['z'][slct]/redmapper_cat['z'][slct], '.')
    plt.ylabel(r'SPIDERS z/Redmapper z')
    plt.xlabel('Center separation [deg]')
    plt.tight_layout()
    
    plt.figure()
    plt.semilogx(distances[slct], spider_cat['z'][slct], '.')
    plt.ylabel(r'SPIDERS z')
    plt.xlabel('Center separation [deg]')
    plt.tight_layout()

    plt.figure()
    plt.semilogx(distances[slct], spider_cat['z'][slct] - redmapper_cat['z'][slct], '.')
    plt.ylabel(r'SPIDERS z -  Redmapper z')
    plt.xlabel('Center separation [deg]')
    plt.tight_layout()

def compare_redmapper_spider():
    spider_cat = load_spider_fits("/data/a/cpac/dkorytov/data/spiders/catCluster-SPIDERS_RASS_CLUS-v2.0.fits")
    spider_bcg_cat = load_spiders_bcg_fits("/data/a/cpac/dkorytov/data/spiders/SpidersXclusterBCGs-v2.0.fits")
    redmapper_cat = load_redmapper_cluster_fits("/data/a/cpac/dkorytov/data/redmapper/redmapper_dr8_public_v6.3_catalog.fits.gz")
    combine_spiders_bcg(spider_cat, spider_bcg_cat)
    # plt.figure()
    # plt.plot(spider_cat['ra_bcg'], spider_cat['dec_bcg'])
    # plt.show()
    spider_cat['ra'] = spider_cat['ra_bcg']
    spider_cat['dec'] = spider_cat['dec_bcg']
    print(spider_cat.keys())
    print(redmapper_cat.keys())
    # plt.figure()
    # plt.plot(redmapper_cat['ra'], redmapper_cat['dec'], '.', alpha=0.3)
    # plt.plot(spider_cat['ra'], spider_cat['dec'], '.r')
    
    find_match(spider_cat, redmapper_cat)


if __name__ == "__main__":
    background.set_cosmology('wmap7')
    compare_redmapper_spider()
    dtk.save_figs("figs/"+__file__+"/", extension='.pdf')
    plt.show()
