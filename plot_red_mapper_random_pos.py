#!/usr/bin/env python2.7
import sys
import numpy as np
import matplotlib.pyplot as plt
import sqlcl
from astropy.io import fits as pyfits
import dtk
from StringIO import StringIO
import background
import time
import h5py 

param = dtk.Param(sys.argv[1])
query_data_folder = param.get_string("query_data_folder")
cluster_size_max  = param.get_bool("cluster_size_max")
cluster_size      = param.get_int("cluster_size")
cluster_start     = param.get_int("cluster_start")
random_size_max   = param.get_bool("random_size_max")
random_size       = param.get_int("random_size")
random_start      = param.get_int("random_start")
query_galaxy_only = param.get_bool("query_galaxy_only")
r200_factor       = param.get_float("r200_factor") 
cosmology_name    = param.get_string("cosmology_name")

if "richness_mass_author" in param:
    richness_mass_author = param.get_string("richness_mass_author")
else:
    richness_mass_author = "Rykoff_crit"
if "mass_type" in param:
    mass_type = param.get_string("mass_type")
    assert (mass_type in ['critical', 'crit', 'mean']), "Mass type ("+mass_type+") not understood"
else:
    mass_type = "critical"

if "cluster_use_random_positions" in param:
    cluster_use_random_positions = param.get_bool("cluster_use_random_positions")
else:
    cluster_use_random_positions = False

if "query_random" in param:
    query_random = param.get_bool("query_random")
else:
    query_random = True

if "query_cluster" in param:
    query_cluster = param.get_bool("query_cluster")
else:
    query_cluster = True

cosmo = background.set_cosmology(cosmology_name)
dtk.ensure_dir(query_data_folder)

hdulist = pyfits.open("redmapper_dr8_public_v6.3_catalog.fits")
hdurnd  = pyfits.open("redmapper_dr8_public_v6.3_randoms.fits")



#header list can be found in http://arxiv.org/pdf/1303.3562v2.pdf

tdata =  hdulist[1].data
rdata =  hdurnd[1].data

red_ra = tdata.field('ra')
red_dec= tdata.field('dec')
red_z  = tdata.field('z_lambda')
red_lambda=tdata.field('lambda')
red_pcen=tdata.field('p_cen')

rnd_ra = rdata.field('ra')
rnd_dec= rdata.field('dec')
rnd_z =  rdata.field('z')
rnd_lambda=rdata.field('lambda')


plt.figure()
plt.plot(rnd_ra[:10000])
plt.title("Rnd_ra")

plt.figure()
plt.plot(rnd_dec[:10000])
plt.title("rnd_dec")

for limit in (100, 1000, 10000, 100000,1000000):
    unique_ra = np.unique(rnd_ra[:limit])
    count = unique_ra.size
    print limit, ": ",count
    print "\t", float(count)/float(limit)
    plt.figure()
    plt.plot(rnd_ra[:limit], rnd_dec[:limit], 'x', alpha=0.3)
plt.show()
exit()
# Use the background.py specified  mass-richness definition
# def lambda_to_m200(l):
#     #mass richness relation from http://arxiv.org/abs/1603.06953
#     #m200 is relative to the mean density
# #    m200 = 10**14.344*(l/40.0)**1.33
#     # Hu & Kravtsov 2003
#     m200 = 1e14*np.exp(1.48+1.06*np.log(l/60.0))*0.7 #the the 0.7 is for Msun/h70 not Msun/h100
#     return m200

# def crit_density(z): #Msun/h /kpc^3
#     gcm3_to_msunkpc3 = 1.477543e31
#     density = cosmo.critical_density(z).value*gcm3_to_msunkpc3
#     #print "crit desity(%f): %f Msun/kpc^3"%(z,density)
#     #print "crit desity(%f): %e Msun/kpc^3"%(z,density*1000**3/cosmo.h**2)
#     return density
    
# def m200_to_r200(m200,z): #r200 in kpc
#     r200 = (3.0*m200/(4.0*200.0*np.pi*crit_density(z)))**(1.0/3.0)
#     return r200

# def r200_to_m200(r200,z):
#     m200 = 4.0/3.0*np.pi*crit_denisty(z)*200*r200**3

# def r200_to_arcmin(r200,z):
#     arcmin = r200/cosmo.kpc_proper_per_arcmin(z).value
#     return arcmin

# def lambda_to_arcmin(l,z, richness_mass_author = "Rykoff"):
#     return r200_to_arcmin(m200_to_r200(background.lambda_to_m200c(l, z, richness_mass_author=richness_mass_author),z),z)

#print "\n\ntest:-----------: ", m200_to_r200(1e14,0)

print rnd_lambda.shape
print np.unique(rnd_lambda).shape

if query_cluster:
    print "Querying redmapper clusters..."
    if(cluster_size_max):
        cluster_size = red_ra.size
    if(cluster_use_random_positions ):
        red_ra[:cluster_size]  = rnd_ra[:cluster_size]
        red_dec[:cluster_size] = rnd_dec[:cluster_size]
    query(query_data_folder,red_ra,red_dec,red_z,red_lambda,"gal",cluster_size,start=cluster_start,plot=False)
else:
    print "Not querying redmapper clusters..."

if query_random:
    print "Querying random fields..."
    if(random_size_max):
        random_size = rnd_ra.size
    query(query_data_folder,rnd_ra,rnd_dec,rnd_z,rnd_lambda,"rnd",random_size,start=random_start,plot=False)
else:
    print "Not quering random fields..."


