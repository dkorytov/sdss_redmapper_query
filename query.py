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


def query_sdss_culster(file_loc, cat_ra, cat_dec, cat_z, cat_lambda,
                       name, num, start=0, plot=False, spider_rad = None,
                       query_galaxy_only=True,
                       r200_factor=None,
                       richness_mass_author=None):
    fails = []
    if(query_galaxy_only):
        query_table = "Galaxy"
    else:
        query_table = "PhotoObj"
    print "querying..."
    hfile = h5py.File(file_loc+"query_results.hdf5",mode='a')
    for i in range(start,num):
        #try:
        start = time.time()
        print "%d/%d"%(i+1,num)
        keys = hfile.keys()
        if("%s%d"%(name,i) in keys and "%s_prop%d"%(name,i) in keys):
            continue;

        if "%s%d"%(name,i) in keys:
            del hfile["%s%d"%(name,i)]

        if "%s_prop%d"%(name,i) in keys:
            del hfile["%s_prop%d"%(name,i)]
 
        #query columns are defined here:
        #http://skyserver.sdss.org/dr8/en/help/browser/browser.asp?n=PhotoObjAll&t=U
           
        ra = cat_ra[i]
        dec = cat_dec[i]
        z   = cat_z[i]
        richness = cat_lambda[i]

        # Xray Spiders have their own r200c, so we don't need to compute it. 
        if spider_rad is None: 
            mass, r200 = background.lambda_to_m200_r200(richness,z, richness_mass_author=richness_mass_author)
            rad = background.r200_to_arcmin(r200, z)
        else:
            r200c_deg = spider_rad[i]
            rad = r200c_deg * 60
            r200 = background.arcmin_to_r200(rad, z)
            mass = background.r200c_to_m200c(r200, z)

        hgroup = hfile.create_group("%s_prop%d"%(name,i))
        hgroup.create_dataset("ra",data=ra)
        hgroup.create_dataset("dec",data=dec)
        hgroup.create_dataset("z",data=z)
        hgroup.create_dataset("mass",data=mass)
        hgroup.create_dataset("rad",data=rad)
        hgroup.create_dataset("r200",data=r200)
        hgroup.create_dataset("richness", data = richness)
        ## Query and save objects around target

        query_str = "select  p.ra, p.dec, p.type,p.insidemask,p.flags_g,p.flags_i,p.flags_r,p.cModelMag_u, p.cModelMagErr_u,p.cModelMag_g, p.cModelMagErr_g,p.cModelMag_r, p.cModelMagErr_r,p.cModelMag_i, p.cModelMagErr_i,p.cModelMag_z, p.cModelMagErr_z from "+query_table+" p join dbo.fGetNearbyObjEq(%f,%f,%f) r on p.ObjID = r.ObjID"%(ra, dec, r200_factor*rad)

        result = sqlcl.query(query_str).read()
        # datagal = np.genfromtxt(StringIO(result),names=True,delimiter=', ',dtype=['f8','f8','i2','i1','i8','i8','i8', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4'])
        try:
            datagal = np.genfromtxt(StringIO(result),names=True,skip_header=1,delimiter=',',dtype=['f8','f8','i2','i1','i8','i8','i8', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4'])
        except ValueError as e:
            print(query_str)
            continue

        hgroup = hfile.create_group("%s%d"%(name,i))
        hgroup.create_dataset("ra",data=datagal['ra'])
        hgroup.create_dataset("dec",data=datagal['dec'])
        hgroup.create_dataset("type",data=datagal['type'])
        hgroup.create_dataset("insidemask",data=datagal['insidemask'])
        hgroup.create_dataset("flags_g",data=datagal['flags_g'])
        hgroup.create_dataset("flags_i",data=datagal['flags_i'])
        hgroup.create_dataset("flags_r",data=datagal['flags_r'])
        hgroup.create_dataset("mag_u",data=datagal['cModelMag_u'])
        hgroup.create_dataset("mag_err_u",data=datagal['cModelMagErr_u'])
        hgroup.create_dataset("mag_g",data=datagal['cModelMag_g'])
        hgroup.create_dataset("mag_err_g",data=datagal['cModelMagErr_g'])
        hgroup.create_dataset("mag_r",data=datagal['cModelMag_r'])
        hgroup.create_dataset("mag_err_r",data=datagal['cModelMagErr_r'])
        hgroup.create_dataset("mag_i",data=datagal['cModelMag_i'])
        hgroup.create_dataset("mag_err_i",data=datagal['cModelMagErr_i'])
        hgroup.create_dataset("mag_z",data=datagal['cModelMag_z'])
        hgroup.create_dataset("mag_err_z",data=datagal['cModelMagErr_z'])

        end = time.time()
        print " time: %.2f"%float(end-start)
        if(plot):
            plt.figure()
            legends = ["uknown","cosmic_ray","defect","galaxy","ghost","knownobj","star","trail","sky","notatype"]
            slct1 = datagal['insidemask']==0
            for i in range(0,10):
                slct = (datagal["type"] == i) & slct1
                plt.plot(datagal['ra'][slct],datagal['dec'][slct],'x',label=legends[i])
            plt.legend(loc='best')
            plt.xlabel('ra')
            plt.ylabel('dec')
            plt.show()
        # except ValueError as ie:
        #     print ie
        #     fails.append(i)
        #     print "Failure"
        np.save(file_loc+"fail_indexs.npy",fails)
    hfile.close();

def load_redmapper_cluster_fits(fname):
    #header list can be found in http://arxiv.org/pdf/1303.3562v2.pdf
    cat = {}
    hdulist = pyfits.open(fname)
    tdata =  hdulist[1].data
    cat['ra'] = tdata.field('ra')
    cat['dec']= tdata.field('dec')
    cat['z']  = tdata.field('z_lambda')
    cat['lambda']=tdata.field('lambda')
    return cat

def load_redmapper_randoms_fits(fname):
    #header list can be found in http://arxiv.org/pdf/1303.3562v2.pdf
    cat = {}
    hdulist = pyfits.open(fname)
    tdata =  hdulist[1].data
    cat['ra'] = tdata.field('ra')
    cat['dec']= tdata.field('dec')
    cat['z']  = tdata.field('z')
    cat['lambda']=tdata.field('lambda')
    return cat
    
def load_spider_fits(fname):
    cat = {}
    hdu = pyfits.open(fname)
    data = hdu[1].data
    header = hdu[1].header
    cat['ra'] = data.field('ra')
    cat['dec'] = data.field('dec')
    cat['z'] = data.field('z_lambda')
    cat['lambda'] = np.zeros_like(cat['z'])
    cat['r200c_deg']  = data.field('R200C_DEG')
    cat = spider_nan_clean_up(cat)
    return cat

def spider_nan_clean_up(cat):
    slct = np.ones_like(cat['ra'], dtype=bool)
    for key in cat.keys():
        slct = slct & np.isfinite(cat[key])
    result = {}
    for key in cat.keys():
        result[key] = cat[key][slct]
    return result

def query(param_fname):
    param = dtk.Param(param_fname)
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

    if "spider_clusters" in param:
        spider_clusters = param.get_bool("spider_clusters")
    else:
        spider_clusters = False
    cosmo = background.set_cosmology(cosmology_name)
    dtk.ensure_dir(query_data_folder)


    if query_random:
        print "Querying random fields..."
        rand_cat = load_redmapper_randoms_fits("redmapper_dr8_public_v6.3_randoms.fits")
        if(random_size_max):
            random_size = cat['ra'].size
        query_sdss_culster(query_data_folder, rand_cat['ra'],
                           rand_cat['dec'], rand_cat['z'],
                           rand_cat['lambda'], "rnd", random_size,
                           start=random_start, plot=False,
                           query_galaxy_only=query_galaxy_only,
                           r200_factor=r200_factor,
                           richness_mass_author=richness_mass_author)
    else:
        print "Not quering random fields..."
    
    if query_cluster:
        print "Querying redmapper clusters..."
        print(spider_clusters)
        if not spider_clusters:
            cluster_cat = load_redmapper_cluster_fits("redmapper_dr8_public_v6.3_catalog.fits")
            spider_rad = None
        else:
            cluster_cat = load_spider_fits("/media/luna1/dkorytov/data/spider_xray/catCluster-SPIDERS_RASS_CLUS-v2.0.fits")
            spider_rad = cluster_cat['r200c_deg']
        if(cluster_size_max):
            cluster_size = cluster_cat['ra'].size
        if(cluster_use_random_positions ):
            raise
            cluster_cat['ra'][:cluster_size]  = rand_cat['ra'][:cluster_size]
            cluster_cat['dec'][:cluster_size] = rand_cat['dec'][:cluster_size]
        # query_cat(cat, 
        query_sdss_culster(query_data_folder, cluster_cat['ra'],
                           cluster_cat['dec'], cluster_cat['z'],
                           cluster_cat['lambda'], "gal", cluster_size,
                           start=cluster_start, plot=False,
                           spider_rad=spider_rad,
                           r200_factor=r200_factor,
                           richness_mass_author=richness_mass_author,
        )
    else:
        print "Not querying redmapper clusters..."


if __name__ == "__main__":
    query(sys.argv[1])
