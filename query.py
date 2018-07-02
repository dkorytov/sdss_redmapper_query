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


def lambda_to_m200(l):
    #mass richness relation from http://arxiv.org/abs/1603.06953
    #m200 is relative to the mean density
#    m200 = 10**14.344*(l/40.0)**1.33
    # Hu & Kravtsov 2003
    m200 = 1e14*np.exp(1.48+1.06*np.log(l/60.0))*0.7 #the the 0.7 is for Msun/h70 not Msun/h100
    return m200
def crit_density(z): #Msun/h /kpc^3
    gcm3_to_msunkpc3 = 1.477543e31
    density = cosmo.critical_density(z).value*gcm3_to_msunkpc3
    #print "crit desity(%f): %f Msun/kpc^3"%(z,density)
    #print "crit desity(%f): %e Msun/kpc^3"%(z,density*1000**3/cosmo.h**2)
    return density
    
def m200_to_r200(m200,z): #r200 in kpc
    r200 = (3.0*m200/(4.0*200.0*np.pi*crit_density(z)))**(1.0/3.0)
    return r200
def r200_to_m200(r200,z):
    m200 = 4.0/3.0*np.pi*crit_denisty(z)*200*r200**3
def r200_to_arcmin(r200,z):
    arcmin = r200/cosmo.kpc_proper_per_arcmin(z).value
    return arcmin
def lambda_to_arcmin(l,z):
    return r200_to_arcmin(m200_to_r200(lambda_to_m200(l),z),z)

#print "\n\ntest:-----------: ", m200_to_r200(1e14,0)

print rnd_lambda.shape
print np.unique(rnd_lambda).shape

def query(file_loc,cat_ra,cat_dec,cat_z,cat_lambda,name,num,start=0,plot=False):
    fails = []
    if(query_galaxy_only):
        query_table = "Galaxy"
    else:
        query_table = "PhotoObj"
    print "querying..."
    hfile = h5py.File(file_loc+"query_results.hdf5",mode='a')
    hgroup = hfile.require_group('test')
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
        rad = lambda_to_arcmin(richness,z)
        mass = lambda_to_m200(richness)
        r200 = m200_to_r200(lambda_to_m200(richness),z)#                                                         1                               2                               3                               4                               5
        ## save target properties 
        #cat = np.array([(ra,cat_dec[i],cat_z[i],mass,rad,r200)],
        #dtype=[('ra','f8'),('dec','f8'),('z','f8'),('mass','f4'),('rad_arcmin','f4'),('r200','f4')])
        #np.save(file_loc+"%s_prop%d.npy"%(name,i),cat)

        hgroup = hfile.create_group("%s_prop%d"%(name,i))
        hgroup.create_dataset("ra",data=ra)
        hgroup.create_dataset("dec",data=dec)
        hgroup.create_dataset("z",data=z)
        hgroup.create_dataset("mass",data=mass)
        hgroup.create_dataset("rad",data=rad)
        hgroup.create_dataset("r200",data=r200)
        
        ## Query and save objects around target
        result = sqlcl.query("select  p.ra, p.dec, p.type,p.insidemask,p.flags_g,p.flags_i,p.flags_r,p.cModelMag_u, p.cModelMagErr_u,p.cModelMag_g, p.cModelMagErr_g,p.cModelMag_r, p.cModelMagErr_r,p.cModelMag_i, p.cModelMagErr_i,p.cModelMag_z, p.cModelMagErr_z from "+query_table+" p join dbo.fGetNearbyObjEq(%f,%f,%f) r on p.ObjID = r.ObjID"%(ra,dec,r200_factor*rad)).read()
        result = result.replace(',',', ')
        output = open("test.txt","w")
        output.write(result)
        output.close()
#        datagal = np.genfromtxt(StringIO(result),names=True,delimiter=', ',dtype=['f8','f8','i2','i1','i8','i8','i8', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4'])
        datagal = np.genfromtxt(StringIO(result),names=True,skip_header=1,delimiter=', ',dtype=['f8','f8','i2','i1','i8','i8','i8', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4', 'f4','f4'])
        #np.save(file_loc+"%s%d.npy"%(name,i),datagal)

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


# print "Querying redmapper clusters..."
# if(cluster_size_max):
#     cluster_size = red_ra.size
# query(query_data_folder,red_ra,red_dec,red_z,red_lambda,"gal",cluster_size,start=cluster_start,plot=False)


print "Querying random fields..."
if(random_size_max):
    random_size = rnd_ra.size
query(query_data_folder,rnd_ra,rnd_dec,rnd_z,rnd_lambda,"rnd",random_size,start=30000,plot=False)
#random_start

