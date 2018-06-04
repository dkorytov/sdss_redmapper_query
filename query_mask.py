#!/usr/bin/env python2.7
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import sqlcl
import dtk
from astropy.io import fits as pyfits
from StringIO import StringIO
#from astropy.cosmology import WMAP9 as cosmo
import time
import h5py 
from background import *

param = dtk.Param(sys.argv[1])
query_data_folder = param.get_string("query_data_folder")
cosmology_name    = param.get_string("cosmology_name")
cluster_size_max  = param.get_bool("cluster_size_max")
cluster_size      = param.get_int("cluster_size")
cluster_start     = param.get_int("cluster_start")
random_size_max   = param.get_bool("random_size_max")
random_size       = param.get_int("random_size")
random_start      = param.get_int("random_start")

query_galaxy_only = param.get_bool("query_galaxy_only")
r200_factor       = param.get_float("r200_factor") 

dtk.ensure_dir(query_data_folder)

hdulist = pyfits.open("redmapper_dr8_public_v6.3_catalog.fits")
hdurnd  = pyfits.open("redmapper_dr8_public_v6.3_randoms.fits")
set_cosmology(cosmology_name)


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


print "\n\ntest:-----------: ", m200_to_r200(1e14,0)

num_pass = 0
num_fail = 0

def area_str_to_num(s):

    data = s.split(" ")
    return float(data[2]),float(data[3]),float(data[4]),float(data[5]),float(data[6]),float(data[7]),float(data[8]),float(data[9])

def area_str_to_lines(s):
    x = []
    y = []
    data = s.tostring().split()
    for i in range(1,len(data)/2):
        x.append(np.float(data[2*i]))
        y.append(np.float(data[2*i+1].replace("\x00",""))) #somtimes it has trailing nulls for some reason? Minimum buffer size?
    x.append(np.float(data[2]))
    y.append(np.float(data[3]))
    return x,y
def area_to_lenghts(ra1,dec1,ra2,dec2):
    ra_diff = np.abs(ra1-ra2)
    ra_avg = (ra1+ra2)/2.0
    dec_diff = np.abs(dec1-dec2)*np.cos(ra_avg)
    return ra_diff,dec_diff

def mask_outside_r200(ra,dec,r200,mask_area,mask_type,extra_space=0):
    result = True
    if(mask_area.size == 1):
        result = single_mask_outside_r200(ra,dec,r200,mask_area,mask_type,extra_space=extra_space)
        return result
    for i in range(0,len(mask_area)):
        if(not single_mask_outside_r200(ra,dec,r200,mask_area[i],mask_type[i],extra_space)):
           result = False
           break
    return result

def single_mask_outside_r200(ra,dec,r200,mask_area,mask_type,extra_space):
    if(mask_type >=3): # bleeding,bright,
        #we don't cut on this type of mask
        return True
    x,y =  area_str_to_lines(mask_area)
    result = True
    for i in range(0,len(x)):
        if(rad_dist2(ra,dec,x[i],y[i])<r200+extra_space):
            #mask corner is inside r200(+buffer space)
            #print i,"type:",mask_type,"cen:",ra,dec,"pos:",x[i],y[i],"dist:",rad_dist2(ra,dec,x[i],y[i]),"rad: ",r200
            result = False
            break

    return result
    
def query(file_loc,cat_ra,cat_dec,cat_z,cat_lambda,name,num,start=0,plot=False,save_data=True):
    global num_pass
    global num_fail
    fails = []
    print file_loc+name+"_mask.hdf5"
    if(save_data):
        hfile = h5py.File(file_loc+name+"_mask.hdf5",'w')

    for i in range(start,num):
        start = time.time()
        #query columns are defined here:
        #http://skyserver.sdss.org/dr8/en/help/browser/browser.asp?n=PhotoObjAll&t=U
        print "%d/%d "%(i,num),
        
        ra = cat_ra[i]
        dec = cat_dec[i]
        z   = cat_z[i]
        richness = cat_lambda[i]
        rad = lambda_to_arcmin(richness,z)
        mass = lambda_to_m200(richness)
        r200 = m200_to_r200(lambda_to_m200(richness),z)
        r200_deg = r200_to_arcmin(r200,z)/60.0
        print "ra",ra,"dec",dec
        print "z",z,"l",richness,"mass",mass,"r200",r200,"r200 deg",r200_deg
        ## Query and save objects around target
        rad_deg = rad/60.0
        rad_mask = rad_deg*3
        a = np.cos(np.pi*dec/180.0)
        result = sqlcl.query("SELECT ra, dec, radius, type,area FROM Mask where ra < %f and ra > %f and dec < %f and dec > %f and type != 4"%(ra+rad_mask/a,ra-rad_mask/a,dec+rad_mask,dec-rad_mask))
        result.readline()
        data_mask = np.genfromtxt(StringIO(result.read()),names=True,delimiter=",",dtype=['f4','f4','f4','i4','S500'])
        mask_pass = mask_outside_r200(ra,dec,rad_deg,data_mask['area'],data_mask['type'])
        if(mask_pass):
            num_pass +=1
        else:
            num_fail +=1
        print "\tpass:",mask_pass,"\tfract total: %.3f "%(float(num_pass)/float(num_pass+num_fail)),
        if(save_data):
            hgroup = hfile.create_group(str(i))
            dataset = hgroup.require_dataset("mask_pass",(1,1),'u1',mask_pass)
            if(mask_pass):
                dataset[0,0]=True
            else:
                dataset[0,0]=False
            dt_str = h5py.special_dtype(vlen=str)
            dataset = hgroup.require_dataset("mask_points",(data_mask["area"].size,),dtype=dt_str)
            if(data_mask["area"].size==1):
                dataset[0]=data_mask['area']
            else:
                for j in range(0,data_mask["area"].size):
                    dataset[j]=data_mask["area"][j]
        
        if(plot):
            plt.figure()
            plt.title("pass: "+str(mask_pass))
            #plot the current cluster radial bins
            a = np.cos(np.pi*dec/180.0)
            circle = pat.Ellipse((ra,dec),rad_deg*2/a,rad_deg*2,ec='k',fc='none')
            plt.gcf().gca().add_artist(circle)
            #ras,decs = circle_radec(ra,dec,rad_deg/60.0)
            # print ra,dec,rad_deg
            # for j in range(0,ras.size):
            #     dist = rad_dist2(ra,dec,ras[j],decs[j])
            #     print "cirlce test: dist: ",dist,"ra,dec: ",ras[j],decs[j]
            # plt.plot(ras,decs,'-k')
            for j in range(0,data_mask['area'].size):
                if(data_mask['area'].size == 1):
                    x,y = area_str_to_lines(str(data_mask['area']))
                    mask_type = data_mask['type']
                else:
                    x,y = area_str_to_lines(data_mask['area'][j])
                    mask_type = data_mask['type'][j]
                    if(mask_type == 0):
                        plot_type = 'r'
                    elif(mask_type == 1):
                        plot_type = 'b'
                    elif(mask_type == 2):
                        plot_type = 'g--'
                    elif(mask_type == 3):
                        plot_type = 'k:'
                    plt.plot(x,y,plot_type)
            plt.plot([],[],'r',label='bleeding')
            plt.plot([],[],'b',label='bright star')
            plt.plot([],[],'g--',label='trail')
            plt.plot([],[],'k:',label='quality hole')

            plt.legend(loc='best')
            plt.xlabel('ra')
            plt.ylabel('dec')
            plt.xlim([ra-rad_mask,ra+rad_mask])
            plt.ylim([dec-rad_mask,dec+rad_mask])
            plt.show()
        #close the file we have been writing to. 
        end = time.time()
        print " time: %.2f"%float(end-start)
    if(save_data):
        hfile.close()

print "Querying redmapper clusters..."
if(cluster_size_max):
    cluster_size = red_ra.size
query(query_data_folder,red_ra,red_dec,red_z,red_lambda,"gal",cluster_size,start=cluster_start,plot=False,save_data=True)

print "Querying random fields..."
if(random_size_max):
    random_size = rnd_ra.size
query(query_data_folder,rnd_ra,rnd_dec,rnd_z,rnd_lambda,"rnd",random_size,start=random_start,plot=False)

