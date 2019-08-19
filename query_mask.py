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
import background

from query import load_redmapper_cluster_fits, load_redmapper_randoms_fits, load_spider_fits, spider_mean_from_crit

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
        if(background.rad_dist2(ra,dec,x[i],y[i])<r200+extra_space):
            #mask corner is inside r200(+buffer space)
            #print i,"type:",mask_type,"cen:",ra,dec,"pos:",x[i],y[i],"dist:",rad_dist2(ra,dec,x[i],y[i]),"rad: ",r200
            result = False
            break

    return result
    
def query_sdss_mask(file_loc, cat_ra, cat_dec, cat_z, cat_lambda,
                    name, num, r200_factor=1.0, start=0, plot=False,
                    save_data=True, spider_rad=None,
                    spider_mean=False, richness_mass_author='Rykoff'):
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
        mass, r200 = background.lambda_to_m200_r200(richness, z, richness_mass_author=richness_mass_author)
        if spider_rad is None:
            r200_deg = background.r200_to_arcmin(r200,z)/60.0
            rad = background.r200_to_arcmin(r200, z)*r200_factor
        else:
            r200c_deg = spider_rad[i]
            rad = r200c_deg * 60
            r200 = background.arcmin_to_r200(rad, z)
            mass = background.r200c_to_m200c(r200, z)
            if spider_mean:
                m200m, r200m, c200m =mass_adv.changeMassDefinitionCModel(mass, z, 
                                                                         "200c", "200m",
                                                                         c_model='child18')
                mass = m200m
                r200m_r200c_ratio = r200m/r200
                rad *= r200m_r200c_ratio
                r200 *= r200m_r200c_ratio
   
        print "ra",ra,"dec",dec
        print "z", z, "l", richness, "mass", mass, "r200", r200,
        ## Query and save objects around target
        rad_deg = rad/60.0
        rad_mask = rad_deg*3
        a = np.cos(np.pi*dec/180.0)
        query_str = "SELECT ra, dec, radius, type,area FROM Mask where ra < %f and ra > %f and dec < %f and dec > %f and type != 4"%(ra+rad_mask/a,ra-rad_mask/a,dec+rad_mask,dec-rad_mask)
        result = sqlcl.query(query_str)
        result.readline()
        try:
            data_mask = np.genfromtxt(StringIO(result.read()),names=True,delimiter=",",dtype=['f4','f4','f4','i4','S500'])
            mask_pass = mask_outside_r200(ra,dec,rad_deg,data_mask['area'],data_mask['type'])
        except ValueError as e:
            print("failed query, auto fail")
            mask_pass = False
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

def query_mask(param_fname):
    global num_pass
    global num_fail

    param = dtk.Param(param_fname)
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

    if "richness_mass_author" in param:
        richness_mass_author = param.get_string("richness_mass_author")
    else:
        richness_mass_author = "Rykoff"
    print "Richness mass author: ", richness_mass_author

    if "clusters_query" in param:
        clusters_query = param.get_bool("clusters_query")
    else:
        clusters_query = True

    if "randoms_query" in param:
        randoms_query = param.get_bool("randoms_query")
    else:
        randoms_query = True
    if "spider_clusters" in param:
        spider_clusters = param.get_bool("spider_clusters")
    else:
        spider_clusters = False
    if "spider_mean" in param:
        spider_mean = param.get_bool("spider_mean")
    else:
        spider_mean = False

    dtk.ensure_dir(query_data_folder)
    background.set_cosmology(cosmology_name)

    num_pass = 0
    num_fail = 0
    
    
    if clusters_query:
        print "Querying redmapper clusters..."
        if not spider_clusters:
            cluster_cat = load_redmapper_cluster_fits("redmapper_dr8_public_v6.3_catalog.fits")
            spider_rad = None
        else:
            cluster_cat = load_spider_fits("/media/luna1/dkorytov/data/spider_xray/catCluster-SPIDERS_RASS_CLUS-v2.0.fits")
            spider_rad = cluster_cat['r200c_deg']

        if(cluster_size_max):
            cluster_size = cluster_cat['ra'].size
        query_sdss_mask(query_data_folder, cluster_cat['ra'],
                        cluster_cat['dec'], cluster_cat['z'],
                        cluster_cat['lambda'], "gal", cluster_size,
                        r200_factor=r200_factor, start=cluster_start,
                        plot=False, save_data=True,
                        spider_rad=spider_rad,
                        richness_mass_author=richness_mass_author)

    if randoms_query:
        print "Querying random fields..."
        rand_cat = load_redmapper_randoms_fits("redmapper_dr8_public_v6.3_randoms.fits")
        if(random_size_max):
             random_size = rand_cat['ra'].size
        query_sdss_mask(query_data_folder, rand_cat['ra'],
                        rand_cat['dec'], rand_cat['z'], rand_cat['lambda'],
                        "rnd", random_size, r200_factor=r200_factor,
                        start=random_start, plot=False, save_data=True,
                        richness_mass_author=richness_mass_author)

   
if __name__ == "__main__":
    query_mask(sys.argv[1])
