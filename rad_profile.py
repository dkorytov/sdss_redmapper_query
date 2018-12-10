#!/usr/bin/env python2.7
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.cosmology import WMAP9 as cosmo
from background import *
import dtk
import sys

param = dtk.Param(sys.argv[1])
######################
## Params from file ##
######################
query_results     = param.get_string("query_data_folder")
query_cluster_all = param.get_bool('query_cluster_all')
query_cluster_num = param.get_int('query_cluster_num')
query_type        = param.get_string('query_type')


cosmology_name    = param.get_string("cosmology_name")
galaxy_type       = param.get_int("galaxy_type")
galaxy_weight     = param.get_int("galaxy_weight")
galaxy_color_type = param.get_int("galaxy_color_type")
galaxy_mag_type   = param.get_int("galaxy_mag_type")
galaxy_color_lim  = param.get_float_list("galaxy_color_lim%d"%galaxy_color_type)
galaxy_mag_lim    = param.get_float_list("galaxy_mag_lim%d"%galaxy_mag_type)

background_folder = param.get_string("background_folder")
background_force  = param.get_bool("background_force")
background_all    = param.get_bool("background_all")
background_num    = param.get_int("background_num")
# z bins for background estimate
background_z_start = param.get_float('background_z_start')
background_z_end   = param.get_float('background_z_end')
background_z_num   = param.get_float('background_z_num')

# external data such as spt
spt_plot           = param.get_bool("spt_plot")
spt_file           = param.get_string("spt_file")

#params for the final cluster bining
z_bins_num        = param.get_int("z_bins")
z_bins_start      = param.get_float("z_bins_start")
z_bins_end        = param.get_float("z_bins_end")

mass_bins_num     = param.get_int("mass_bins")
mass_bins_start   = param.get_float("mass_bins_start")
mass_bins_end     = param.get_float("mass_bins_end")

radial_r200_rescale = param.get_bool ('radial_r200_rescale')
radial_bin_num      = param.get_int  ('radial_bin_num')
radial_bin_kpc_min  = param.get_float('radial_bin_kpc_min')
radial_bin_kpc_max  = param.get_float('radial_bin_kpc_max')
radial_bin_r200_min = param.get_float('radial_bin_r200_min')
radial_bin_r200_max = param.get_float('radial_bin_r200_max')

RS_line_plot        = param.get_bool('RS_line_plot')

#params for saving the resulting profiles

result_folder       = param.get_string("result_folder")

if('mstar_cut_offset' in param.data.keys()):
    mstar_cut_offset = param.get_float('mstar_cut_offset')
else:
    mstar_cut_offset = 1.0
#######################
## Processing Params ##
#######################

background_z_bins = np.linspace(background_z_start, background_z_end, background_z_num)

if(galaxy_type == 1):
    galaxy_type_name = "all"
elif(galaxy_type == 2):
    galaxy_type_name = "red"
elif(galaxy_type == 3):
    galaxy_type_name = "non-red"
else:
    print "unknwon galaxy type"
    raise 
if(galaxy_weight == 1):
    galaxy_weight_name = "normal"
elif(galaxy_weight == 2):
    galaxy_weight_name = "red squence"
elif(galaxy_weight == 3):
    galaxy_weight_name = 'red squence rykoff'
else:
    print "unknown galaxy weight"
    raise

if(galaxy_color_type ==1):
    galaxy_color_name = 'g-r'
elif(galaxy_color_type==2):
    galaxy_color_name = 'g-r - RS(z)'
if(galaxy_mag_type ==1):
    galaxy_mag_name = 'i'
elif(galaxy_mag_type==2):
    galaxy_mag_name = 'i-m*(z)'

set_gal_clr_type(galaxy_color_type,galaxy_color_lim)
set_gal_mag_type(galaxy_mag_type,  galaxy_mag_lim)
set_gal_type_weight(galaxy_type,galaxy_weight)
set_cosmology(cosmology_name)
set_mstar_cut(mstar_cut_offset)
################################
## Gettting SPT/external data ##
################################

if(spt_plot):
    dataspt=load_spt_clusters(spt_file)

##############################
## Get background estimates ##
##############################

make_background_estimates(query_results,background_folder,galaxy_type,galaxy_weight,background_z_bins,force=background_force,use_all=background_all,use_num=background_num)
[bkgnd_sqkpc,H_bkgnd_sqkpc] = get_background_estimate_sqkpc(background_folder,galaxy_type,galaxy_weight)
[bkgnd_sqdeg,H_bkgnd_sqdeg] = get_background_estimate_sqdeg(background_folder,galaxy_type,galaxy_weight)


def rad_dist(ra1,dec1,ra2,dec2):
    del_dec = dec1-dec2
    avg_dec = (dec1+dec2)/2.0
    del_ra = (ra1-ra2)*np.cos(dec2)
    return np.sqrt(del_dec**2+del_ra**2)

def deg2rad(deg):
    return (deg/360.0)*(2.0*np.pi)

def rad2deg(rad):
    return (rad/(2.0*np.pi))*360.0

def rad_dist2(ra1,dec1,ra2,dec2):
    ra1=deg2rad(ra1)
    ra2=deg2rad(ra2)
    dec1=deg2rad(dec1)
    dec2=deg2rad(dec2)
    del_dec = dec1-dec2
    del_ra  = ra1-ra2
    radial_dist= 2.0*np.arcsin(np.sqrt(np.sin(del_dec/2.0)**2+np.cos(dec1)*np.cos(dec2)*np.sin(del_ra/2.0)**2))
    return rad2deg(radial_dist)

########################
# Getting cluster data #
########################

print "Loading clusters..."
print query_cluster_num
[dataclstr,datagal,data_pass_mask,clstr_num] = get_clstr(query_type,query_results,query_cluster_num,till_end=query_cluster_all)
#[dataclstr,datagal,clstr_num] = get_clstr("gal",query_results,10)
print clstr_num
print "done."




# Checking the radial profile of the galaxies     
if(True):
    bins = np.linspace(0,3000,20)
    [h,xbins]=np.histogram([],bins)
    [h2,xbins]=np.histogram([],bins)
    bin_avg = (xbins[1:]+xbins[0:-1])/2.0
    bin_area = np.pi*(xbins[1:]**2-xbins[0:-1]**2)
    print(xbins)
    for i in range(0,clstr_num):
        slct = select_gal_old(dataclstr[i]['z'],datagal[i]['m_i'],datagal[i]['m_i_err'],datagal[i]['clr_g-r'],datagal[i]['clr_g-r_err'])
        [hh,_]=np.histogram(datagal[i]['rad_kpc'][slct],bins)
        hh=hh/bin_area
        h2=h2+hh #no background sub
        hhh=hh-bkgnd_sqkpc(dataclstr[i]['z'])
        h=h+hhh #with background sub
        if(False):
            est_bkgnd = bkgnd_sqkpc(dataclstr[i]['z'])
            plt.plot()
            plt.plot(bin_avg,hh,'-x',label='raw')
            plt.plot(bin_avg,hhh,'-x',label='bkgnd sub')
            plt.axvline(dataclstr[i]['r200_kpc'],c='k',label='r200')
            plt.plot(bin_avg,np.ones_like(bin_avg)*est_bkgnd,label='background')
            #plt.yscale('log')
            plt.title('background level: %f'%est_bkgnd)
            plt.legend()
            plt.grid()

            rrad = datagal[i]['rad_kpc']/dataclstr[i]['r200_kpc']
            yy = np.random.rand(rrad.size)
            plt.figure()
            plt.plot(rrad, yy, 'x')
            plt.show()

#################################
## Radial profiles by z & mass ##
#################################

if(radial_r200_rescale):
    radbins = np.linspace(radial_bin_r200_min,
                       radial_bin_r200_max,
                       radial_bin_num)
else:
    radbins = np.linspace(radial_bin_kpc_min,
                       radial_bin_kpc_max,
                       radial_bin_num)

# stacked radial density
[stk_rden,radbins]=np.histogram([],radbins)
stk_bksub = np.zeros_like(stk_rden,dtype='f8')
bin_avg = (radbins[1:]+radbins[0:-1])/2.0
bin_area = np.pi*(radbins[1:]**2-radbins[0:-1]**2)

#z_num = 5
#mass_num = 15
z_bins = np.linspace(z_bins_start,z_bins_end,z_bins_num)
z_bins_avg = (z_bins[0:-1]+z_bins[1:])/2.0
mass_bins = np.logspace(mass_bins_start,mass_bins_end,mass_bins_num)
mass_bins_avg = (mass_bins[0:-1]+mass_bins[1:])/2.0

#clstr_z_bin = np.digitize(dataclstr['z'][:,0],z_bins)
#clstr_m200_bin = np.digitize(dataclstr['m200'][:,0],mass_bins)
#assign a z and mass bin to each cluster 
clstr_z_bin = np.ones(clstr_num,dtype='i4')
clstr_m200_bin = np.ones(clstr_num,dtype='i4')
print clstr_num
for i in range(0,clstr_num):
    clstr_z_bin[i] = np.digitize(np.atleast_1d(dataclstr[i]['z']),z_bins)-1
    clstr_m200_bin[i] = np.digitize(np.atleast_1d(dataclstr[i]['m200']),mass_bins)-1
if(spt_plot):
    dataspt['z_bin_i'] = np.digitize(dataspt['z'],z_bins)-1
    dataspt['mass_bin_i']=np.digitize(dataspt['m200'],mass_bins)-1

def valid_zi_mi(z_i,m_i):
    return (z_i >= 0) and (m_i >= 0) and (z_i < z_bins_num) and (m_i < mass_bins_num)

def sqr200_to_sqkpc(sq_r200,r200_kpc):
    sqkpc_per_sqr200 = r200_kpc**2 
    return sq_r200*sqkpc_per_sqr200
    
# color-mag data
[mag_bins,clr_bins] = get_color_mag_bins()
[h,mtmp,clrtmp]=histogram2d_wrapper([],[],bins=(mag_bins,clr_bins))
[pclr_mi_grid,pclr_gr_grid]=np.meshgrid(mtmp,clrtmp)

## TODO: for kpc radial scale, account for r200 ending in different radial bins
h_radial_profile = {} #
h_radial_profile_err = {} #
h_radial_profile_err2 = {} #
h_radial_profile_cnt = {} #
h_radial_profile_var = {} #
h_cnt = {}
h_Ngal     = {}
h_bkgd_Ngal= {}
h_bkgd_Ngal_var = {}
h_bkgd_Ngal_err = {}
h_mass     = {}
h_clr_mg   = {}
h_zm_rad_clr_mg = {}
h_zm_rad_clr_mg_cnt = {}
h_z_rs_dist = {}
h_z_rs_mag = {}
h_z_rs_clr_err = {}
h_clr_mag_clr_err = {}
h_clr_mag_mag_err = {}
for i in range(0,z_bins_num):
    h_z_rs_dist[i] = []
    h_z_rs_mag[i] = []
    h_z_rs_clr_err[i] =[]
    for j in range(0,mass_bins_num):
        h_radial_profile[i,j]=np.zeros_like(stk_rden)
        h_radial_profile_err[i,j]=np.zeros_like(stk_rden)
        h_radial_profile_err2[i,j]=np.zeros_like(stk_rden,dtype='f8')
        h_radial_profile_cnt[i,j]=np.zeros_like(stk_rden)
        h_cnt[i,j]=0.0
        h_Ngal[i,j]=0.0
        h_bkgd_Ngal[i,j]=0.0
        h_bkgd_Ngal_var[i,j]=[]
        h_bkgd_Ngal_err[i,j]=0.0
        h_mass[i,j]=[]
        h_clr_mg[i,j]=np.zeros((100,40),dtype='f8') # 100 axis: m_i 12->22.   40 axis: clr_g-r -1->3
        for k in range(0,len(radbins)-1):
            h_zm_rad_clr_mg[i,j,k]=np.zeros_like(h,dtype='f8')
            h_zm_rad_clr_mg_cnt[i,j,k]=0.0
            h_radial_profile_var[i,j,k]=[]

for i in range(0,len(mag_bins)):
    for j in range(0,len(clr_bins)):
        h_clr_mag_clr_err[i,j]=[]
        h_clr_mag_mag_err[i,j]=[]

for i in range(0,clstr_num):
    # compute the actual area of the cluster in each radial bin. 
    # inner bins are the same, but the bins further out are set to zero
    # and r200 limit bin is some fraction dep. on far r200 extends into it


    #check if the cluster has no masks in it. If it does, do not use it. 
    if(data_pass_mask[i]==0):
        #does have a mask in the cluster region -> ignore
        continue
    if(i%1000==0):
        print i,"/",clstr_num

    # calculate the area of each radial bin
    bin_area_clstr = np.copy(bin_area)
    if(radial_r200_rescale):
        clstr_r200 = 1.0
    else:
        clstr_r200 = dataclstr[i]['r200_kpc']

    #what does this set of lines do?
    for j in range(0,len(radbins)-1): 
        if(radbins[j] < clstr_r200):
            if(radbins[j+1]>clstr_r200):
                bin_area_clstr[j] = np.pi*(clstr_r200**2-radbins[j]**2)
        else:
            bin_area_clstr[j]=0.0
    # cluster area in what ever rescaling selected
    clstr_area = clstr_r200**2
    # select galaxies for what ever qualifications used 
    slct = select_gal_old(dataclstr[i]['z'],datagal[i]['m_i'],datagal[i]['m_i_err'],datagal[i]['clr_g-r'],datagal[i]['clr_g-r_err'])
    # make sure to select only galaxies within r200
    slct_r200 = select_gal_r200(datagal[i]['rad_kpc'],dataclstr[i]['r200_kpc'])
    slct_r200 = slct_r200 & slct
    weights = weigh_galaxies(dataclstr[i],datagal[i],galaxy_weight)
    # find the radial distribution of galaxies
    if(radial_r200_rescale):
        [rad_count,_]=np.histogram(datagal[i]['rad_kpc'][slct]/dataclstr[i]['r200_kpc'],radbins,weights=weights[slct])
        slct_by_bin = select_gal_radial_bins(datagal[i]['rad_kpc'][slct]/dataclstr[i]['r200_kpc'],radbins)
    else:
        [rad_count,_]=np.histogram(datagal[i]['rad_kpc'][slct],radbins,weights=weights[slct])
        slct_by_bin = select_gal_radial_bins(datagal[i]['rad_kpc'][slct],radbins)
    rad_count = rad_count.astype('f8') # convert int histogram in floats
    rad_density = np.zeros_like(rad_count)
    # calculate the radial density
    for j in range(0,len(rad_count)):
        if(bin_area_clstr[j]>0.0):
            rad_density[j] =rad_count[j]/bin_area_clstr[j]
        else:
            if(rad_density[j] > 0.0):
                raise "This shouldn't happen"

    # add the current radial density to the stack
    stk_rden=stk_rden+rad_density
    
    # calcluated the expected background and get the over density
    if(radial_r200_rescale):
        rad_density_bksub = rad_density-bkgnd_sqkpc(dataclstr[i]['z'])*(dataclstr[i]['r200_kpc']**2)
    else:
        rad_density_bksub = rad_density-bkgnd_sqkpc(dataclstr[i]['z'])
    
    #add the radial over density to the stack
    stk_bksub = stk_bksub+rad_density_bksub
    
    # find the redshift and mass bin the cluster belongs to, add the histogram there
    z_i =clstr_z_bin[i] 
    m_i = clstr_m200_bin[i]
    if(valid_zi_mi(z_i,m_i)):
        # add the background subtracted radial density
        h_radial_profile[z_i,m_i] = rad_density_bksub*bin_area_clstr + h_radial_profile[z_i,m_i]
        h_radial_profile_cnt[z_i,m_i] = h_radial_profile_cnt[z_i,m_i] + rad_count
        # keep track of how many galaxies fall into this z/mass bins
        h_cnt[z_i,m_i] = h_cnt[z_i,m_i] + 1.0
        # accumalate the raw Ngal
        h_Ngal[z_i,m_i]=h_Ngal[z_i,m_i]+np.sum(weights[slct_r200])
        # accumalate the Ngal - background
        Ngal_bkgd_sub  = np.sum(weights[slct_r200])-bkgnd_sqkpc(dataclstr[i]['z'])*(np.pi*dataclstr[i]['r200_kpc']**2)
        h_bkgd_Ngal[z_i,m_i]=h_bkgd_Ngal[z_i,m_i]+Ngal_bkgd_sub
        h_bkgd_Ngal_var[z_i,m_i].append(Ngal_bkgd_sub)
        # remember the mass in this z/mass bins
        h_mass[z_i,m_i].append(dataclstr[i]['m200'])
        # find the clr/mag distribution of the galaxies
        [h,_,_]=histogram2d_wrapper(datagal[i]['mag'][slct],datagal[i]['clr'][slct],
                                    weights=weights[slct],bins=(mag_bins,clr_bins))
        # calucate distribution in cnt/kpc^2/mag^2
    
        h=h/dataclstr[i]['sq_kpc']/get_color_mag_bin_area() - H_bkgnd_sqkpc(dataclstr[i]['z'])
        # correct to selected radial scale
        if(radial_r200_rescale):
            r200_area = np.pi #in r200^2 units
            h=h*dataclstr[i]['sq_kpc']/r200_area
        h_clr_mg[z_i,m_i] += h 
        # find how many sigmas away from the distribution
        if(galaxy_color_type==2):
            sigmas = datagal[i]['clr'][slct]/np.sqrt(datagal[i]['clr_err'][slct]**2+0.05**2)
            h_z_rs_dist[z_i]+=sigmas.tolist()
            h_z_rs_mag[z_i]+= datagal[i]['m_i'][slct].tolist()
            h_z_rs_clr_err[z_i]+= datagal[i]['clr_err'][slct].tolist()

        # select galaxies that belong to particular radial bin, and make a clr/mag distribution for them
        for k in range(0,len(radbins)-1):
            h_radial_profile_var[z_i,m_i,k].append(rad_density_bksub[k]*bin_area_clstr[k])
            slct_bin = slct_by_bin[k]
            # find the current bin area in kpc^2
            if(radial_r200_rescale):
                r200_area = np.pi # in r200^2 units
                bin_area_kpc = bin_area_clstr[k]/r200_area*dataclstr[i]['sq_kpc']
            else:
                bin_area_kpc = bin_area_clstr[k] # already in kpc
            # find clr/mag distribution of this radial bin
            [h_zmr,_,_] = histogram2d_wrapper(datagal[i]['mag'][slct][slct_bin], datagal[i]['clr'][slct][slct_bin], 
                                              weights=weights[slct][slct_bin],   bins=(mag_bins,clr_bins))
            a = h_zmr/bin_area_kpc/get_color_mag_bin_area() - H_bkgnd_sqkpc(dataclstr[i]['z'])
            # correct to /kpc^2/mag^2 and subtract expected background
            if(radial_r200_rescale):
                r200_area = np.pi #in r200^2 units
                a=a*dataclstr[i]['sq_kpc']/r200_area
            # append to stack
            h_zm_rad_clr_mg[z_i,m_i,k] += a
            # the number of clusters that hit this radial bin, partial overlaps count as fractions
            h_zm_rad_clr_mg_cnt[z_i,m_i,k] +=bin_area_clstr[k]/bin_area[k] 
        
        # Finding the mean error in color as function of color/mag
        gal_mi_clrgr_bins = galaxy_mi_clrgr_bins(datagal[i]['m_i'][slct],datagal[i]['clr_g-r'][slct])
        for k in range(0,np.sum(slct)): #iterat over all selected galaxies and record their mag/clr errors
            if(valid_mi_clrgr_bin(gal_mi_clrgr_bins[k])):
                #if(datagal[i]['clr_err'][slct][k] > 3.0):
                #    print datagal[i]['clr_err'][slct][k],datagal[i]['m_i'][slct][k],datagal[i]['m_i_err'][slct][k],datagal[i]['m_r'][slct][k],datagal[i]['m_r_err'][slct][k],datagal[i]['m_g'][slct][k],datagal[i]['m_g_err'][slct][k],datagal[i]['ra'][slct][k],datagal[i]['dec'][slct][k]
                h_clr_mag_clr_err[gal_mi_clrgr_bins[k,0],gal_mi_clrgr_bins[k,1]].append(datagal[i]['clr_err'][slct][k])
                h_clr_mag_mag_err[gal_mi_clrgr_bins[k,0],gal_mi_clrgr_bins[k,1]].append(datagal[i]['mag_err'][slct][k])


################################################################
# Adjust the accumlated stats to correct (i.e. divide by cnts) #
################################################################

for i in range(0,z_bins_num-1):
    for j in range(0,mass_bins_num-1):
        if(h_cnt[i,j]!=0):
            h_radial_profile[i,j]=h_radial_profile[i,j]/(h_cnt[i,j]*bin_area)
            h_bkgd_Ngal[i,j]= h_bkgd_Ngal[i,j]/np.float(h_cnt[i,j])
            h_bkgd_Ngal_var[i,j]=np.std(h_bkgd_Ngal_var[i,j])
            h_bkgd_Ngal_err[i,j]=h_bkgd_Ngal_var[i,j]/np.sqrt(h_cnt[i,j])
            h_clr_mg[i,j]=h_clr_mg[i,j]/np.float(h_cnt[i,j])
            h_radial_profile_err2[i,j]=np.sqrt(h_radial_profile_cnt[i,j])/(np.float(h_cnt[i,j])*bin_area)
            #h_bkgd_Ngal_err[i,j] = h_radial_profile_err2[i,j]
            for k in range(0,radial_bin_num-1):
                h_zm_rad_clr_mg[i,j,k]=h_zm_rad_clr_mg[i,j,k]/h_zm_rad_clr_mg_cnt[i,j,k]
                h_radial_profile_var[i,j,k]= np.std(np.array(h_radial_profile_var[i,j,k])/(h_cnt[i,j]*bin_area[k]))
                h_radial_profile_err[i,j,k]=h_radial_profile_var[i,j,k] #/h_cnt[i,j]
        else:
            h_radial_profile[i,j]=h_radial_profile[i,j]/np.nan
            h_bkgd_Ngal_var[i,j]=0.0
            h_bkgd_Ngal_err[i,j]=0.0
            for k in range(0,radial_bin_num-1):
                h_radial_profile_var[i,j,k]=0.0
                h_radial_profile_err[i,j,k]=0.0



clr_err = {}
mag_err = {}

for i in range(0,len(mag_bins)-1):
    # print "\nL ",               
    for j in range(0,len(clr_bins)-1):
        # print len(h_clr_mag_clr_err[i,j]),
        if(len(h_clr_mag_clr_err[i,j])>0):
            clr_err[i,j]= np.average(np.array(h_clr_mag_clr_err[i,j],dtype='f4'))
            mag_err[i,j]= np.average(np.array(h_clr_mag_mag_err[i,j],dtype='f4'))
        else:
            clr_err[i,j]=-1.0
            mag_err[i,j]=-1.0

h_clr_mag_clr_err = np.ones((len(mag_bins),len(clr_bins)))
h_clr_mag_mag_err = np.zeros_like(h_clr_mag_clr_err)

for i in range(0,len(mag_bins)-1):
    for j in range(0,len(clr_bins)-1):
        h_clr_mag_clr_err[i,j]=clr_err[i,j]
        h_clr_mag_mag_err[i,j]=mag_err[i,j]

mi_bins,clrgr_bins= get_mi_clrgr_bins()
save_error(background_folder,h_clr_mag_clr_err,h_clr_mag_mag_err,mi_bins,clrgr_bins)

plt.figure()
plt.pcolor(mi_bins,clrgr_bins,h_clr_mag_clr_err.T,cmap=plt.cm.BuPu,norm=LogNorm())
plt.colorbar()

plt.title("Color Errors")
#plt.figure()
#plt.pcolor(mag_bins,clr_bins,h_clr_mag_mag_err.T,cmap=plt.cm.BuPu,norm=LogNorm())
#plt.colorbar()
#plt.title("Magnitude Errors")





################################################
# save the clr/mag radial profile into a file  #
################################################
save_radial_profile(result_folder,
                    galaxy_type,
                    galaxy_weight,
                    z_bins,
                    mass_bins,
                    h_cnt, #z mass count histogram
                    radbins,
                    radial_r200_rescale,
                    h_zm_rad_clr_mg,
                    h_zm_rad_clr_mg_cnt,
                    h_bkgd_Ngal,
                    h_bkgd_Ngal_err,
                    h_bkgd_Ngal_var,
                    h_radial_profile,
                    h_radial_profile_err2,
                    h_radial_profile_var)
# #save_radial_profile(result_folder,
#                     galaxy_type,
#                     galaxy_weight,
#                     z_bins,
#                     m_bins,
#                     r_bins,
#                     zm_Ngal,
#                     zm_Ngal_err,
#                     zm_Ngal_var,
#                     zmr_gal_density,
#                     zmr_gal_density_err,
#                     zmr_gal_density_var,
#                     zm_counts)

print clstr_num        
h=h/clstr_num
stk_bksub=stk_bksub/clstr_num
stk_rden/=clstr_num
plt.figure()
plt.plot(bin_avg,stk_rden,label='galaxies raw')
plt.plot(bin_avg,stk_bksub,label='galaxies bkgnd sub')
plt.legend()
plt.yscale('log')
if(radial_r200_rescale):
    plt.ylabel('galaxies/r200^2')
    plt.xlabel('radius [r200]')
else:
    plt.ylabel('galaxies/kpc^2')
    plt.xlabel('radius [kpc]')

plt.grid()


##############################
# make all the summery plots #
##############################
for i in range(0,z_bins_num-1):
    z_h = np.zeros_like(h)
    plt.figure(figsize=(12,12))

    # Radial Profile
    plt.subplot(2,2,1)
    plt.title('%.2f<z<%.2f'%(z_bins[i],z_bins[i+1]))
    mass_limit=[]
    log_scale = False
    for j in range(0,mass_bins_num-1):
        if(h_cnt[i,j]!=0):
            mass_limit.append(mass_bins_avg[j])
    if(len(mass_limit)==0):
        mass_clrs = dtk.get_colors(mass_bins_avg,plt.cm.BuPu,log=False,vmax=1.0,vmin=0.0)
    else:
        mass_clrs = dtk.get_colors(mass_bins_avg,plt.cm.copper_r,log=True,vmax=np.max(mass_limit),vmin=np.min(mass_limit))
    for j in range(0,mass_bins_num-1):
        if(h_cnt[i,j]!=0):
            plt.plot(bin_avg,h_radial_profile[i,j],label='%.2e'%mass_bins_avg[j],c=mass_clrs[j])
            error = h_radial_profile_err2[i,j]
            #error = np.zeros(bin_avg.size)
            #for k in range(0,bin_avg.size):
            #    error[k]=h_radial_profile_err2[i,j,k]
            plt.fill_between(bin_avg,h_radial_profile[i,j]-error+0.01,h_radial_profile[i,j]+error,color=mass_clrs[j],alpha=0.5)
            if(np.sum(h_radial_profile[i,j]>0) > 0):
                log_scale = True
    plt.legend(prop={'size':12},framealpha=0.5)
    if(log_scale):
        plt.yscale('log')
        
    if(radial_r200_rescale):
        plt.ylim(ymin=1e-1)
    else:
        plt.ylim(ymin=1e-6)
    plt.grid()
    if(radial_r200_rescale):
        plt.ylabel('galaxies/r200^2')
        plt.xlabel('radius [r200]')
    else:
        plt.ylabel('galaxies/kpc^2')
        plt.xlabel('radius [kpc]')

    # Ngal plot
    plt.subplot(2,2,2)
    for j in range(0,mass_bins_num-1):
        if(h_cnt[i,j]!=0):
            plt.plot(mass_bins_avg[j],np.nansum(h_radial_profile[i,j]*bin_area),'bo')
            plt.plot(mass_bins_avg[j],h_Ngal[i,j]/h_cnt[i,j],'ro')
            plt.plot(mass_bins_avg[j],h_bkgd_Ngal[i,j],'g^',mfc='none')
    
    plt.plot([],[],'bo',label='Ngal density')
    plt.plot([],[],'ro',label='Ngal raw')
    plt.plot([],[],'g^',label='Ngal -bkgnd')
    if(spt_plot):
        spt_slct = dataspt['z_bin_i']==i
        plt.plot(dataspt['m200'][spt_slct],dataspt['ngal200'][spt_slct],'*c',label='SPT Ngal200')
        plt.plot(dataspt['m200'][spt_slct],dataspt['ngalRS200'][spt_slct],'xm',label='SPT NgalRS200')
    plt.ylabel('Ngal,200')
    plt.legend(loc=2,framealpha=0.5)
    plt.xlabel('cluster M200 [Msun/h]')
    plt.xscale('log')
    xlim = plt.xlim()
    plt.grid()

    # Mass hist
    plt.subplot(2,2,3)
    masses = []
    for j in range(0,mass_bins_num-1):
        masses.extend(h_mass[i,j])
    [hist_mass,_]=np.histogram(np.array(masses),bins=mass_bins)
    plt.title('cluster sum = %d'%np.sum(hist_mass))
    plt.step(mass_bins_avg,hist_mass+1.0, where='mid')
    plt.xscale('log')
    plt.xlabel('cluster M200 [Msun/h]')
    plt.grid()
    plt.yscale('log')
    plt.ylabel('count + 1')
    plt.xlim(xlim)

    # Mag-Clr plot
    plt.subplot(2,2,4)
    plt.title('stacked clr,i-mag galaxy distribution')
    h_all = np.zeros_like(h_clr_mg[i,0])
    for j in range(0,mass_bins_num-1):
        if(hist_mass[j]>0):
            h_all +=h_clr_mg[i,j]*hist_mass[j]/np.sum(hist_mass)
    #print  h_all.max(),h_all.max()*1e-5 ,h_all.min()
    vmin = np.max((h_all.max()*1e-5,h_all.min()))
    if(np.sum(h_all>0)>0): # if there positive values -> log scale them
        plt.pcolor(mag_bins,clr_bins,h_all.T,cmap=plt.cm.BuPu,norm=LogNorm(vmin=vmin))
    else:
        plt.pcolor(mag_bins,clr_bins,h_all.T,cmap=plt.cm.BuPu)

    if(np.max(h_all.clip(min=0.0)) != np.min(h_all.clip(min=0.0))):
        if(radial_r200_rescale):
            a=1
            plt.colorbar().set_label('cnt/r200^2/mag^2')
        else:
            a=1
            plt.colorbar().set_label('cnt/kpc^2/mag^2')
    if(True):#np.sum(h_all)>0):
        contain_lvls = dtk.get_sigmas([0.5,1,2])#[0.68,0.86,0.95,0.97]
        smooth_h = dtk.smoothen_H(h_all)
        c_lvls=dtk.conf_interval(smooth_h,contain_lvls)
        tmp2 = dtk.contour_labels(contain_lvls,c_lvls)
        [m_avgs,clr_avgs]= get_color_mag_bin_avgs()
        cs = plt.contour(m_avgs,clr_avgs,smooth_h.T,c_lvls,colors='k')
        #plt.clabel(cs,fmt = tmp2,colors='k') <- This doesn't work anymore?!
    if(RS_line_plot):
        RS_scatter = 0.05
        if(galaxy_color_type==1):
            plt.plot(mag_bins,redsequence_line(mag_bins,z_bins[i]),'r--',label='Min/Max RS(z) mean')
            plt.plot(mag_bins,redsequence_line(mag_bins,z_bins[i+1]),'r--')
            plt.plot(mag_bins,redsequence_line(mag_bins,z_bins[i])-RS_scatter,'r:',label='RS scatter')
            plt.plot(mag_bins,redsequence_line(mag_bins,z_bins[i+1])+RS_scatter,'r:')
        elif(galaxy_color_type==2):
            plt.axhline(RS_scatter,color='r',ls='--')
            plt.axhline(-RS_scatter,color='r',ls='--')
            #plt.axhline(RS_scatter,color='r',ls='..',label="RS Scatter")
            #plt.axhline(-RS_scatter,color='r',ls='..')
                        #plt.axhline(0.0,'r--',label="RS")
    plt.ylabel(galaxy_color_name)
    plt.xlabel(galaxy_mag_name)
    plt.legend()
    plt.tight_layout()
    if(False):
        plt.figure(figsize=(15,4))
   	plt.subplot(1,4,1)
   	H_radial_sum = np.zeros_like(h_all)
   	total_radial_area = np.pi
   	for j in range(0,mass_bins_num-1):
   	    H_radial_mass_sum = np.zeros_like(H_radial_sum)
   	    for k in range(0,radial_bin_num-1):
   	        H_radial_mass_sum +=h_zm_rad_clr_mg[i,j,k]*bin_area[k]/np.pi
   	    H_radial_sum += H_radial_mass_sum*hist_mass[j]/np.sum(hist_mass)
   	plt.pcolor(mag_bins,clr_bins,H_radial_sum.T,cmap=plt.cm.BuPu,norm=LogNorm())
   	if(np.max(H_radial_sum) != np.min(H_radial_sum)):
   	    plt.colorbar()
   	plt.title("radial sum: %e"%np.sum(H_radial_sum))
   	plt.subplot(1,4,2)
   	plt.pcolor(mag_bins,clr_bins,h_all.T,cmap=plt.cm.BuPu,norm=LogNorm())
   	if(np.max(h_all) != np.min(h_all)):
   	    plt.colorbar()
   	plt.title("normal sum: %e"%np.sum(h_all))
   	plt.subplot(1,4,3)
   	h_diff = h_all - H_radial_sum
   	print np.nanmean(h_all/H_radial_sum), 1.0/ np.nanmean(h_all/H_radial_sum)
   	plt.pcolor(mag_bins,clr_bins,h_diff.T,cmap=plt.cm.RdBu)
   	plt.title("diff sum: %e"%np.sum(h_diff))
   	plt.subplot(1,4,4)
   	plt.pcolor(mag_bins,clr_bins,H_bkgnd_sqkpc([z_bins[i+1]]).T,cmap=plt.cm.BuPu,norm=LogNorm())
   	plt.colorbar()
   	plt.title("bkgnd sum: %e"%np.sum(H_bkgnd_sqkpc([z_bins[i+1]])))
   	plt.tight_layout()






dtk.save_figs("figs/"+param.file+"/rad_profile/")
plt.show()


radial_profile = load_radial_profile(result_folder,galaxy_type,galaxy_weight)
