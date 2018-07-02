#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from astropy.cosmology import WMAP9 as cosmowmap9
from astropy.cosmology import Planck13 as cosmoplanck15
from astropy.cosmology import WMAP7 as cosmowmap7
import scipy.interpolate
import dtk
import os
import h5py

cosmo=None
def set_cosmology(name):
    #sets the cosmology related for the angluar distance -> kpc Phys
    #conversion
    global cosmo
    if(name == "planck13"):
        cosmo = cosmoplanck15
    elif(name == "wmap9"):
        cosmo = cosmowmap9
    elif(name == 'wmap7'):
        cosmo = cosmowmap7
    else:
        print "Cosmology not defined"
        raise()
    return cosmo
def lambda_to_m200(l):
    #mass richness relation from http://arxiv.org/abs/1603.06953 
    #m200 is relative to the mean density not crit
    #    m200 = 10**14.344*(l/40.0)**1.33

    #mass richness relation from http://iopscience.iop.org/article/10.1088/0004-637X/746/2/178/pdf
    # Hu & Kravtsov 2003 
    # m200 relative to crit
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
    return m200

def r200_to_arcmin(r200,z):
    arcmin = r200/cosmo.kpc_proper_per_arcmin(z).value
    return arcmin

def lambda_to_arcmin(l,z):
    return r200_to_arcmin(m200_to_r200(lambda_to_m200(l),z),z)


def get_clstr(name,folder,num,start=0,till_end=False):
    # till_end = true will make the loop continue until 
    # it fails the find a file. 
    print "\tloading query data..."
    print "\t\tstart: ",start
    if(till_end):
        print "\t\ttill_end: ",till_end
    else:
        print "\t\tnum: ", num

    # the final outputs. Both are lists with one dictionary for each
    # cluster. The dataclstr dictionary has properties of the cluster 
    # (r200, mass etc).The datagal dictionary has np.arrays for galaxy
    # properties fo r
    all_datagal   = [] 
    all_dataclstr = []
    all_mask_pass = []
    # if till_end = false, i=j at all times
    i = start
    j = start
    if(till_end):
        j = num-1
    print folder+'query_results.hdf5'
    print folder+name+'_mask.hdf5'
    hquery_results = h5py.File(folder+'query_results.hdf5','r')
    hmask_results  = h5py.File(folder+name+'_mask.hdf5','r')
    print "\t"+folder+'query_results.hdf5'
    while ((j<num+1) | till_end):
        j = i
        i=i+1
        if(i%1000==0):
            print "\t",i
        j=j+1
        # if(till_end):
        #        j=num-2
        # for i in range(start,num):
        try:
            #adding support for hdf5 query outputs
            #datagal = np.load(folder+"%s%d.npy"%(name,i))
            #dataclstr = np.load(folder+"%s_prop%d.npy"%(name,i))
            datagal = {}
            dataclstr = {}
            galgroup = hquery_results[name+'%d'%j] #individual galaxy properties
            datagal['dec']       = galgroup['dec'][:]
            datagal['flags_g']   = galgroup['flags_g'][:]   
            datagal['flags_i']   = galgroup['flags_i'][:]   
            datagal['flags_r']   = galgroup['flags_r'][:]   
            datagal['insidemask']= galgroup['insidemask'][:]
            datagal['m_g_err']   = galgroup['mag_err_g'][:] 
            datagal['m_i_err']   = galgroup['mag_err_i'][:] 
            datagal['m_r_err']   = galgroup['mag_err_r'][:] 
            datagal['m_u_err']   = galgroup['mag_err_u'][:] 
            datagal['m_z_err']   = galgroup['mag_err_z'][:] 
            datagal['m_g']       = galgroup['mag_g'][:]     
            datagal['m_i']       = galgroup['mag_i'][:]     
            datagal['m_r']       = galgroup['mag_r'][:]     
            datagal['m_u']       = galgroup['mag_u'][:]  
            datagal['m_z']       = galgroup['mag_z'][:]     
            datagal['ra']        = galgroup['ra'][:]        
            datagal['type']      = galgroup['type'][:]      

            gpgroup = hquery_results[name+'_prop%d'%j] #cluster properties 
            dataclstr['dec']  = gpgroup['dec'].value
            dataclstr['mass'] = gpgroup['mass'].value
            dataclstr['r200'] = gpgroup['r200'].value
            dataclstr['ra']   = gpgroup['ra'].value
            dataclstr['rad']  = gpgroup['rad'].value
            dataclstr['z']    = gpgroup['z'].value
            
            mask_pass = hmask_results['%d'%j]['mask_pass'][:]
            #what does this do? Skip the thing if it has no galaxies?
            #will comment out
            #if(datagal.size ==0):
            #     continue
        except KeyError as ie:
            print "\tget_clstr: read %d files in %s type: %s"%((i-1),folder,name)
            i-=1
            break
        datagal_i   ={}
        dataclstr_i ={}
        datagal_names = datagal.keys()
        dataclstr_names = dataclstr.keys()

        # flags defined at http://skyserver.sdss.org/dr8/en/help/browser/enum.asp?n=PhotoFlags
        #          bright               saturated             satur_center         nopetro             deblended_as_moving
        cut_flgs = 0x0000000000000002 | 0x0000000000040000  | 0x0000080000000000 | 0x0000000000000100 | 0x0000000100000000
        flgs = datagal['flags_i'] | datagal['flags_r'] | datagal['flags_g']
        slct1 = datagal['type']!=31
        slct2 = np.logical_not(flgs&cut_flgs)
        slct = slct1 & slct2
        # add all galaxy data saved to the output
        for nm in datagal_names:
            datagal_i[nm]=datagal[nm][slct]
        # rename a the r200 fields to be more clear
        for nm in dataclstr_names:
            if(nm=='r200'):
                nm2 ='r200_kpc'
            #elif(nm=='rad_arcmin'): old
            elif(nm=='rad'):
                nm2 ='r200_arcmin'
            elif(nm=='mass'):
                nm2 ='m200'
            else:
                nm2=nm
            dataclstr_i[nm2]=dataclstr[nm]

        # Give mags easier aliases to use
        # this is done above
        #datagal_i['m_u']=datagal_i['mag_u']
        #datagal_i['m_g']=datagal_i['mag_g']
        #datagal_i['m_r']=datagal_i['mag_r']
        #datagal_i['m_i']=datagal_i['mag_i']
        #datagal_i['m_z']=datagal_i['mag_z']
        #datagal_i['m_u_err']=datagal_i['mag_err_u']
        #datagal_i['m_g_err']=datagal_i['mag_err_g']
        #datagal_i['m_r_err']=datagal_i['mag_err_r']
        #datagal_i['m_i_err']=datagal_i['mag_err_i']
        #datagal_i['m_z_err']=datagal_i['mag_err_z']
        datagal_i['clr_g-r']=datagal_i['m_g']-datagal_i['m_r']
        datagal_i['clr_g-r_err']=np.sqrt(datagal_i['m_u_err']**2+datagal_i['m_r_err']**2)

        set_datagal_z_dep(dataclstr['z'],datagal_i)
        # calculate distance 
        ra1 = datagal_i['ra']
        dec1 = datagal_i['dec']
        ra2 = dataclstr['ra']
        dec2 = dataclstr['dec']
        a_rad = rad_dist2(ra1,dec1,ra2,dec2)
        datagal_i['rad_deg']=a_rad
        datagal_i['rad_kpc']=a_rad*cosmo.kpc_proper_per_arcmin(dataclstr['z']).value*60.0

        
        dataclstr_i['sq_deg']= np.pi*dataclstr_i['r200_arcmin']**2/3600.0
        dataclstr_i['sq_kpc']= np.pi*dataclstr_i['r200_kpc']**2
        dataclstr_i['sq_arcmin']=np.pi*dataclstr_i['r200_arcmin']**2
        dataclstr_i['gal_cnt']=-1.0

        all_datagal.append(datagal_i)
        all_dataclstr.append(dataclstr_i)
        all_mask_pass.append(mask_pass)
        clstr_num = i # the number of clusters in the output
    print "\tdone."
    return [all_dataclstr,all_datagal,all_mask_pass,clstr_num]

def set_datagal_z_dep(z,datagal_i):
    # selecting galaxy mag to selected type
    if(_gal_mag_type==1): # i-band mag
        datagal_i['mag']=datagal_i['m_i']
        datagal_i['mag_err']=datagal_i['m_i_err']
    elif(_gal_mag_type==2): # m* normed mag
        datagal_i['mag']= datagal_i['m_i'] - mstar(z) 
        datagal_i['mag_err']= datagal_i['m_i_err']
        
    # Setting galaxy clr to selected type
    if(_gal_clr_type==1): # g-r color
        datagal_i['clr']=datagal_i['clr_g-r']
        datagal_i['clr_err']=datagal_i['clr_g-r_err']
    elif(_gal_clr_type==2): # g-r - Rs(z)   Red squence normalized color
        expected_Rs = redsequence_line(datagal_i['m_i'],z)
        datagal_i['clr']=datagal_i['clr_g-r']-expected_Rs
        datagal_i['clr_err']=datagal_i['clr_g-r_err']


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


def circle_radec(ra1,dec1,radius):
    print "circle_radec"
    ra1 = deg2rad(ra1)
    dec1 = deg2rad(dec1)
    size = 12
    del_decs = np.sin(np.pi*np.linspace(0,1,size))*radius
    del_ras1 = del_ra(ra1,dec1,radius,del_decs)
    del_ras2 = -del_ras1
    decs = np.zeros(size*2)
    ras  = np.zeros(size*2)
    decs[:size]=dec1 + del_decs
    decs[size:]=dec1 + del_decs
    ras [:size]=ra1  + del_ras1
    ras [size:]=ra1  + del_ras2
    return rad2deg(ras),rad2deg(decs)

def del_ra(ra,dec,rad,del_dec):
    return 2*np.arcsin(np.sqrt((np.sin(rad/2.0)**2-np.sin(del_dec/2.0)**2)/(np.cos(dec)*np.cos(dec+del_dec))))


def histogram2d_wrapper(datax,datay,weights=None,bins=None):
    #return np.histogram2d(datax,datay,bins=bins)
    return cell_in_cloud(datax,datay,bins=bins,weights=weights)

def deposite(xi,yi,val,result,num_x,num_y):
    if(xi >= 0 and xi < num_x and yi >= 0 and yi < num_y):
        result[xi,yi]+=val

def cell_in_cloud(datax,datay,bins=None,weights=None):
    # bins must be equal spaced
    min_x = bins[0].min()
    max_x = bins[0].max()
    min_y = bins[1].min()
    max_y = bins[1].max()
    num_x = bins[0].size-1
    num_y = bins[1].size-1
    del_x = (max_x-min_x)/num_x
    del_y = (max_y-min_y)/num_y
    A = del_x*del_y
    result = np.zeros([num_x,num_y],dtype='f8')
    if(weights is None):
        weights = np.ones_like(datax)
    for i in range(0,len(datax)):
        xi = int((datax[i]-min_x)/del_x)
        yi = int((datay[i]-min_y)/del_y)
        a1=  (datax[i]-(xi*del_x+min_x))*(datay[i]-(yi*del_y+min_y))/A
        a2=  (datax[i]-(xi*del_x+min_x))*(-datay[i]+((yi+1)*del_y+min_y))/A
        a3= (-datax[i]+((xi+1)*del_x+min_x))*(datay[i]-(yi*del_y+min_y))/A
        a4= (-datax[i]+((xi+1)*del_x+min_x))*(-datay[i]+((yi+1)*del_y+min_y))/A
        w = weights[i]
        deposite(xi  , yi  , a1*w, result, num_x, num_y)
        deposite(xi  , yi+1, a2*w, result, num_x, num_y)
        deposite(xi+1, yi  , a3*w, result, num_x, num_y)
        deposite(xi+1, yi+1, a4*w, result, num_x, num_y)
    return [result,bins[0],bins[1]]

_mstar_cut = 1
def set_mstar_cut(mstar_cut):
    global _mstar_cut
    _mstar_cut = mstar_cut

def mstar(z):
    # defined in Rykoff 2011 http://arxiv.org/abs/1104.2089 eq. 11 
    # only good for 0.05 < z < 0.35
    return 12.27+62.36*z-289.79*z**2+729.69*z**3-709.42*z**4

def m_cut(z):
    return mstar(z)+_mstar_cut



def select_gal_old(z,gal_i,gal_i_err,gal_clr,gal_clr_err):
    m17 = -0.0701*z-0.008969
    b17 = 3.2982*z+0.5907
    mean_clr = (m17*(gal_i-17.0)+b17)
    clr_diff = mean_clr-gal_clr
    #scatter =np.sqrt(0.15**2+gal_clr_err**2)
    #slct2 = (np.abs(clr_diff) <= 3.0*scatter)
    slct1 = (gal_i<m_cut(z)) & (12.0<gal_i) & (gal_clr > -1.0) & (gal_clr < 3.0)
    #slct = slct1 & slct2
    return slct1 #no red squence selection at the moment, just all galaxie
    

def calc_guassian_overlap(u1,s1,u2,s2):
    # calculates the integral of the product of two guassians
    # equation from http://www.tina-vision.net/docs/memos/2003-003.pdf
    # checked against numerical integration
    return np.exp(-(u1-u2)**2/(2*(s1**2+s2**2)))/np.sqrt(2*np.pi*(s1**2+s2**2))
    

def redsequence_line(m_i,z):
    # defined in Rykoff 2011 http://arxiv.org/abs/1104.2089 eq. 21-23
    m17 = -0.0701*z-0.008969
    b17 = 3.2982*z+0.5907
    mean_clr = (m17*(m_i-17.0)+b17)
    return mean_clr

def redsequence_line_scatter(m_i_err,z):
    m17 = -0.0701*z-0.008969
    mi_scatter = m17*m_i_err
    return np.sqrt(0.05**2+mi_scatter**2)

def select_gal_luminonsity(z,gal_i,gal_i_err,gal_clr,gal_clr_err):
    slct = select_gal(z,gal_i,gal_i_err,gal_clr,gal_clr_err)
    return slct.astype('float')*gal_i

def select_gal_RS_Pm(z,gal_i,gal_i_err,gal_clr,gal_clr_err):
    slct = select_gal(z,gal_i,gal_i_err,gal_clr,gal_clr_err)
    return slct

def select_gal_r200(clstr_gal_radius,clstr_r200):
    return clstr_gal_radius < clstr_r200

def select_gal_radial_bins(gal_radius,radial_bins):
    result = {}
    if(len(gal_radius)==0):
        for i in range(0,len(radial_bins)):
            result[i]=np.zeros(0,dtype=bool)
        return result
    dig = np.digitize(gal_radius,radial_bins)
    for i in range(0,len(radial_bins)):
        result[i]=dig==(i+1)
    return result

def weigh_galaxies(dataclstr_i,datagal_i,gal_weight):
    if(gal_weight == 1):
        return np.ones(len(datagal_i['m_i']))
    if(gal_weight == 2):
        return weigh_red_squence(dataclstr_i,datagal_i)
    if(gal_weight == 3):
        return weigh_red_squence_rykoff(dataclstr_i,datagal_i)

def weigh_red_squence(dataclstr_i,datagal_i):
    z = dataclstr_i['z']
    [m_avg,gr_avg] = get_color_mag_bin_avgs()
    # guassian pdf of galaxy color
    u1s = datagal_i['clr_g-r'] # mean value measured
    s1s = datagal_i['clr_g-r_err'] # 1 sigma scatter
    # guassian pdf of red squence by i-band mag
    u2s = redsequence_line(datagal_i['m_i'],z) 
    s2s = redsequence_line_scatter(datagal_i['m_i_err'],z)
    return calc_guassian_overlap(u1s,s1s,u2s,s2s)

def weigh_red_squence_rykoff(dataclstr_i,datagal_i):
    z = dataclstr_i['z']
    [m_avg,gr_avg] = get_color_mag_bin_avgs()
    # guassian pdf of galaxy color
    u1s = datagal_i['clr_g-r'] # mean value measured
    s1s = datagal_i['clr_g-r_err'] # 1 sigma scatter
    # guassian pdf of red squence by i-band mag
    u2s = redsequence_line(datagal_i['m_i'],z) 
    s2s = redsequence_line_scatter(datagal_i['m_i_err'],z)
    sigma = np.sqrt(s2s**2+s1s**2)
    dist = u1s-u2s
    return np.exp(-dist**2/(2.0*sigma**2))/(np.sqrt(2.0*np.pi)*sigma)

def galaxy_mi_clrgr_bins(gal_mag,gal_clr):
    mi_bins ,  clrgr_bins = get_mi_clrgr_bins()
    clr_mag_bin = np.zeros((len(gal_clr),2))
    if(len(gal_clr)>0):
        clr_mag_bin[:,0]=np.digitize(gal_mag,bins=mi_bins)-1
        clr_mag_bin[:,1]=np.digitize(gal_clr,bins=clrgr_bins)-1
    return clr_mag_bin

def valid_mi_clrgr_bin(clr_mag_bins):
    [mi_bins ,  clrgr_bins] = get_mi_clrgr_bins()
    clr_lim = len(clrgr_bins)-1
    mag_lim = len(mi_bins)-1
    mag_i = clr_mag_bins[0]
    clr_i = clr_mag_bins[1]
    if(clr_i>0 and clr_i<clr_lim and mag_i > 0 and mag_i < mag_lim):
        return True
    else:
        return False

def make_background_estimates(query_results,background_folder,
                              galaxy_type,galaxy_weight,
                              zbins,
                              force=False,
                              use_all=True,
                              use_num=1000):
    background_file = background_file_name(background_folder,galaxy_type,galaxy_weight)
    print background_file
    if(not os.path.isfile(background_file) or force):
        print "Calculating new background estimates..."
        if(use_all):
            [dataclstr,datagal,rnd_mask_pass,rnd_num]= get_clstr("rnd",query_results,-1, till_end=True)
        else:
            [dataclstr,datagal,rnd_mask_pass,rnd_num]= get_clstr("rnd",query_results,use_num)
        backgrounds_counts = np.zeros_like(zbins)
        [m_i_bins,clr_gr_bins] =  get_color_mag_bins()
        [H,_,_]=np.histogram2d([],[],bins=(m_i_bins,clr_gr_bins))
        # make a 2d histogram space for each redshift bin
        H_z_backgrounds = np.zeros((len(zbins),H[:,0].size,H[0,:].size))
        for i in range(0,len(zbins)):
            z = zbins[i]
            H_background = np.zeros_like(H)
            counts = []
            used_number = 0
            for j in range(0,rnd_num):
                #if there is a mask in the area, don't use that area
                if(rnd_mask_pass == False): 
                    continue 
                used_number +=1
                #print j,len(datagal), datagal[j]['m_i'].shape
                #a = datagal[j]['m_i']
                #b = datagal[j]['m_i_err']
                #c = datagal[j]['clr_g-r']
                #d = datagal[j]['clr_g-r_err']
                slct = select_gal_old(z,datagal[j]['m_i'],datagal[j]['m_i_err'],datagal[j]['clr_g-r'],datagal[j]['clr_g-r_err'])
                weights = weigh_galaxies(dataclstr[j],datagal[j],galaxy_weight)
                counts.append(np.sum(weights[slct])/dataclstr[j]['sq_deg'])
                #change the z-dep features of the galaxies (poor man's k-correction, etc)
                set_datagal_z_dep(z,datagal[j])
                [h,_,_]=histogram2d_wrapper(datagal[j]['mag'][slct], datagal[j]['clr'][slct],
                                            weights=weights[slct], bins=(m_i_bins,clr_gr_bins))
                if(np.abs(1.0-np.sum(weights[slct])/np.sum(h)) > 1e-3):
                    print "cnt: %d, h: %d diff:%f"%(np.sum(weights[slct]),np.sum(h),np.abs(1.0-np.sum(weights[slct])/np.sum(h)))
                H_background+=h/dataclstr[j]['sq_deg']/get_color_mag_bin_area()
                #print "[%d] z:%f counts: %f counts/sqdeg: %f counts/sqkpc: %f"%(j,clstr_z[j],tmp,tmp2,tmp3)
            H_z_backgrounds[i,:,:]=H_background/float(used_number)
            backgrounds_counts[i] =np.average(counts)
            print "z=%.2f, cnt: %f, H_z: %f "%(z,backgrounds_counts[i],np.sum(H_z_backgrounds[i])*get_color_mag_bin_area())
            if(False):
                plt.figure()
                plt.pcolor(m_i_bins,clr_gr_bins,H_z_backgrounds[i].T,cmap=plt.cm.BuPu,norm=LogNorm())
                plt.axvline(1)
                plt.colorbar()
                plt.grid()
                plt.show()
        dtk.ensure_dir(background_folder)
        np.savez(background_file,z=zbins,cnt=backgrounds_counts,H=H_z_backgrounds)
        print "done."
    return 

def get_background_estimate_sqdeg(background_folder,galaxy_type,galaxy_weight):
    """Backgrounds is a funciton of redshift and in units galaxies/sq deg"""
    from scipy.interpolate import interp1d
    fname = background_file_name(background_folder,galaxy_type,galaxy_weight)
    npzfile = np.load(fname)
    z       = npzfile['z']
    cnt     = npzfile['cnt'] # per sq deg
    H_sqdeg = npzfile['H']  # per sq deg, per mag^2
    cnt_inter= scipy.interpolate.interp1d(z,cnt,kind='cubic')
    #[m_i_avg,clr_gr_avg]= get_color_mag_bin_avgs()
    #H_inter = interpol_H(z,m_i_avg,clr_gr_avg,H)
    H_inter = H_Interpol_fast(H_sqdeg,z)
    return [cnt_inter,H_inter]

def get_background_estimate_sqkpc(background_folder,galaxy_type,galaxy_weight):
    """get background as a funciton of redshift and in unites galaxies/sq kpc"""
    from scipy.interpolate import interp1d
    fname = background_file_name(background_folder,galaxy_type,galaxy_weight)
    npzfile = np.load(fname)
    z = npzfile['z']
    cnt=npzfile['cnt'] #/sq deg
    H  = npzfile['H']  #/sq deg, per mag^2
    
    sqdeg2sqkpc = 1.0/3600.0/(cosmo.kpc_proper_per_arcmin(z).value)**2
    cnt_sqkpc = cnt*sqdeg2sqkpc
    H_sqkpc = np.zeros_like(H)
    for i in range(0,len(z)):
        H_sqkpc[i,:,:] = H[i,:,:]*sqdeg2sqkpc[i]
    cnt_inter= scipy.interpolate.interp1d(z,cnt_sqkpc,kind='cubic')
    H_inter = H_Interpol_fast(H_sqkpc,z)
    zz=[0.05]
    #print "z=%.2f, cnt_inter: %e, H_inter %e" %(zz[0],cnt_inter(zz[0]),np.sum(H_inter(zz)*get_color_mag_bin_area()))
    #exit()
    return [cnt_inter,H_inter]

def background_file_name(bkgnd_folder,gal_type,gal_weight):
    return "%stype%d_weight%d_clr%d_mag%d_mstar%f_background.npz"%(bkgnd_folder,gal_type,gal_weight,_gal_clr_type,_gal_mag_type,_mstar_cut)


def save_error(file_base,h_mag_clr_clr_err,h_mag_clr_mag_err,mi_bins,gr_bins):
    fname = error_file_name(file_base)
    np.savez(fname,
             mi_bins = mi_bins,
             clrgr_bins = gr_bins,
             clr_err = h_mag_clr_clr_err,
             mag_err = h_mag_clr_mag_err);
    return 

def load_error(file_base):
    fname = error_file_name(file_base)
    return np.load(fname)

def error_file_name(bkgnd_folder): 
    return "%stype%d_weight%d_mag%d_clr%d_error.npz"%(bkgnd_folder,_gal_type,_gal_weight,_gal_mag_type,_gal_clr_type)

    
_gal_type = None
_gal_weight = None
def set_gal_type_weight(gal_type,gal_weight):
    global _gal_type
    global _gal_weight 
    _gal_type=gal_type
    _gal_weight = gal_weight

_gal_mag_type = None
_gal_mag_lim  = None
def set_gal_mag_type(mag_type,mag_lim):
    global _gal_mag_type
    global _gal_mag_lim
    _gal_mag_type =mag_type
    _gal_mag_lim  =mag_lim

_gal_clr_type = None
_gal_clr_lim  = None
def set_gal_clr_type(clr_type,clr_lim):
    global _gal_clr_type
    global _gal_clr_lim
    _gal_clr_type =clr_type
    _gal_clr_lim  =clr_lim

def get_mi_clrgr_bins():
    mi_bins = np.linspace(12,22,101)
    clrgr_bins = np.linspace(-1,3,41)
    return [mi_bins, clrgr_bins]

def get_color_mag_bins():
    mag = _gal_mag_lim
    clr = _gal_clr_lim
    m_bins    = np.linspace(mag[0],mag[1],np.int(mag[2]))
    clr_bins = np.linspace(clr[0],clr[1],np.int(clr[2]))
    return [m_bins,clr_bins]

def get_color_mag_bin_avgs():
    [m_i,clr]= get_color_mag_bins()
    mi_avg = (m_i[:-1]+m_i[1:])/2.0
    clr_avg = (clr[:-1]+clr[1:])/2.0
    return [mi_avg,clr_avg]

def get_color_mag_bin_area():
    # mag/clr bins must be linearily spaced
    [m_i,clr]= get_color_mag_bins()
    del_m = m_i[1]-m_i[0]    
    del_clr =clr[1]-clr[0]
    return del_m*del_clr

def interpol_H(z,mi,clr,H):
    result = {}
    for mi_i in range(0,len(mi)):
        for clr_i in range(0,len(clr)):
            data = H[:,mi_i,clr_i]
            result[mi_i,clr_i]=scipy.interpolate.interp1d(z,data,kind='linear')
    
    return H_Interpol(result)

class H_Interpol():

    def __init__(self,interpol):
        self.interpol = interpol
        [self.m_i,self.clr]=get_color_mag_bin_avgs()

    def __getitem__(self,z):
        print "interpolating z:",z
        [mi_avg,clr_avg]= get_color_mag_bin_avgs()
        H_of_z = np.zeros((len(mi_avg),len(clr_avg)))
        for m_i in range(0,len(mi_avg)):
            for c_i in range(0,len(clr_avg)):
                H_of_z[m_i,c_i]=self.interpol[m_i,c_i](z)
        return H_of_z

class H_Interpol_fast():
    def __init__(self,h_data,z_bins):
        self.h_data   = h_data 
        self.z_bins = z_bins
    def __call__(self,z):

        z_i = np.digitize(np.atleast_1d(z),self.z_bins)-1
        del_z = self.z_bins[z_i+1]-self.z_bins[z_i]
        z_1_weight = (self.z_bins[z_i+1]-z)/del_z
        z_2_weight = (z-self.z_bins[z_i])/del_z
        h_1 = self.h_data[z_i][0]*z_1_weight
        h_2 = self.h_data[z_i+1][0]*z_2_weight
        #print "z:",z
        #print "z_i",z_i,self.z_bins[z_i]
        #print "z_i+1",z_i+1,self.z_bins[z_i+1]
        #print "weights: ",z_1_weight,z_2_weight,z_1_weight+z_2_weight
        #print "h1,h2 shapes: ",h_1.shape,h_2.shape
        return (h_1+h_2)

def load_spt_clusters(spt_file):
    data = np.loadtxt(spt_file)
    # hubble ~ 0.7
    mass = data[:,0]*0.7
    z = data[:,1]
    ngal200 = data[:,2]
    ngalRS200 = data[:,3]
    
    mass = 1e14*mass
    result = {}
    result['m200']=mass
    result['z']=z
    result['ngal200']=ngal200
    result['ngalRS200']=ngalRS200
    return result

def result_file_name(result_folder,galaxy_type,galaxy_weight,extension='npz'):
    return "%stype%d_weight%d_mag%d_clr%d_result.%s"%(result_folder,galaxy_type,galaxy_weight,_gal_mag_type,_gal_clr_type,extension)

def save_radial_profile(base_folder, 
                        gal_type, 
                        gal_weight,
                        z_bins, 
                        mass_bins, 
                        z_mass_cnt,
                        rad_bins,
                        r200_rescale,
                        h_zm_rad_clr_mg,
                        h_zm_rad_clr_mg_cnt,
                        h_bkgd_Ngal,
                        h_bkgd_Ngal_err,
                        h_bkgd_Ngal_var,
                        h_radial_profile,
                        h_radial_profile_err,
                        h_radial_profile_var):
    dtk.ensure_dir(base_folder)
    fname = result_file_name(base_folder,gal_type,gal_weight,"npz")
    [mg_bins,clr_bins ]= get_color_mag_bins()
    zm_shp = (len(z_bins)-1,len(mass_bins)-1) # the zm matrix shape
    zmr_shp = (len(z_bins)-1,len(mass_bins)-1,len(rad_bins)-1)
    zm_cnt = np.zeros(zm_shp)
    Ngal     = np.zeros(zm_shp)
    Ngal_err = np.zeros(zm_shp)
    Ngal_var = np.zeros(zm_shp)
    rad_prof     = np.zeros(zmr_shp)
    rad_prof_err = np.zeros(zmr_shp)
    rad_prof_var = np.zeros(zmr_shp)
    zmr_cnt = np.zeros((len(z_bins)-1,len(mass_bins)-1,len(rad_bins)-1))
    zmr_clr_mg = np.zeros((len(z_bins)-1,len(mass_bins)-1,len(rad_bins)-1,len(mg_bins)-1,len(clr_bins)-1))
    for i in range(0,len(z_bins)-1):
        for j in range(0,len(mass_bins)-1):
            zm_cnt[i,j]=z_mass_cnt[i,j]
            Ngal[i,j] = h_bkgd_Ngal[i,j]
            Ngal_err[i,j] = h_bkgd_Ngal_err[i,j]
            Ngal_var[i,j] = h_bkgd_Ngal_var[i,j]
            for k in range(0,len(rad_bins)-1):
                zmr_cnt[i,j,k]    = h_zm_rad_clr_mg_cnt[i,j,k]
                zmr_clr_mg[i,j,k] = h_zm_rad_clr_mg[i,j,k]
                rad_prof[i,j,k]     = h_radial_profile[i,j][k]
                rad_prof_err[i,j,k] = h_radial_profile_err[i,j][k]
                rad_prof_var[i,j,k] = h_radial_profile_var[i,j,k]
    np.savez(fname,
             z_bins      = z_bins,
             mass_bins   = mass_bins,
             zm_cnt      = zm_cnt,
             rad_bins    = rad_bins,
             r200_rescale= r200_rescale,
             zmr_clr_mg  = zmr_clr_mg,
             zmr_cnt     = zmr_cnt,
             clr_bins    = clr_bins,
             mg_bins     = mg_bins,
             Ngal        = Ngal,
             Ngal_err    = Ngal_err,
             Ngal_var    = Ngal_var,
             rad_prof    = rad_prof,
             rad_prof_err= rad_prof_err,
             rad_prof_var= rad_prof_var)
    save_zmr_hdf5(base_folder,gal_type,gal_weight,
                  z_bins = z_bins,
                  m_bins = mass_bins,
                  r_bins = rad_bins,
                  zm_Ngal = Ngal,
                  zm_Ngal_err = Ngal_err,
                  zm_Ngal_var = Ngal_var,
                  zmr_gal_density = rad_prof,
                  zmr_gal_density_err = rad_prof_err,
                  zmr_gal_density_var = rad_prof_var,
                  zm_counts           = zm_cnt,
                  zmr_counts          = zmr_cnt)
                  
    #load_radial_profile(base_folder,gal_type,gal_weight)

def save_zmr_hdf5(result_folder,gal_type,gal_weight,
                  z_bins,
                  m_bins,
                  r_bins,
                  zm_Ngal,
                  zm_Ngal_err,
                  zm_Ngal_var,
                  zmr_gal_density,
                  zmr_gal_density_err,
                  zmr_gal_density_var,
                  zm_counts,
                  zmr_counts):
    hfname = result_file_name(result_folder,gal_type,gal_weight,'hdf5')
    hf = h5py.File(hfname,'w')
    hf.create_dataset('z_bins',data=z_bins)
    hf.create_dataset('m_bins',data=m_bins)
    hf.create_dataset('r_bins',data=r_bins)
    hf.create_dataset('zm_Ngal',data=zm_Ngal)
    hf.create_dataset('zm_Ngal_err',data=zm_Ngal_err)
    hf.create_dataset('zm_Ngal_var',data=zm_Ngal_var)
    hf.create_dataset('zmr_gal_density',data=zmr_gal_density)
    hf.create_dataset('zmr_gal_density_err',data=zmr_gal_density_err)
    hf.create_dataset('zmr_gal_density_var',data=zmr_gal_density_var)
    hf.create_dataset('zm_counts',data=zm_counts)
    hf.create_dataset('zmr_counts',data=zmr_counts)
    hf.flush()
    hf.close()
    return

def load_zmr_hdf5():
    hfname = result_file_name(result_folder,gal_type,gal_weight,'hdf5')
    hf = hdf5.file(hfname,'r')
    return hf

def load_radial_profile(base_folder,gal_type,gal_weight):
    fname = result_file_name(base_folder,gal_type,gal_weight)
    npzfile = np.load(fname)
    return npzfile # use the named variables in save radial_profile such as npz['mass_bins']=mass_bins
    

    
def print_out_clr_mg (clr_mg):
    [clr,mg]= get_color_mag_bins()
    for i in range(0,len(mg)-1):
        for j in range(0,len(clr)-1):
            print "%1.1f"%clr_mg[j,i],
        print "\n"
    print "\n==========\n\n"
    return


def red_squence_select(z1,z2,h):
    mg_avg,clr_avg = get_color_mag_bin_avgs()
    mg_bins,clr_bins = get_color_mag_bins()
    rs_clr1 = redsequence_line(mg_avg,z1)
    rs_clr2 = redsequence_line(mg_avg,z2)
    slct = np.zeros_like(h)
    #print rs_clr1
    #print rs_clr2
    for i in range(0,len(mg_avg)):
        for j in range(0,len(clr_avg)-1):
            clr = clr_avg[j]
            rs1 = rs_clr1[i] #lower clr
            rs2 = rs_clr2[i] #higher clr
            upper = clr_avg[j+1]
            lower = clr_avg[j]
            if(upper < rs2 and lower > rs1):
                slct[i,j]=1#slct[i,j]
            elif(upper > rs2 and lower < rs2):
                del_clr = upper-lower
                del_rs = rs2-lower
                slct[i,j] = del_rs/del_clr
            elif(upper > rs1 and lower < rs1):
                del_clr = upper-lower
                del_rs =  upper-rs1
                slct[i,j] = del_rs/del_clr
    return slct
    
def non_red_sqeuence_select(z1,z2,h):
    slct = np.ones_like(h)-red_squence_select(z1,z2,h)
    return slct

def red_squence_redder_select(z1,z2,h):
    mg_avg,clr_avg = get_color_mag_bin_avgs()
    mg_bins,clr_bins = get_color_mag_bins()
    rs_clr1 = redsequence_line(mg_avg,z1)
    rs_clr2 = redsequence_line(mg_avg,z2)
    slct = np.zeros_like(h)
    #print rs_clr1
    #print rs_clr2
    for i in range(0,len(mg_avg)):
        for j in range(0,len(clr_avg)-1):
            clr = clr_avg[j]
            rs1 = rs_clr1[i] #lower clr
            upper = clr_avg[j+1]
            lower = clr_avg[j]
            if(lower > rs1):
                slct[i,j]=1#slct[i,j]
            elif(upper > rs1 and lower < rs1):
                del_clr = upper-lower
                del_rs =  upper-rs1
                slct[i,j] = del_rs/del_clr
    return slct

def non_red_squence_redder_select(z1,z2,h):
    slct = np.ones_like(h)-red_squence_redder_select(z1,z2,h)
    return slct
