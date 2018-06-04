#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.cosmology import WMAP9 as cosmo


def get_clstr(name,num,start=0):
    clstr_i=[]
    clstr_i_err =[]
    clstr_clr = []
    clstr_clr_err = []
    clstr_m200 = []
    clstr_r200 = []
    clstr_sq_deg = []
    clstr_gal_cnt =[]
    clstr_z = []
    clstr_rad = []
    for i in range(start,num):
        datagal = np.genfromtxt("sdss_pulls/%s%d.txt"%(name,i),names=True,delimiter=',',dtype="f4,f4,i4,i1,i8,i8,i8,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4")
        dataclust = np.genfromtxt("sdss_pulls/%s_prop%d.txt"%(name,i),names=True,delimiter=',')
                # flags defined at http://skyserver.sdss.org/dr8/en/help/browser/enum.asp?n=PhotoFlags
        #          bright               saturated             satur_center         nopetro             deblended_as_moving
        cut_flgs = 0x0000000000000002 | 0x0000000000040000  | 0x0000080000000000 | 0x0000000000000100 | 0x0000000100000000
        flgs = datagal['flags_i'] | datagal['flags_r'] | datagal['flags_g']
        slct1 = datagal['type']==3
        slct2 = np.logical_not(flgs&cut_flgs)
        slct = slct1 & slct2
        a_i = datagal['cModelMag_i'][slct]
        a_clr = datagal['cModelMag_g'][slct]-datagal['cModelMag_r'][slct]
        a_i_err = datagal['cModelMagErr_i'][slct]
        a_clr_err = np.sqrt(datagal['cModelMagErr_r'][slct]**2+datagal['cModelMagErr_g'][slct]**2)
        a_rad = rad_from_ra_dec(dataclust['ra'],dataclust['dec'],datagal['ra'],datagal['dec'])
        a_m200 =dataclust['mass']
        a_sq_deg = np.pi*dataclust['rad_arcmin']**2/3600.0
        a_r200 = dataclust['r200']
        clstr_i.append(a_i)
        clstr_i_err.append(a_i_err)
        clstr_clr.append(a_clr)
        clstr_clr_err.append(a_clr_err)
        clstr_rad.append(a_rad)
        
        clstr_m200.append(a_m200)
        clstr_r200.append(a_r200)
        clstr_sq_deg.append(a_sq_deg)
        clstr_gal_cnt.append(-1.0)
        clstr_z.append(dataclust['z'])
        
    return [np.array(clstr_i),np.array(clstr_i_err),np.array(clstr_clr),np.array(clstr_clr_err),np.array(clstr_rad),np.array(clstr_m200),np.array(clstr_r200),np.array(clstr_sq_deg),np.array(clstr_gal_cnt),np.array(clstr_z)]

def kpc_from_arcmin(arcmin,z):
    return arcmin*cosmo.kpc_proper_per_arcmin(z).value

def rad_from_ra_dec(ra1,dec1,ra2,dec2):
    mean_dec = (dec1+dec2)/2.0
    del_dec = dec1-dec2
    del_ra = (ra1-ra2)*np.cos(mean_dec)
    dist = np.sqrt(del_dec**2+del_ra**2)
    return dist

def deposite(xi,yi,val,result,num_x,num_y):
    if(xi >= 0 and xi < num_x and yi >= 0 and yi < num_y):
        result[xi,yi]+=val

def cell_in_cloud(datax,datay,bins=None):
    min_x = bins[0].min()
    max_x = bins[0].max()
    min_y = bins[1].min()
    max_y = bins[1].max()
    num_x = bins[0].size-1
    num_y = bins[1].size-1
    del_x = (max_x-min_x)/num_x
    del_y = (max_y-min_y)/num_y
    A = del_x*del_y
    result = np.zeros([num_x,num_y])
    for i in range(0,len(datax)):
        xi = int((datax[i]-min_x)/del_x)
        yi = int((datay[i]-min_y)/del_y)
        a1=  (datax[i]-(xi*del_x+min_x))*(datay[i]-(yi*del_y+min_y))/A
        a2=  (datax[i]-(xi*del_x+min_x))*(-datay[i]+((yi+1)*del_y+min_y))/A
        a3= (-datax[i]+((xi+1)*del_x+min_x))*(datay[i]-(yi*del_y+min_y))/A
        a4= (-datax[i]+((xi+1)*del_x+min_x))*(-datay[i]+((yi+1)*del_y+min_y))/A
        deposite(xi,yi,a1,     result,num_x,num_y)
        deposite(xi,yi+1,a2,   result,num_x,num_y)
        deposite(xi+1,yi,a3,   result,num_x,num_y)
        deposite(xi+1,yi+1,a4, result,num_x,num_y)
    return [result,bins[0],bins[1]]

def mstar(z):
    return 12.27+62.36*z-289.79*z**2+729.69*z**3-709.42*z**4

def m_cut(z):
    return mstar(z)+1.0

backgrounds_redshifts = np.linspace(0.05,0.35,30)

def select_background(z,gal_i,gal_i_err,gal_clr,gal_clr_err):
    m17 = -0.0701*z-0.008969
    b17 = 3.2982*z+0.5907
    mean_clr = (m17*(gal_i-17.0)+b17)
    clr_diff = mean_clr-gal_clr
    scatter =np.sqrt(0.15**2+gal_clr_err**2)
    slct2 = (np.abs(clr_diff) <= 3.0*scatter)
    slct1 = (gal_i<m_cut(z)) & (12.0<gal_i)
    slct = slct1 & slct2
    return slct

def make_background_estimates(backgrounds_redshifts,clstr_i,clstr_i_err,clstr_clr,clstr_clr_err,clstr_sq_deg):
    from scipy.interpolate import interp1d
    backgrounds_counts = np.zeros_like(backgrounds_redshifts)
    #print "clstr_i.shape: ",clstr_i.shape
    for i in range(0,backgrounds_counts.size):
        z = backgrounds_redshifts[i]
        counts = []
        for j in range(0,rnd_num):
            tmp = np.sum(select_background(z,clstr_i[j],clstr_i_err[j],clstr_clr[j],clstr_clr_err[j]))
            tmp2=tmp/clstr_sq_deg[j]
            counts.append(tmp2)
        backgrounds_counts[i] +=np.average(counts)
    return interp1d(backgrounds_redshifts,backgrounds_counts,kind='cubic')
                                       
def cluster_bckgnd_sub_rs(z,clstr_i,clstr_clr,clstr_sq_deg,H,xbins,ybins):
    expected = background_rs(z,H,xbins,ybins)*clstr_sq_deg
    clstr_raw = cluster_mcut_rs(z,clstr_i,clstr_clr)
    return clstr_raw-expected

def background_rs(z,H,xbins,ybins):
    min_i = 12
    max_i = m_cut(z)
    mask = background_mask(H,xbins,ybins,z,min_i,max_i)
    #print "background rs check, all should be the same: ",mask.shape,H.shape,(mask*H).shape
    return np.sum(mask*H)

def cluster_mcut_rs(z,gal_i,gal_i_err,gal_clr,gal_clr_err):
    mstar_cut =m_cut(z)
    m17 = -0.0701*z-0.008969
    b17 = 3.2982*z+0.5907
    mean_clr = (m17*(gal_i-17.0)+b17)
    gal_clr_diff = mean_clr-gal_clr
    scatter = np.sqrt(0.15**2+gal_clr_err**2)
    slct1 = (gal_i<m_cut(z)) & (12.0<gal_i)
    slct2 = np.abs(gal_clr_diff) <= 3.0*scatter
    return np.sum(slct1&slct2)
    

def mstar(z):
    ms = 12.27+62.36*z-289.79*z**2+729.69*z**3-709.42*z**4
    return ms

def cluster_gal_count(clstr_z,clstr_i,clstr_i_err,clstr_clr,clstr_clr_err):
    clstr_gal_cnt = np.zeros(clstr_num)
    for i in range(0,clstr_num):
        clstr_gal_cnt[i] = cluster_mcut_rs(clstr_z[i],
                                           clstr_i[i],
                                           clstr_i_err[i],
                                           clstr_clr[i],
                                           clstr_clr_err[i])
    return clstr_gal_cnt

        
clstr_num = 100
rnd_num =300
[clstr_i,clstr_i_err,clstr_clr,clstr_clr_err,clstr_rad,clstr_m200,clstr_r200,clstr_sq_deg,clstr_gal_cnt,clstr_z] = get_clstr("gal",clstr_num)
[rnd2_i,rnd2_i_err,rnd2_clr,rnd2_clr_err,rnd2_rad,rnd2_m200,rnd2_r200,rnd2_sq_deg,rnd2_gal_cnt,rnd_z] = get_clstr("rnd",rnd_num)
clstr_all_gal_cnt = np.zeros(len(clstr_i))

clstr_slct = (clstr_z < 0.35) & (clstr_z > 0.05)
bkgnd_per_sq_deg = make_background_estimates(backgrounds_redshifts,rnd2_i,rnd2_i_err,rnd2_clr,rnd2_clr_err,rnd2_sq_deg)
clstr_gal_count = cluster_gal_count(clstr_z,clstr_i,clstr_i_err,clstr_clr,clstr_clr_err)
clstr_bkgnd_count = bkgnd_per_sq_deg(clstr_z[clstr_slct])*clstr_sq_deg[clstr_slct]
NgalRS= clstr_gal_count[clstr_slct] - clstr_bkgnd_count
M200_ngal = clstr_m200[clstr_slct]
z_ngal    = clstr_z[clstr_slct]

def plot_background(z):
    H_sum,xbins,ybins=np.histogram2d([],[],bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))
    Hs = []
    for i in range(0,rnd_num):
        slct = select_background(z,rnd2_i[i],rnd2_i_err[i],rnd2_clr[i],rnd2_clr_err[i])
        #print "slct sum: %f, %f/sqdeg"%(np.sum(slct),np.sum(slct)/rnd2_sq_deg[i])
        H,xbins,ybins=np.histogram2d(rnd2_i[i][slct],rnd2_clr[i][slct],bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))
        H_sum +=H/rnd2_sq_deg[i]
    print "H_sum.shape: ", H_sum.shape
    plt.figure()
    [X,Y]=np.meshgrid(xbins,ybins)
    plt.pcolor(X,Y,H_sum.T,cmap=plt.cm.PuBu,norm=LogNorm())
    plt.colorbar()
    plt.title('z=%.1f background, sum = %.1f'%(z,H_sum.sum()/rnd_num))
    plt.xlabel('i-mag')
    plt.ylabel('g-r clr')

plt.figure()
plt.hist(clstr_rad.flatten())
    
plt.figure()
plt.plot(backgrounds_redshifts,bkgnd_per_sq_deg(backgrounds_redshifts))
plt.title("Bkgnd/sq deg at diff. redshifts")
plt.ylabel("bkgnd galaxy count/sq.deg.")
plt.xlabel('redshift')
plt.grid()

plt.figure()
plt.scatter(M200_ngal,NgalRS,c=z_ngal)
plt.colorbar().set_label('redshift')
plt.ylabel('NgalRS')
plt.xlabel('M200')
plt.xscale('log')
plt.grid()

plt.figure()
plt.scatter(M200_ngal,clstr_gal_count[clstr_slct],c=z_ngal)
plt.colorbar().set_label('redshift')
plt.ylabel('RS galaxy count')
plt.xlabel('M200')
plt.xscale('log')
plt.grid()

plt.figure()
plt.scatter(M200_ngal,clstr_bkgnd_count,c=z_ngal)
plt.colorbar().set_label('redshift')
plt.ylabel('estimated background')
plt.xlabel('M200')
plt.xscale('log')
plt.grid()

plot_background(0.1)
plot_background(0.2)
plot_background(0.3)


plt.show()
exit()
