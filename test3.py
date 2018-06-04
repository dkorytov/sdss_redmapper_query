#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def get(name,num,start=0):
    gal_i = []
    gal_clr = []
    ster_rad = 0.0
    sq_deg = np.array(0.0,dtype='f8')
    num_per_sq_arcmin = []
    for i in range(start,num):
        datagal = np.genfromtxt("sdss_pulls/%s%d.txt"%(name,i),names=True,delimiter=',',dtype="f4,f4,i4,i1,i8,i8,i8,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4")
        dataclust = np.genfromtxt("sdss_pulls/%s_prop%d.txt"%(name,i),names=True,delimiter=',')
        area = np.pi*dataclust['rad_arcmin']**2/3600.0
        obj_num = float(datagal['cModelMag_i'].size)
        #    print obj_num/area, obj_num,area,"      ",40/(obj_num/area)
        #datarnd = np.genfromtxt("output/rnd%d.txt"%i,names=True,delimiter=',')
        strange = np.argmin(datagal['cModelMag_u'])
        #print "u:%f clr:%f"%(datagal['cModelMag_u'][strange],datagal['cModelMag_g'][strange]-datagal['cModelMag_r'][strange])
        gal_i_i = datagal['cModelMag_i']
        gal_clr_i =datagal['cModelMag_g']-datagal['cModelMag_r']
        slct3 = datagal['type'] ==3
        flgs = datagal['flags_i'] | datagal['flags_r'] | datagal['flags_g']
        
        # flags defined at http://skyserver.sdss.org/dr8/en/help/browser/enum.asp?n=PhotoFlags
        #          bright               saturated             satur_center         nopetro             deblended_as_moving
        cut_flgs = 0x0000000000000002 | 0x0000000000040000  | 0x0000080000000000 | 0x0000000000000100 | 0x0000000100000000
        slct4 = np.logical_not(flgs&cut_flgs)
        #print slct4
        slct =  slct3 & slct4
        num_per_sq_arcmin.append(float(np.sum(slct))/area)
        #print "avg #gal/sq arcmin", np.average(num_per_sq_arcmin)
        
        gal_i = np.append(gal_i,gal_i_i[slct])
        gal_clr =  np.append(gal_clr,gal_clr_i[slct])
        #ster_rad += 2.0*np.pi*(1.0-np.cos(dataclust['rad_arcmin']*np.pi/(360.0*60.0)))
        sq_deg = np.append(sq_deg,np.pi*dataclust['rad_arcmin']**2/3600.0)
    return np.array(gal_i),np.array(gal_clr),np.array(ster_rad),np.array(sq_deg),np.array(num_per_sq_arcmin)

def get_clstr(name,num,start=0):
    clstr_i=[]
    clstr_clr = []
    clstr_m200 = []
    clstr_r200 = []
    clstr_sq_deg = []
    clstr_gal_cnt =[]
    clstr_z = []
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
        a_m200 =dataclust['mass']
        a_sq_deg = np.pi*dataclust['rad_arcmin']**2/3600.0
        a_r200 = dataclust['r200']
        clstr_i.append(a_i)
        clstr_clr.append(a_clr)
        clstr_m200.append(a_m200)
        clstr_r200.append(a_r200)
        clstr_sq_deg.append(a_sq_deg)
        clstr_gal_cnt.append(-1.0)
        clstr_z.append(dataclust['z'])
    return [np.array(clstr_i),np.array(clstr_clr),np.array(clstr_m200),np.array(clstr_r200),np.array(clstr_sq_deg),np.array(clstr_gal_cnt),np.array(clstr_z)]
        
clstr_num = 1000
rnd_num =3000
[rnd_i,rnd_clr,ster_rad,rnd_sq_deg,rnd_num_per_sq_arcmin] = get("rnd",rnd_num)
[clstr_i,clstr_clr,clstr_m200,clstr_r200,clstr_sq_deg,clstr_gal_cnt,clstr_z] = get_clstr("gal",clstr_num)
[rnd2_i,rnd2_clr,rnd2_m200,rnd2_r200,rnd2_sq_deg,rnd2_gal_cnt,rnd_z] = get_clstr("rnd",rnd_num)


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

def background(z,H2,xbins,ybins):
    mstar_cut = m_cut(z)
    x_select = np.searchsorted(xbins,mstar_cut)
    print H2[:x_select,:].shape
    sq_deg = np.sum(H2[:x_select,:]) #bins that are clearly above m*
    sq_deg_sliver = np.sum(H2[x_select:(x_select+1),:]) #the bin that contains
    #m*, we will take a linear cut
    sliver_fraction = (xbins[x_select]-mstar_cut)/(xbins[x_select]-xbins[x_select-1])
    #print "raw sum: ",np.sum(H2)
    #print "within sum: ", sq_deg
    #print "sq_deg_sliver: ",sq_deg_sliver
    #print "sliver_fraction: ",sliver_fraction
    #print "(%f-%f)/(%f-%f) = %f" % ( xbins[x_select],mstar_cut,xbins[x_select],xbins[x_select-1],sliver_fraction)

    return sq_deg + sq_deg_sliver*sliver_fraction

def check_background_mask(H,xbins,ybins,z):
    from matplotlib.colors import LogNorm
    min_i = 12
    max_i = m_cut(z)
    mask = background_mask(H,xbins,ybins,z,min_i,max_i)
    print "z: %f, i=%f->%f"%(z,min_i,max_i)
    [X,Y] = np.meshgrid(xbins,ybins)
    plt.figure()
    plt.pcolor(X,Y,mask.T,cmap=plt.cm.PuBu,norm=LogNorm())
    plt.title('Mask value')
    plt.figure()
    plt.pcolor(X,Y,(mask*H).T,cmap=plt.cm.PuBu,norm=LogNorm())
    plt.show()
    
def rs_color_cut(z):
    avg_clr = 0.625+3.149*z
    scatter = 0.05
    return [avg_clr-3.0*scatter,avg_clr+3.0*scatter]

def background_mask(H,i_bins,clr_bins,z,min_i,max_i):
    mask = np.ones_like(H)
    for i in range(0,i_bins.size-1):
        #too far out on either side
        if((xbins[i] > max_i) or (xbins[i+1]<min_i)):
           mask[i,:]*=0.0
        #if on the left edge, fractional weight 
        elif((xbins[i]<min_i) and (xbins[i+1]>min_i)):
            weight = (xbins[i+1]-min_i)/(xbins[i+1]-xbins[i])
            mask[i,:]*=weight
        #if on the right edge, fractional weight 
        elif((xbins[i]<max_i) and (xbins[i+1]>max_i)):
            weight = (max_i-xbins[i])/(xbins[i+1]-xbins[i])
            mask[i,:]*=weight
        m17 = -0.0701*z-0.008969
        b17 = 3.2982*z+0.5907
        mean_clr = (m17*(xbins[i]-17.0)+b17)
        min_clr = mean_clr-0.15
        max_clr = mean_clr+0.15
        #print "min_clr %f, max_clr %f"%(min_clr,max_clr)
        for j in range(0,clr_bins.size-1):
            if((ybins[j]>max_clr ) or (ybins[j+1]<min_clr)):
                mask[i,j]*=0.0
                #print "clr[%f-%f]: %f"%(ybins[j],ybins[j+1],0)
            elif((ybins[j]<min_clr) and(ybins[j+1]>min_clr)):
                weight = (ybins[j+1]-min_clr)/(ybins[j+1]-ybins[j])
                mask[i,j]*=weight
                #print "clr[%f-%f]: %f"%(ybins[j],ybins[j+1],weight)
            elif((ybins[j]<max_clr) and (ybins[j+1] > max_clr)):
                weight = (max_clr - ybins[j])/(ybins[j+1]-ybins[j])
                mask[i,j]*=weight
                #print "clr[%f-%f]: %f"%(ybins[j],ybins[j+1],weight)
        #print "\n\n"
    return mask

def cluster_mcut(z,clstr_i,clstr_clr):
    mstar_cut =m_cut(z)
    slct1 = (clstr_i<m_cut(z)) & (12.0<clstr_i)
    slct2 = (clstr_clr<3.0) & (-1.0<clstr_clr)
    return np.sum(slct1&slct2)
    
def cluster_bckgnd_sub(z,clstr_i,clstr_clr,clstr_sq_deg,H,xbins,ybins):
    expected = background(z,H,xbins,ybins)*clstr_sq_deg
    clstr_raw = cluster_mcut(z,clstr_i,clstr_clr)
    return clstr_raw-expected

def cluster_bckgnd_sub_rs(z,clstr_i,clstr_clr,clstr_sq_deg,H,xbins,ybins):
    expected = background_rs(z,H,xbins,ybins)*clstr_sq_deg
    clstr_raw = cluster_mcut_rs(z,clstr_i,clstr_clr)
    return clstr_raw-expected

def background_rs(z,H,xbins,ybins):
    [min_clr,max_clr] = rs_color_cut(z)
    min_i = 12
    max_i = m_cut(z)
    mask = background_mask(H,xbins,ybins,z,min_i,max_i)
    #print "background rs check, all should be the same: ",mask.shape,H.shape,(mask*H).shape
    return np.sum(mask*H)

def cluster_mcut_rs(z,clstr_i,clstr_clr):
    mstar_cut =m_cut(z)
    m17 = -0.0701*z-0.008969
    b17 = 3.2982*z+0.5907
    mean_clr = (m17*(clstr_i-17.0)+b17)
    clstr_clr_diff = mean_clr-clstr_clr
    max_diff = 0.15
    min_diff = -0.15
    slct1 = (clstr_i<m_cut(z)) & (12.0<clstr_i)
    [min_clr,max_clr] = rs_color_cut(z)
    slct2 = (clstr_clr_diff > min_diff) & (clstr_clr_diff < max_diff)
    return np.sum(slct1&slct2)
    
[rnd_H2,xbins,ybins] = cell_in_cloud(rnd_i,rnd_clr,bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))
print "shapes: ", rnd_H2.shape,xbins.shape,ybins.shape
rnd_H = rnd_H2/np.sum(rnd_sq_deg) #the x100 is to get /mag^2 bins
clstr_all_gal_cnt = np.zeros(len(clstr_i))
for i in range(0,len(clstr_i)):
    [clust_H,_,_] = cell_in_cloud(clstr_i[i],clstr_clr[i],bins=(xbins,ybins))

    clust_H/clstr_sq_deg[i]
    diff = clust_H-rnd_H
    [X,Y]=np.meshgrid(xbins,ybins)
    #  plt.figure()
    #  plt.pcolor(X,Y,diff.T.clip(min=0.0),cmap=plt.cm.PuBu)
    #  plt.colorbar()
    #  plt.ylabel('g-r')
    #  plt.xlabel('i')
    #  plt.title('gal diff /sq deg')
    #  plt.show()
    clstr_gal_cnt[i] = np.sum(diff[diff>0.0])*clstr_sq_deg[i]
    clstr_all_gal_cnt[i] = np.sum(diff)*clstr_sq_deg[i]
    # print "m200: %e sq_deg: %f gal cnt: %f gal all cnt: %f"%(clstr_m200[i],clstr_sq_deg[i],clstr_gal_cnt[i],    clstr_all_gal_cnt[i])

def mstar(z):
    ms = 12.27+62.36*z-289.79*z**2+729.69*z**3-709.42*z**4
    return ms



plt.figure()
plt.plot(clstr_m200,clstr_gal_cnt,'.')
plt.ylabel('gal cnt')
plt.xlabel('m200 [Msun/h]')
plt.yscale('log')
plt.xscale('log')

plt.figure()
plt.plot(clstr_m200,clstr_all_gal_cnt,'.')
plt.ylabel('gal cnt all')
plt.xlabel('m200 [Msun/h]')
plt.xscale('log')


plt.figure()
plt.plot(clstr_sq_deg,clstr_all_gal_cnt,'.')
plt.ylabel('gal cnt all')
plt.xlabel('sq deg')
plt.xscale('log')

clstr_gal_raw = np.array([l.size for l in clstr_i])
rnd2_gal_raw = np.array([l.size for l in rnd2_i])

rnd2_per_sq_deg = np.sum(rnd2_gal_raw)/np.sum(rnd2_sq_deg)
clstr_per_sq_deg = np.sum(clstr_gal_raw)/np.sum(clstr_sq_deg)
clstr_slct = np.zeros(clstr_num)
rnd_slct   = np.zeros(rnd_num)
clstr_correct = np.zeros(clstr_num)
for i in range(0,rnd2_i.size):
    slct = (rnd2_i[i]<22) & (rnd2_i[i]>12) & (rnd2_clr[i] > -1) & (rnd2_clr[i]<3)
    rnd_slct[i] = np.sum(slct)
for i in range(0,clstr_i.size):
    slct = (clstr_i[i]<22) & (clstr_i[i]>12) & (clstr_clr[i] > -1) & (clstr_clr[i]<3)
    clstr_slct[i] = np.sum(slct)
    clstr_correct[i] = cluster_bckgnd_sub_rs(clstr_z[i],clstr_i[i],clstr_clr[i],clstr_sq_deg[i],rnd_H,xbins,ybins)
    if clstr_z[i]>0.35:
        clstr_correct[i]=np.nan
plt.figure()
plt.plot(clstr_sq_deg,clstr_gal_raw,'.',label='redmapper cluster')
plt.plot(rnd2_sq_deg,rnd2_gal_raw,'.',label='redmapper randoms')
plt.plot([0,0.5],[0,0.5*clstr_per_sq_deg],'-',label='expected clstr')
plt.plot([0,0.5],[0,0.5*rnd2_per_sq_deg],'-',label='expected rnd')
plt.xlabel('sq.deg.')
plt.ylabel('galaxy count')
plt.legend(loc='best')

plt.figure()
plt.plot(rnd2_sq_deg,rnd2_gal_raw-rnd2_sq_deg*rnd2_per_sq_deg,'.',label='redmapper randoms')
plt.plot(clstr_sq_deg,clstr_gal_raw-rnd2_per_sq_deg*clstr_sq_deg,'.',label='redmapper cluster')
plt.legend(loc='best')

print rnd2_sq_deg.size
print rnd_slct.size
plt.figure()
plt.plot(rnd2_sq_deg,rnd_slct,'.',label='redmapper randoms')
plt.plot(clstr_sq_deg,clstr_slct,'.',label='redmapper clusters')
plt.plot([0.0,0.5],[0,np.sum(rnd_slct)/np.sum(rnd2_sq_deg)*0.5],'-',label='randoms line')
plt.title('selected in i-band & g-r color')
plt.ylabel('gal num')
plt.xlabel('sq deg')

background_per_sq_deg = np.sum(rnd_slct)/np.sum(rnd2_sq_deg)

plt.figure()
plt.scatter(clstr_m200,clstr_slct-background_per_sq_deg*clstr_sq_deg,c=clstr_z)
plt.title("Ngal w/ background subtraction")
plt.colorbar().set_label('redshift')
plt.ylabel('Ngal200')
plt.xlabel('m200')
plt.xscale('log')
plt.grid()

plt.figure()
plt.scatter(clstr_m200,clstr_correct,c=clstr_z)
plt.colorbar().set_label('redshift')
plt.ylabel('RS Ngal,200')
plt.xlabel('M200')
plt.xscale('log')
plt.title("RS Ngal")
plt.grid()

plt.show()

check_background_mask(rnd_H,xbins,ybins,0.1)
