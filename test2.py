#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

expected = np.genfromtxt("expected.csv",delimiter=',')

print expected
def get(name,num,start=0):
    gal_i = []
    gal_clr = []
    ster_rad = 0.0
    sq_deg = np.array(0.0,dtype='f8')
    num_per_sq_arcmin = []
    for i in range(start,num):
 #       datagal = np.load("sdss_pulls/%s%d.npy"%(name,i))
        datagal = np.genfromtxt("sdss_pulls/%s%d.txt"%(name,i),names=True,delimiter=',',dtype="f4,f4,i4,i1,i8,i8,i8,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4")
#        dataclust = np.load("sdss_pulls/%s_prop%d.npy"%(name,i))
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
    return gal_i,gal_clr,ster_rad,sq_deg,num_per_sq_arcmin

def get_clstr(name,num,start=0):
    clstr_i=[]
    clstr_clr = []
    clstr_m200 = []
    clstr_r200 = []
    clstr_sq_deg = []
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
        clstr_i.append(a_i)
        clstr_clr.append(a_clr)

    return [clstr_i,clstr_clr,clstr_m200,clstr_r200,clstr_sq_deg]
        
        
[gal_i,gal_clr,ster_rad,sq_deg,num_per_sq_arcmin] = get("gal",1000)
[rnd_i,rnd_clr,ster_rad,rnd_sq_deg,rnd_num_per_sq_arcmin] = get("rnd",3000)
[clstr_i,clstr_clr,clstr_m200,clstr_r200,clstr_sq_deg] = get_clstr("gal",1000)
plt.figure()
plt.plot(gal_i,gal_clr,'.',alpha=0.01)

def smoothen(data,xsize,ysize):
    result = np.zeros(data.shape)
    for i in range(-1,1):
        for j in range(-1,1):
            print i, j
            result[0+max(i,0):xsize+min(i,0),0+max(j,0):ysize+min(j,0)] += data[0-min(i,0):xsize-max(i,0),0-min(j,0):ysize-max(j,0)]
    result[0,:] += 3.0*data[0,:]
    result[-1,:] += 3.0*data[-1,:]
    result[1:,0] += 3.0*data[1:,0]
    result[1:,-1] += 3.0*data[1:,-1]
    return result/9.0

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
#        print i, a1,a2,a3,a4,result.sum()
    return result

print gal_i.size, gal_clr.size
print rnd_i.size, rnd_clr.size
plt.figure()
print "hist 1"
red_H,xbins,ybins=np.histogram2d(gal_i,gal_clr,bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))#,np.linspace(10,30,40),np.linspace(-1,1,40))
rnd_H,xbins,ybins=np.histogram2d(rnd_i,rnd_clr,bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))#,np.linspace(10,30,40),np.linspace(-1,1,40))
print "bin shaples: ", np.linspace(12,22,101).shape, np.linspace(12,22,101).min(),np.linspace(12,22,101).max(),xbins.shape,xbins.min(),xbins.max()
print "cic"
red_H2 = cell_in_cloud(gal_i,gal_clr,bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))
rnd_H2 = cell_in_cloud(rnd_i,rnd_clr,bins=(np.linspace(12,22,101),np.linspace(-1,3,41)))
print rnd_H2.sum(), rnd_H.sum()

#H = smoothen(H,100,40)
print "sq deg rnd:%f clst:%f"%(np.sum(sq_deg),np.sum(rnd_sq_deg))
print "avg clstr: %f"%np.average(sq_deg)
red_H = red_H2/np.sum(sq_deg)*10.0*10.0
rnd_H = rnd_H2/np.sum(rnd_sq_deg)*10.0*10.0
np.savetxt("background.txt",rnd_H.T)
print red_H.sum(),rnd_H.sum(),(red_H.sum()-rnd_H.sum())*np.sum(sq_deg)/100.0
#red_H = red_H/np.sum(red_H)
#rnd_H = rnd_H/np.sum(rnd_H)
expected = expected
val_min = 1e-1#rnd_H[rnd_H>0].min()
val_max = max([red_H.max(), expected.max()])
print red_H.sum(), rnd_H.sum(), expected.sum()
[X,Y]=np.meshgrid(xbins,ybins)
plt.pcolor(X,Y,rnd_H.T,cmap=plt.cm.PuBu,norm=LogNorm(vmin=val_min,vmax=val_max))
plt.colorbar()
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('SDSS RedMapper Random Galaxies. Sum: %f %f'%(np.sum(rnd_H)/100.0,np.sum(rnd_H)/100.0*2.5))

plt.figure()
plt.pcolor(X,Y,red_H.T,cmap=plt.cm.PuBu,norm=LogNorm(vmin=val_min,vmax=val_max))
plt.colorbar()
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('SDSS RedMapper Cluster Galaxies. Sum: %f'%(np.sum(red_H)/100.0))

plt.figure()
plt.pcolor(X,Y,red_H.T-rnd_H.T,cmap=plt.cm.PuBu,norm=LogNorm(vmin=val_min,vmax=val_max))
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('RedMapper Cluster Gal - Random Gal')
plt.colorbar()

plt.figure()
plt.pcolor(X,Y,rnd_H.T-red_H.T,cmap=plt.cm.PuBu,norm=LogNorm(vmin=val_min,vmax=val_max))
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('RedMapper Random Gal - Cluster Gal')
plt.colorbar()

plt.figure()
plt.pcolor(X,Y,red_H.T/np.sum(red_H)-rnd_H.T/np.sum(rnd_H),cmap=plt.cm.PuBu,norm=LogNorm())
plt.colorbar()
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('normed RedMapper Cluster Gal - normed Random Gal')

plt.figure()
plt.pcolor(X,Y,rnd_H.T/np.sum(rnd_H)-red_H.T/np.sum(red_H),cmap=plt.cm.PuBu,norm=LogNorm())
plt.colorbar()
plt.ylabel('g-r')
plt.xlabel('i')
plt.title('normed Random Gal - normed Cluster Gal')

plt.figure()
plt.title('Given Data. Sum: %f'%(np.sum(expected)/100.0))
plt.pcolor(X,Y,expected.T,cmap=plt.cm.PuBu,norm=LogNorm(vmin=val_min,vmax=val_max))
plt.colorbar()
plt.ylabel('g-r')
plt.xlabel('i')

plt.show()

