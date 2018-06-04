#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt

import h5py 

import sys
import dtk

print "Testing if there is any bias in the mask cut out"
#I would expect that exclusion rates will scale roughly linearly 
#with area. 

param = dtk.Param(sys.argv[1])
mask_hfilename = param.get_string('query_data_folder')+'gal_mask.hdf5'
gal_hfilename  =param.get_string('query_data_folder')+'query_results.hdf5'
mask_hfile = h5py.File(mask_hfilename,'r')
gal_hfile  = h5py.File(gal_hfilename,'r')

z_bins = np.linspace(param.get_float('z_bins_start'),
                     param.get_float('z_bins_end'),
                     20)
z_bins_avg = (z_bins[:-1]+z_bins[1:])/2.0

mass_bins = np.logspace(param.get_float('mass_bins_start'),
                        param.get_float('mass_bins_end'),
                        20)

rad_bins = np.logspace(-3,2,100)
rad_bins_avg = (rad_bins[:-1]+rad_bins[1:])/2.0

mass_bins_avg = (mass_bins[:-1]+mass_bins[1:])/2.0
print z_bins_avg
print mass_bins_avg

z_pass = np.zeros_like(z_bins)
z_fail = np.zeros_like(z_bins)

mass_pass = np.zeros_like(mass_bins)
mass_fail = np.zeros_like(mass_bins)

z_mass_pass = np.zeros((z_bins.size,mass_bins.size))
z_mass_fail = np.zeros((z_bins.size,mass_bins.size))

rad_pass = np.zeros_like(rad_bins)
rad_fail = np.zeros_like(rad_bins)

for i in range(0,2000):
    mask_pass = mask_hfile['%d'%i]['mask_pass'][:]
    mass = gal_hfile['gal_prop%d'%i]['mass'] #m200 mass
    rad = gal_hfile['gal_prop%d'%i]['rad'] #rad arcmin radius on sky
    z    = gal_hfile['gal_prop%d'%i]['z']
    z_i = np.searchsorted(z_bins,z)-1
    m_i = np.searchsorted(mass_bins,mass)-1
    r_i = np.searchsorted(rad_bins,rad)-1
    if(mask_pass):
        z_pass[z_i]+=1
        mass_pass[m_i]+=1
        z_mass_pass[z_i,m_i]+=1
        rad_pass[r_i]+=1
    else:
        z_fail[z_i]+=1
        mass_fail[m_i]+=1
        z_mass_fail[z_i,m_i]+=1
        rad_fail[r_i]+=1

z_pass_rate = z_pass/(z_pass+z_fail)
mass_pass_rate = mass_pass/(mass_pass + mass_fail)
z_mass_pass_rate = z_mass_pass/(z_mass_pass+z_mass_fail)
rad_pass_rate = rad_pass/(rad_pass+rad_fail)

plt.figure()
plt.plot(z_bins_avg,z_pass_rate[:-1],'x-')
err = np.sqrt(1.0/z_pass[:-1]+1.0/(z_pass[:-1]+z_fail[:-1]))
err = (1.0/np.sqrt(z_pass[:-1]))*z_pass[:-1]/(z_pass[:-1]+z_fail[:-1])
plt.fill_between(z_bins_avg,z_pass_rate[:-1]-err,z_pass_rate[:-1]+err,edgecolor=None,alpha=0.3)
plt.ylim([0,1.1])
plt.xlabel('red shift bin')
plt.ylabel('pass rate')
plt.grid()

plt.figure()
plt.plot(mass_bins_avg,mass_pass_rate[:-1],'x-')
err = 1.0/np.sqrt(mass_pass[:-1]+mass_fail[:-1])
err = np.sqrt(1.0/mass_pass[:-1]+1.0/(mass_pass[:-1]+mass_fail[:-1]))
err = (1.0/np.sqrt(mass_pass[:-1]))*mass_pass[:-1]/(mass_pass[:-1]+mass_fail[:-1])
plt.fill_between(mass_bins_avg,mass_pass_rate[:-1]-err,mass_pass_rate[:-1]+err,edgecolor=None,alpha=0.3)
plt.ylim([0,1.1])
plt.xscale('log')
plt.xlabel('cluster m200c [Msun/h]')
plt.ylabel('pass rate')
plt.grid()

plt.figure()
plt.plot(np.pi*rad_bins_avg**2,rad_pass_rate[:-1],'x-')
err = 1.0/np.sqrt(rad_pass[:-1]+rad_fail[:-1])
err = np.sqrt(1.0/rad_pass[:-1]+1.0/(rad_pass[:-1]+rad_fail[:-1]))
err = (1.0/np.sqrt(rad_pass[:-1]))*rad_pass[:-1]/(rad_pass[:-1]+rad_fail[:-1])
plt.fill_between(np.pi*rad_bins_avg**2,rad_pass_rate[:-1]-err,rad_pass_rate[:-1]+err,edgecolor=None,alpha=0.3)
plt.xlabel('cluster r200 area [arcmin^2]')
plt.ylabel('pass rate')
plt.ylim([0,1.1])
plt.grid()

import numpy.ma as ma
z_mass_pass_rate= ma.array(z_mass_pass_rate,mask=np.isnan(z_mass_pass_rate))
plt.figure()
mycmap = plt.cm.gray

mycmap.set_bad('lightblue',2)

plt.pcolormesh(z_bins,mass_bins,z_mass_pass_rate[:-1,:-1].T,cmap=mycmap,vmin=0.0,vmax=1.0)
plt.colorbar().set_label('pass rate')
plt.yscale('log')
plt.xlabel('z')
plt.ylabel('cluster m200c [Msun/h]')
plt.grid()

plt.figure()
import matplotlib.colors as colors
plt.pcolormesh(z_bins,mass_bins,z_mass_pass.T+z_mass_fail.T,cmap='PuBu',norm=colors.LogNorm())
plt.colorbar().set_label('Number of Clusters')
plt.yscale('log')
plt.xlabel('z')
plt.ylabel('cluster m200c [Msun/h]')
plt.grid()

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')
print z_mass_pass_rate.shape
Z,M = np.meshgrid(z_bins,np.log10(mass_bins))
ax.plot_surface(Z,M,np.log10(z_mass_pass+z_mass_fail+1).T,rstride=1,cstride=1,linewidth=0,facecolors=cm.gray(np.nan_to_num(z_mass_pass_rate.T)))
ax.set_ylabel('Log10 mass [Msun/h]')
ax.set_xlabel('z')
ax.set_zlabel('Log10 Freq')
dtk.save_figs("figs/"+param.file+"/mask_bias_check/")

plt.show()
