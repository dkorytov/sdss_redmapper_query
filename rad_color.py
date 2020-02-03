#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from background import *
import sys

param = dtk.Param(sys.argv[1])

fname = param.get_string("result_folder")#"data/normal/result/"
gal_type = 1
gal_weight=1

data = load_radial_profile(fname,gal_type,gal_weight)

z_bins      = data['z_bins']
mass_bins   = data['mass_bins']
rad_bins    = data['rad_bins']
zm_cnt      = data['zm_cnt']
zmr_cnt     = data['zmr_cnt']
zmr_clr_mag = data['zmr_clr_mg']
clr_bins    = data['clr_bins']
mg_bins     = data['mg_bins']    

z_num = len(z_bins)-1
mass_num = len(mass_bins)-1
rad_num  = len(rad_bins)-1
z_bin_avg = (rad_bins[1:]+rad_bins[0:-1])/2.0
clr_bin_avg = (clr_bins[1:]+clr_bins[0:-1])/2.0
rad_bin_avg = (rad_bins[1:]+rad_bins[0:-1])/2.0
mass_bin_avg = (mass_bins[1:]+mass_bins[0:-1])/2.0
zmr_avg_color = np.zeros_like(zmr_cnt,dtype='f4')

plt.figure()
z1 = 0
z2 = 3
mg_avg,clr_avg = get_color_mag_bin_avgs()
slct = red_squence_redder_select((z_bins[z1]+z_bins[z1+1])/2.0,z_bins[z2],zmr_clr_mag[3,3,0])
rs_clr1 = redsequence_line(mg_avg,z_bins[z1])
rs_clr2 = redsequence_line(mg_avg,z_bins[z2])
plt.plot(mg_avg,rs_clr1,'rx-')
plt.plot(mg_avg,rs_clr2,'rx-')
print np.sum(slct)
plt.title("red squence selection\nz~%.2f,m~%.2e"%(z_bin_avg[3],mass_bin_avg[3]))
plt.pcolor(mg_bins,clr_bin_avg,slct.T,cmap=plt.cm.BuPu)
plt.colorbar()
plt.tight_layout()
for z_i in range(0,z_num):
    plt.figure()
    plt.title('%.2f < z < %.2f'%(z_bins[z_i],z_bins[z_i+1]))
    for m_i in range(0,mass_num):
        red_fract = np.zeros(rad_num)
        for r_i in range(0,rad_num):
            hclrmg = zmr_clr_mag[z_i,m_i,r_i].clip(min=0.0)
            red_gal = np.sum(hclrmg*red_squence_redder_select((z_bins[z_i]+z_bins[z_i+1])/2.0,z_bins[z_i+1],hclrmg))
            tot_gal = np.sum(hclrmg)
            red_fract[r_i]=red_gal/tot_gal
        if(zm_cnt[z_i,m_i] != 0):
            plt.plot(rad_bin_avg,red_fract,label='m~%1.2e (n=%d)'%(mass_bin_avg[m_i],zm_cnt[z_i,m_i]))
    plt.legend(loc='best')
    plt.ylabel('Red Fraction')
    plt.xlabel('radius [R200]')
    plt.grid()
    plt.legend()
    plt.tight_layout()


print zm_cnt.shape
for z_i in range(0,z_num):
    for m_i in range(0,mass_num):
        for r_i in range(0,rad_num):
            data_flattened_clr = np.sum(zmr_clr_mag[z_i,m_i,r_i],axis=0)
            slct = data_flattened_clr>0
            if(np.sum(data_flattened_clr[slct]) != 0):
               average_color = np.average(clr_bin_avg[slct],weights=data_flattened_clr[slct])
               if(np.abs(average_color) > 5):
                   plt.figure()
                   plt.title('avg=%.2f, z=%.1f\nm=%1.1e,rad=%.2f'%(average_color,z_bin_avg[z_i],mass_bin_avg[m_i],rad_bin_avg[r_i]))
                   mmax = np.max(np.abs(zmr_clr_mag[z_i,m_i,r_i]))
                   plt.pcolor(mg_bins,clr_bins,zmr_clr_mag[z_i,m_i,r_i].T,cmap=plt.cm.RdBu,vmax=mmax,vmin=-mmax)
                   plt.colorbar()
                   plt.tight_layout()
               zmr_avg_color[z_i,m_i,r_i] = average_color        
            else:
                zmr_avg_color[z_i,m_i,r_i]=np.nan
                
print zmr_avg_color.shape

for z_i in range(0,z_num):
    plt.figure()
    plt.title("%.2f<z<%.2f"%(z_bins[z_i],z_bins[z_i+1]))
    for m_i in range(0,mass_num):
        tmp =zm_cnt[z_i,m_i]
        if(tmp>0):
            plt.plot(rad_bin_avg,zmr_avg_color[z_i,m_i],"x-",label="m~%.2e (n=%d)"%(np.log10(mass_bin_avg[m_i]),tmp))
    plt.legend(loc='best')
    plt.grid()
    plt.ylim([0.0,3.0])
    plt.tight_layout()
    print "--"
print clr_bin_avg

m_i = 2
z_i = 2
cnt = 0
cnt_limit = 4
data = np.zeros_like(zmr_clr_mag[3,3,0])
for i in range(0,len(rad_bins)-1):
    print i,cnt
    if(cnt<cnt_limit):
        data +=zmr_clr_mag[3,3,i]
        cnt+=1
    if(cnt==cnt_limit or i == len(rad_bins)-2):
        print "*", len(rad_bins)-2
        plt.figure()
        plt.title('%.3f < r < %.3f [r200]\nz~%.2f, m~%1.2e'%(rad_bins[i-cnt+1],rad_bins[i+1],z_bin_avg[z_i],mass_bin_avg[m_i]))
        plt.pcolor(mg_bins,clr_bins,data.T/np.float(cnt),cmap=plt.cm.BuPu,norm=LogNorm())
        plt.colorbar()
        plt.grid()
        data_flattened_clr = np.sum(data,axis=0)
        slct= data_flattened_clr>0.0
        print data_flattened_clr.shape, clr_bin_avg.shape
        average_color = np.average(clr_bin_avg[slct],weights=data_flattened_clr[slct])
        print average_color
        plt.axhline(average_color,c='r')
        plt.text(12.05,average_color+0.1,"avg clr=%1.2f"%(average_color),color='red')
        plt.tight_layout()
        data=np.zeros_like(zmr_clr_mag[3,3,0])
        cnt = 0

plt.figure()
slct = red_squence_select(z_bins[2],z_bins[3],zmr_clr_mag[3,3,0])
plt.title("red squence")
plt.pcolor(mg_bins,clr_bin_avg,slct.T,cmap=plt.cm.BuPu,norm=LogNorm())
plt.tight_layout()



dtk.save_figs("figs/"+param.file+"/rad_color/")
plt.show()

