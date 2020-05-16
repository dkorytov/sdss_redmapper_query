#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import h5py
import dtk

import background

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

def pyplot_zoom_out(all=0.0, horizontal=None, vertical=None, right=None, left=None, top=None, bottom=None):
    """ increases the axes by fractional amount in each direction"""
    # Selects the first non-None value from the 'or' list
    right = right or horizontal or all
    left = left or horizontal or all
    top = top or vertical or all
    bottom = bottom or vertical or all
    
    ylim = plt.ylim()
    ydiff_cen = (ylim[1]-ylim[0])/2
    ycen = np.average(ylim)
    new_ylim =[ycen-ydiff_cen*(1.0+bottom), ycen+ydiff_cen*(1.0+top)]
    # plt.ylim(
    xlim = plt.xlim()
    xdiff_cen = (xlim[1]-xlim[0])/2
    xcen = np.average(xlim)
    new_xlim = [xcen-xdiff_cen*(1.0+left), xcen+xdiff_cen*(1.0+right)]
    plt.plot(new_xlim, new_ylim, ',', alpha=0.0)
    plt.plot(new_xlim, new_ylim[::-1], ',', alpha=0.0)
    print(plt.xlim())

def plot_saved_clusters(param_fname):
    param = dtk.Param(param_fname)
    query_data_file = param.get_string('query_data_folder')+'query_results.hdf5'
    mstar_offset = param.get_int('mstar_cut_offset')

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
        
    if 'background_r200_factor' in param:
        background_r200_factor = param.get_float('background_r200_factor')
    else:
        background_r200_factor = 1.0
        
    if 'richness_mass_author' in param:
        richness_mass_author = param.get_string("richness_mass_author")
    else:
        richness_mass_author = None
        
    if 'convert_m200m_to_m200c' in param:
        convert_m200m_to_m200c = param.get_bool('convert_m200m_to_m200c')
    else:
        convert_m200m_to_m200c = False
    if 'convert_m200c_to_m200m' in param:
        convert_m200c_to_m200m = param.get_bool('convert_m200c_to_m200m')
    else:
        convert_m200c_to_m200m = False
    #######################
    ## Processing Params ##
    #######################

    background_z_bins = np.linspace(background_z_start, background_z_end, int(background_z_num))

    if(galaxy_type == 1):
        galaxy_type_name = "all"
    elif(galaxy_type == 2):
        galaxy_type_name = "red"
    elif(galaxy_type == 3):
        galaxy_type_name = "non-red"
    else:
        print("unknwon galaxy type")
        raise 
    if(galaxy_weight == 1):
        galaxy_weight_name = "normal"
    elif(galaxy_weight == 2):
        galaxy_weight_name = "red squence"
    elif(galaxy_weight == 3):
        galaxy_weight_name = 'red squence rykoff'
    else:
        print("unknown galaxy weight")
        raise

    if(galaxy_color_type ==1):
        galaxy_color_name = 'g-r'
    elif(galaxy_color_type==2):
        galaxy_color_name = 'g-r - RS(z)'
    if(galaxy_mag_type ==1):
        galaxy_mag_name = 'i'
    elif(galaxy_mag_type==2):
        galaxy_mag_name = 'i-m*(z)'
        
    background.set_gal_clr_type(galaxy_color_type,galaxy_color_lim)
    background.set_gal_mag_type(galaxy_mag_type,  galaxy_mag_lim)
    background.set_gal_type_weight(galaxy_type,galaxy_weight)
    background.set_cosmology(cosmology_name)
    background.set_mstar_cut(mstar_cut_offset)
    
    background.set_mstar_cut(param.get_float('mstar_cut_offset'))
    [dataclstr,datagal,data_pass_mask,clstr_num] =  background.get_clstr(query_type,
                                                                         query_results,
                                                                         100,
                                                                         till_end=False,

                                                                         # query_cluster_num,
                                                                         # till_end=query_cluster_all,
                                                                         richness_mass_author=richness_mass_author,
                                                                         convert_m200c_to_m200m=convert_m200c_to_m200m,
                                                                         convert_m200m_to_m200c=convert_m200m_to_m200c)


    for i in range(0, 400):
        print(list(dataclstr[i].keys()))
        print(list(datagal[i].keys()))
        gal_prop=dataclstr[i]
        gal=datagal[i]
        rad = float(gal_prop['r200_arcmin'])/60
        dec = float(gal_prop['dec'])
        ra = float(gal_prop['ra'])
        mass = float(gal_prop['m200'])
        z = float(gal_prop['z'])
        if dec > 10 or mass <2e14 or z>0.35:
            continue
        plt.figure()
        print(ra, dec, rad)
        # c1 = plt.Circle((ra, dec), rad, fill=True, color='r', lw=2)
        e_r200 = matplotlib.patches.Ellipse((ra,dec), rad/np.cos(dec/360)*2.0, rad*2.0, fill=False, color='k', lw=2)
        plt.gca().add_artist(e_r200)
        # plotting radial bins
        for j in range(0, 16):
            rad_bin = float(j)/16*rad
            e_bin = matplotlib.patches.Ellipse((ra,dec),
                                               rad_bin/np.cos(dec/360)*2.0,
                                               rad_bin*2.0, fill=False,
                                               color='k', lw=1,
                                               ls='--')
            plt.gca().add_artist(e_bin)
        
        plt.ylabel('Dec [deg]')
        plt.xlabel('Ra [deg]')
        plt.plot([],[],'k', lw=2, label='$R_{200c}$')
        plt.plot([],[], 'k', ls='--', label='radial bins')
        
        ## Galaxies
        #          bright               saturated             satur_center         nopetro             deblended_as_moving
        cut_flgs = 0x0000000000000002 | 0x0000000000040000  | 0x0000080000000000 | 0x0000000000000100 | 0x0000000100000000
        flgs = gal['flags_i'] | gal['flags_r'] | gal['flags_g']
        slct1 = gal['type']!=31
        slct2 = np.logical_not(flgs&cut_flgs)
        slct3 = gal['m_i']<background.m_cut(z)
        slct = slct1 & slct2 & slct3
        
        plt.plot(gal['ra'][slct], gal['dec'][slct], '.', alpha=0.6, label='galaxies')
        plt.text(0.05, 0.05, 'M$_{{200c}}$=\n{:2.2e}[$h^{{-1}}$M$_{{\odot}}$]'.format(mass),
                 transform=plt.gca().transAxes, verticalalignment='bottom',
                horizontalalignment='left')
        plt.text(0.95, 0.05, 'z={:.2f}'.format(z),
                 transform=plt.gca().transAxes, verticalalignment='bottom',
                 horizontalalignment='right')
        
        plt.legend(framealpha=0.0, labelspacing=0)
        plt.axis('equal')
        pyplot_zoom_out(top=0.2, right=-0.5, left=0.3)
        plt.tight_layout()
        dtk.save_figs('figs/'+__file__+'/','.pdf')
        plt.close('all')

if __name__ == "__main__":
    param_fname = sys.argv[1]
    plot_saved_clusters(param_fname)
    plt.show()
