import matplotlib.pyplot as plt
import numpy as np
from cluster_rotation_analysis import clust_rot_utils as cru
from cluster_rotation_analysis import vector_point_diagram_general as vpdg
import pdb
from astropy.table import Table
from cluster_rotation_analysis import scatter_proper_motion as spm

catalog_filename = '/home/rlewismerrill/cluster_rotation/' +\
                   'quint_catalog_probs_4_Aks_H18b_Pcol.fits'
clust_mem_min = 0.0

# Read in the catalog
catalog = cru.load_cluster_members(fits_file_path = catalog_filename,
                                   clust_mem_min = clust_mem_min)

# Column names
x_col="x_2010_f153m"
y_col="y_2010_f153m"
xerr_col="xe_2010_f153m"
yerr_col="ye_2010_f153m"
pm_x_col="fit_vx"
pm_y_col="fit_vy"
pm_xerr_col="fit_vxe"
pm_yerr_col="fit_vye"

# Cluster center position in pixel units 
radec_center=335
dec_center=499

# Observations pixel scale (mas/pixel)
plate_scale=121.625

# Proper motions of stars (in mas/yr; catalog in pixel units)
x_pm = catalog[pm_x_col].data * plate_scale
y_pm = catalog[pm_y_col].data * plate_scale
xerr_pm = catalog[pm_xerr_col].data * plate_scale
yerr_pm = catalog[pm_yerr_col].data * plate_scale

transformed_catalog = spm.transform_pixel_to_sky(\
    catalog_or_file = catalog_filename,
    x_col="x_2010_f153m",
    y_col="y_2010_f153m",
    x_err_col="xe_2010_f153m",
    y_err_col="ye_2010_f153m",
    pm_x_col="fit_vx",
    pm_y_col="fit_vy",
    pm_x_err_col="fit_vxe",
    pm_y_err_col="fit_vye",
    pa_deg=-42.3,
    pixel_scale=121.625,
    scale_divisor=1000.0,
    pixel_east_left=False,
    clust_mem_min = clust_mem_min)

transformed_catalog_with_gal = spm.transform_pixel_to_sky(\
    catalog_or_file = catalog_filename,
    x_col="x_2010_f153m",
    y_col="y_2010_f153m",
    x_err_col="xe_2010_f153m",
    y_err_col="ye_2010_f153m",
    pm_x_col="fit_vx",
    pm_y_col="fit_vy",
    pm_x_err_col="fit_vxe",
    pm_y_err_col="fit_vye",
    pa_deg=-42.3,
    pixel_scale=121.625,
    scale_divisor=1000.0,
    pixel_east_left=False,
    clust_mem_min = clust_mem_min,
    to_galactic = True,
    ra_zero_deg = 266.56304167, # Q12 from Rui+2019
    dec_zero_deg = -28.82630556) #Q12 from Rui+2019

# Labels

xlabel = 'X Proper Motion [mas/yr]'
ylabel = 'Y Proper Motion [mas/yr]'
fontsize = 20
limits = (-10,10)

fig1, ax_scatter1, ax_histx1, ax_histy1 = vpdg.vector_point_diagram(\
                                      xdata = x_pm,
                                      ydata = y_pm,
                                      xerr = np.zeros(len(x_pm)),
                                      yerr = np.zeros(len(y_pm)),
                                      #xerr = xerr_pm,
                                      #yerr = yerr_pm,
                                      label_fontsize = fontsize,
                                      decimals_round = 2,
                                      xlabel = xlabel,
                                      ylabel = ylabel,
                                      xlim = limits,
                                      ylim = limits,
                                    normed_gaussians = True,
                                    gaussian_label = 'Norm. #', 
                                    errorbar_kwargs = {'fmt': 'o',
                                                        'color': 'k',
                                                        'capsize': 2,
                                                        'elinewidth': 2,
                                                        'alpha': 0.5},
                                    xhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    yhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    fig_kwargs = {'figsize': (12.5,10)})

fig2, ax_scatter2, ax_histx2, ax_histy2 = vpdg.vector_point_diagram(\
                                      xdata = transformed_catalog['pm_radec_sky'],
                                      ydata = transformed_catalog['pm_dec_sky'],
                                      xerr = np.zeros(len(transformed_catalog['pm_radec_sky'])),
                                      yerr = np.zeros(len(transformed_catalog['pm_dec_sky'])),
                                      #xerr = transformed_catalog['pm_radec_sky_err'],
                                      #yerr = transformed_catalog['pm_dec_sky_err'],
                                      label_fontsize = fontsize,
                                      decimals_round = 2,
                                      xlabel = 'RAcosDEC PM [mas/yr]',
                                      ylabel = 'DEC PM [mas/yr]',
                                      xlim = limits,
                                      ylim = limits,
                                    normed_gaussians = True,
                                    gaussian_label = 'Norm. #', 
                                    errorbar_kwargs = {'fmt': 'o',
                                                        'color': 'k',
                                                        'capsize': 2,
                                                        'elinewidth': 2,
                                                        'alpha': 0.5},
                                    xhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    yhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    fig_kwargs = {'figsize': (12.5,10)})

fig3, ax_scatter3, ax_histx3, ax_histy3 = vpdg.vector_point_diagram(\
                                      xdata = transformed_catalog_with_gal['pm_l_gal'],
                                      ydata = transformed_catalog_with_gal['pm_b_gal'],
                                      xerr = np.zeros(len(transformed_catalog_with_gal['pm_l_gal'])),
                                      yerr = np.zeros(len(transformed_catalog_with_gal['pm_b_gal'])),
                                      label_fontsize = fontsize,
                                      decimals_round = 2,
                                      xlabel = 'l cos(b) PM [mas/yr]',
                                      ylabel = 'b PM [mas/yr]',
                                      xlim = limits,
                                      ylim = limits,
                                    normed_gaussians = True,
                                    gaussian_label = 'Norm. #', 
                                    errorbar_kwargs = {'fmt': 'o',
                                                        'color': 'k',
                                                        'capsize': 2,
                                                        'elinewidth': 2,
                                                        'alpha': 0.5},
                                    xhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    yhist_kwargs = {'histtype': 'step',
                                                    'linewidth': 3,
                                                    'color': 'k',
                                                    'bins': 11},
                                    fig_kwargs = {'figsize': (12.5,10)})

fig1.savefig('vpd_QuintIMF_Pmem0_obs_orient.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('vpd_QuintIMF_Pmem0_sky_orient.png', dpi = 300, bbox_inches = 'tight')
fig3.savefig('vpd_QuintIMF_Pmem0_gal_orient.png', dpi = 300, bbox_inches = 'tight')


