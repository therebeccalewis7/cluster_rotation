import matplotlib.pyplot as plt
import numpy as np
from cluster_rotation_analysis import clust_rot_utils as cru
import matplotlib.colors as colors
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pdb

#Temp; Remove after testing

#xdata_random = np.random.normal(loc = 0, scale = 5, size = 1000)
#ydata_random = np.random.normal(loc = 1, scale = 5, size = 1000)
#xdata_err_random = np.random.uniform(size = 1000)
#ydata_err_random = np.random.uniform(size = 1000)

def vector_point_diagram(xdata = None, 
                         ydata = None, 
                         xerr = None, 
                         yerr = None,
                         label_fontsize = 26, 
                         decimals_round = None,
                         xlim = None,
                         ylim = None,
                         xlabel = 'X',
                         ylabel = 'Y',
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
                         fig_kwargs = {'figsize': (12.5,10)}):

    """
    Plots a vector point diagram of the given x and y data
    along with their distributions in x and y in histograms along 
    the sides.
    """

    # Initialize figure and set up the axes
    fig = plt.figure(**fig_kwargs)

    # Sizes of the side plots
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # Add scatter Axes and the histogram Axes
    ax_scatter = fig.add_axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    
    ax_histx = fig.add_axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    
    ax_histy = fig.add_axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # If limits for the x and y axis are given
    if type(xlim) != type(None):
        ax_scatter.set_xlim(xlim)
        ax_histx.set_xlim(xlim)
    if type(ylim) != type(None):
        ax_scatter.set_ylim(ylim)
        ax_histy.set_ylim(ylim)

    ax_scatter.grid(True)

    # Plot scatter plot data
    ax_scatter.errorbar(x = xdata,
                        y = ydata,
                        xerr = xerr,
                        yerr = yerr,
                        **errorbar_kwargs)

    # Plot histogram data
    ax_histx.hist(x = xdata,
                  orientation = 'vertical',
                  density = normed_gaussians,
                  range = xlim,
                  **xhist_kwargs)

    ax_histy.hist(x = ydata,
                  orientation = 'horizontal',
                  density = normed_gaussians,
                  range = ylim,
                  **yhist_kwargs)

    # Determine tick marks for scatter plot x-axis
    xtick_locs, xformatted_labels = cru.auto_space_ticks(\
                                                        ax = ax_scatter, 
                                                        axis='x', 
                                                        padding=0.1, 
                                                        fontsize=label_fontsize, 
                                                        fig=fig,
                                                        decimals = decimals_round)

    # Determine tick marks for scatter plot y-axis
    ytick_locs, yformatted_labels = cru.auto_space_ticks(\
                                                        ax = ax_scatter, 
                                                        axis='y', 
                                                        padding=0.1, 
                                                        fontsize=label_fontsize, 
                                                        fig=fig,
                                                        decimals = decimals_round)

    ax_scatter.set_xticks(xtick_locs, labels=xformatted_labels, fontsize=label_fontsize)
    ax_scatter.set_yticks(ytick_locs, labels=yformatted_labels, fontsize=label_fontsize)

    # Determine tick marks for x-axis histogram
    xhist_tick_locs, xhist_tick_labels = cru.auto_space_ticks(\
                                                        ax = ax_histx, 
                                                        axis='y', 
                                                        padding=0.25, 
                                                        fontsize=label_fontsize, 
                                                        fig=fig,
                                                        decimals = decimals_round)

    # Determine tick marks for x-axis histogram
    yhist_tick_locs, yhist_tick_labels= cru.auto_space_ticks(\
                                                        ax = ax_histy, 
                                                        axis='x', 
                                                        padding=0.25, 
                                                        fontsize=label_fontsize, 
                                                        fig=fig,
                                                        decimals = decimals_round)

    ax_histx.set_xticks(xtick_locs, labels=xformatted_labels, fontsize=label_fontsize)
    ax_histy.set_yticks(ytick_locs, labels=yformatted_labels, fontsize=label_fontsize)
    
    ax_histx.set_yticks(xhist_tick_locs, labels=xhist_tick_labels, fontsize=label_fontsize)
    ax_histy.set_xticks(yhist_tick_locs, labels=yhist_tick_labels, fontsize=label_fontsize)

    ax_histx.set_ylabel(gaussian_label, fontsize = label_fontsize)
    ax_histy.set_xlabel(gaussian_label, fontsize = label_fontsize)
    ax_scatter.set_xlabel(xlabel, fontsize = label_fontsize)
    ax_scatter.set_ylabel(ylabel, fontsize = label_fontsize)

    return fig, ax_scatter, ax_histx, ax_histy

    