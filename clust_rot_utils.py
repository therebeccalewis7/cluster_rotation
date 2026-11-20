from astropy.table import Table
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import pdb
import copy


def load_cluster_members(fits_file_path,
                         clust_mem_min = 0.3,
                         clust_mem_col1 = 'Membership',
                         clust_mem_col2 = 'Membership_color',
                         verbose = True):
    """
    """
    # Load in the fits file containing the catalog
    tab = Table.read(fits_file_path)

    # Determine which stars are cluster members by the given cut
    try:
        if (type(clust_mem_col1) != type(None)) & (type(clust_mem_col2) != type(None)):
                clust_membership = tab[clust_mem_col1].data * tab[clust_mem_col2].data
        elif (type(clust_mem_col1) != type(None)) & (type(clust_mem_col2) == type(None))
                clust_membership = tab[clust_mem_col1].data

    except KeyError:
        clust_membership = np.ones(len(tab))

        if verbose:
            print('No matching membership probability keyword in fits file. No cuts applied.')

    # Trim the catalog
    tab_bool = clust_membership >= clust_mem_min
    tab_trimmed = tab[tab_bool]

    return tab_trimmed

def get_axis_font_sidelength_fraction(ax, fontsize, 
                                      num_digits = 1,
                                      has_decimal = False,
                                      is_negative = False,
                                      axis='y', 
                                      fig = None):
    """
    Query an active Matplotlib Axis object to find the fraction of its 
    width or height covered by a specific font size.

    Arguments:
    ----------
    ax: matplotlib.axes.Axes
        Matplotlib Axes object

    fontsize: int
        Fontsize of the label

    num_digits: int
        The number of digits in the tick mark

    has_decimal: bool
        Whether or not the tick mark has a decimal

    axis: str; 'y' or 'x'
        Which axis we are measuring along

    fig: matplotlib Figure
        Matplotlib Figure object that the axis is 
        presumably on
    """
    if fig is None:
        fig = ax.get_figure()

    # Relative figure coordinates [x0, y0, width, height]
    bbox = ax.get_position() 

    #bbox = ax.get_tightbbox()

    # Allowed forms of axis string values
    allowed_axis_values = ['x', 'y', 'X', 'Y']

    if any(s in axis for s in allowed_axis_values):
        None
    else:
        raise Exception("Axis string value not supported.")

    if (axis == 'x') or (axis == 'X'):
        fig_size_inches = fig.get_size_inches()[0] # Figure width in inches
        axis_side_inches = fig_size_inches * bbox.width

        measured_str = []
        for ii in np.arange(num_digits):
            measured_str.append('0')
        if has_decimal == True:
            measured_str.append('.')
        if is_negative == True:
            measured_str.append('-')

        measured_str = ''.join(measured_str)

        str_length_in = get_str_width_inches(fontsize = fontsize,
                                                 measured_str = measured_str,
                                                 fig = fig)

        return str_length_in/axis_side_inches

    elif (axis == 'y') or (axis == 'Y'):
        fig_size_inches = fig.get_size_inches()[1] # Figure height in inches
        axis_side_inches = fig_size_inches * bbox.height
        
        font_height_inches = fontsize / 72.0
    
        return font_height_inches / axis_side_inches

    else:

        return None


def get_str_width_inches(fontsize, measured_str="0", fig=None):
    """
    Returns the physical width of a string in inches 
    for a given fontsize and active font family.

    Arguments:
    -------------
    fontsize: int
        Fontsize being used on the figure

    measured_str: str
        The string whose width in inches will be measured

    fig: matplotlib Figure
        Matplotlib Figure that this fontsize is on
    """
    if fig is None:
        fig = plt.gcf()
        
    renderer = fig.canvas.get_renderer()
    font_props = FontProperties(size=fontsize)
    
    # Measure the text bounding box using Matplotlib's internal renderer
    text_obj = plt.text(0, 0, str(measured_str), fontproperties=font_props)
    bbox = text_obj.get_window_extent(renderer=renderer)
    text_obj.remove()  # Clean up temporary text object
    
    # Convert pixels to inches (dpi = dots per inch)
    width_inches = bbox.width / fig.dpi

    return width_inches

def auto_space_ticks(ax, axis='x', 
                     padding=0.1, 
                     fontsize=26, 
                     fig=None,
                     decimals = None,
                     verbose = True):
    """
    Calculates and applies the optimal number of ticks for an axis based on
    string size dimensions, padding, and scientific notation compatibility.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        Target Axes object.
    axis : str
        'x' or 'y' axis. 
        Allowed values for the x-axis are 'x' or 'X'. 
        Allowed values for the y-axis are 'y' or 'Y'.
    padding : float
        Fraction of the total axis sidelength to use as spacing between ticks.
        Goes between 0.0 and 1.0.
    fontsize : int
        Font size of the tick labels.
    fig : matplotlib.figure.Figure, optional
        Target Figure object. Inferred from `ax` if None.
    decimals: int
        Number of decimal places to round to on ticks
    verbose: bool
        Optional bool for troubleshooting

    Returns:
    --------
    locs : numpy.ndarray
        The generated tick locations.
    labels : list of str
        The formatted tick string labels (including scientific notation formatting if active).
    """
    if fig is None:
        fig = ax.get_figure()

    # Copy the figure and axis so that the original objects 
    # remain unaffected

    ax_copy = copy.deepcopy(ax)
    fig_copy = ax_copy.get_figure()
    #fig_copy = copy.deepcopy(fig)
    #ax_copy = fig_copy.axes[0]

    # Allowed forms of axis string values
    allowed_axis_values = ['x', 'y', 'X', 'Y']

    if any(s in axis for s in allowed_axis_values):
        None
    else:
        raise Exception("Axis string value not supported.")

    # Retrieve view limits range for the axis
    if (axis == 'x' or axis == 'X'):
        axis_obj = ax_copy.xaxis
    elif (axis == 'y' or axis == 'Y') :
        axis_obj = ax_copy.yaxis

    vmin, vmax = axis_obj.get_view_interval()

    # Create a temporary formatter to generate representative tick strings 
    # (handles scientific notation)
    formatter = axis_obj.get_major_formatter()
    if not isinstance(formatter, ScalarFormatter):
        formatter = ScalarFormatter()
        formatter.set_axis(axis_obj)
    
    # Estimate extreme tick label representations across range to measure max dimensions
    if type(decimals) == type(None):
        sample_values = np.linspace(vmin, vmax, 5)
    else:
        sample_values = np.round(np.linspace(vmin, vmax, 5), decimals = decimals)

    formatter.set_locs(sample_values)
    sample_strings = [formatter(val) for val in sample_values]

    # Determine maximum number of digits and decimal status among candidate tick labels
    max_num_digits = 0
    has_decimal = False
    is_negative = False
    
    for s in sample_strings:
        digits_only = [c for c in s if c.isdigit()]
        max_num_digits = max(max_num_digits, len(digits_only))
        if '.' in s:
            has_decimal = True
        if any(np.array(digits_only, dtype = float) < 0.0):
            is_negative = True

    # Get the fractional coverage of a single label along the target axis
    label_fraction = get_axis_font_sidelength_fraction(
        ax=ax_copy,
        fontsize=fontsize,
        num_digits=max_num_digits,
        has_decimal=has_decimal,
        is_negative = is_negative,
        axis=axis,
        fig=fig_copy
    )
    
    # Calculate maximum number of ticks supported based on label size + padding fraction 
    if label_fraction > 0:
        max_ticks = int(np.floor((1.0 - padding - label_fraction) / label_fraction))
        #if max_ticks <= 3:
        #    max_ticks = 5
        #    if verbose == True:
        #        print('Must be at least two tick marks after pruning ends. Ticks set to 2 after pruning.')
    else:
        raise ZeroDivisionError("Size of label for single tick is zero.")

    # Autogenerate ticks using standard Matplotlib MaxNLocator algorithm
    locator = MaxNLocator(nbins=max(2, max_ticks - 1), 
                          steps=[1, 2, 2.5, 5, 10],
                          prune = 'both')

    
    tick_locs = locator.tick_values(vmin, vmax)

    # Apply new locator and formatters back to the axis object
    axis_obj.set_major_locator(locator)
    axis_obj.set_major_formatter(formatter)

    # Force a redraw to finalize string label formatting
    #fig_copy.canvas.draw_idle()

    formatted_labels = [formatter(loc) for loc in tick_locs]

    plt.close(fig_copy)

    return tick_locs, formatted_labels