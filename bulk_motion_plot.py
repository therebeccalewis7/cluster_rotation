import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.table import Table
from scipy.stats import binned_statistic

def load_cluster_members(fits_file_path,
                         clust_mem_min = 0.3,
                         clust_mem_col1 = 'Membership',
                         clust_mem_col2 = 'Membership_color'):
    """
    """
    # Load in the fits file containing the catalog
    tab = Table.read(fits_file_path)

    # Determine which stars are cluster members by the given cut
    clust_membership = tab['Membership'].data * tab['Membership_color'].data

    # Trim the catalog
    tab_bool = clust_membership >= clust_mem_min
    tab_trimmed = tab[tab_bool]

    return tab_trimmed

def plot_star_velocities(
    fits_file,
    radec_col="x_2010_f153m",
    dec_col="y_2010_f153m",
    pm_radec_col="fit_vx",
    pm_dec_col="fit_vy",
    pm_radec_err_col="fit_vxe",
    pm_dec_err_col="fit_vye",
    clust_mem_min=0.3,
    plate_scale=121.625,
    radec_center=335,
    dec_center=499,
    scale=None,
    key_length=0.1,
    key_label="0.1 mas/yr",
):
    """Plots a vector/quiver field of star proper motions or velocities.

    Parameters:
    -----------
    fits_file : str
        Path to the FITS catalog file.
    radec_col : str
        Name of the column containing RA (or X position).
    dec_col : str
        Name of the column containing Dec (or Y position).
    pm_radec_col : str
        Name of the column containing Proper Motion in RA (or X velocity).
    pm_dec_col : str
        Name of the column containing Proper Motion in Dec (or Y velocity).
    pm_radec_err_col : str
        Name of the column containing Proper Motion RA error.
    pm_dec_err_col : str
        Name of the column containing Proper Motion Dec error.
    clust_mem_min : float, optional
        Minimum cluster membership probability threshold passed to load_cluster_members.
    plate_scale : float, optional
        Scaling factor to convert position and velocity units.
    radec_center, dec_center : float, optional
        Reference origin offsets for the coordinates.
    scale : float, optional
        Scaling factor for vector lengths in quiver display (larger = shorter arrows).
    key_length : float, optional
        Length of the reference arrow in vector units (default is 1.0).
    key_label : str, optional
        Text label for the reference arrow in the quiver key.
    """
    # Load filtered catalog using existing function
    catalog = load_cluster_members(fits_file, clust_mem_min)

    # Extract coordinates and velocities, applying the plate_scale factor
    x = (catalog[radec_col] - radec_center) * plate_scale / 1000.0
    y = (catalog[dec_col] - dec_center) * plate_scale / 1000.0
    u = catalog[pm_radec_col] * plate_scale / 1000.0
    v = catalog[pm_dec_col] * plate_scale / 1000.0
    #u_err = catalog[pm_radec_err_col] * plate_scale / 1000.0
    #v_err = catalog[pm_dec_err_col] * plate_scale / 1000.0
    clust_membership = catalog['Membership'].data * catalog['Membership_color'].data

    u_centered = u - np.average(u, weights=1/clust_membership)
    v_centered = v - np.average(v, weights=1/clust_membership)
    print(np.average(u, weights=1/clust_membership), np.average(v, weights=1/clust_membership))

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create quiver plot
    quiver = ax.quiver(
        x,
        y,
        u_centered,
        v_centered,
        color="teal",
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.8,
        width=0.003,
    )

    # Add Quiver Key (Legend)
    # X=0.9, Y=1.03 places the reference key in axes-relative coordinates just above top-right
    ax.quiverkey(
        quiver,
        X=0.9,
        Y=0.9,
        U=key_length * plate_scale/1000.,
        label=key_label,
        labelpos="N",
        coordinates="axes",
        color = 'r'
    )

    # Set axis labels and plot formatting
    ax.set_xlabel("RAcosDEC")
    ax.set_ylabel("DEC")
    ax.set_title(
        f"Star Vector Field (Membership ≥ {clust_mem_min})", pad=20
    )  # pad added so title doesn't overlap key
    ax.grid(True, linestyle="--", alpha=0.5)

    # Invert X-axis for sky projection alignment
    ax.invert_xaxis()

    plt.tight_layout()
    plt.show()
    fig.savefig('bulk_motion_centered_on_weighted_mean.png', dpi = 300)

from scipy.stats import binned_statistic_2d


def plot_gridded_star_velocities(
    fits_file,
    radec_col="x_2010_f153m",
    dec_col="y_2010_f153m",
    pm_radec_col="fit_vx",
    pm_dec_col="fit_vy",
    pm_radec_err_col="fit_vxe",
    pm_dec_err_col="fit_vye",
    clust_mem_min=0.3,
    plate_scale=121.625,
    radec_center=335,
    dec_center=499,
    scale=None,
    key_length=0.005,
    key_label="0.005 mas/yr",
    nbins_x=10,
    nbins_y=10,
    min_stars_per_bin=1,
    weighted=True,
):
    """Grids spatial coordinates and plots binned average star proper motions.

    Parameters:
    -----------
    nbins_x, nbins_y : int, optional
        Number of grid bins along the X (RA) and Y (Dec) axes (default is 10x10).
    min_stars_per_bin : int, optional
        Minimum number of stars required in a grid cell to plot a vector (default is 1).
    weighted : bool, optional
        If True, computes membership-weighted averages for each grid bin.
        If False, computes standard arithmetic means per bin.
    """
    # Load filtered catalog using existing function
    catalog = load_cluster_members(fits_file, clust_mem_min)

    # Extract coordinates and velocities, applying the plate_scale factor
    conversion_factor = plate_scale / 1000.0
    x = (catalog[radec_col] - radec_center) * conversion_factor
    y = (catalog[dec_col] - dec_center) * conversion_factor
    u = catalog[pm_radec_col] * conversion_factor
    v = catalog[pm_dec_col] * conversion_factor

    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )

    # Center bulk motion relative to overall weighted mean
    weights_overall = 1.0 / clust_membership
    u_centered = u - np.average(u, weights=weights_overall)
    v_centered = v - np.average(v, weights=weights_overall)

    # Define 2D spatial grid bins
    x_bins = np.linspace(x.min(), x.max(), nbins_x + 1)
    y_bins = np.linspace(y.min(), y.max(), nbins_y + 1)

    if weighted:
        # Bin weighted velocities and weights using 2D binned statistics
        u_weighted_sum, _, _, _ = binned_statistic_2d(
            x, y, u_centered * weights_overall, statistic="sum", bins=[x_bins, y_bins]
        )
        v_weighted_sum, _, _, _ = binned_statistic_2d(
            x, y, v_centered * weights_overall, statistic="sum", bins=[x_bins, y_bins]
        )
        weight_sum, _, _, _ = binned_statistic_2d(
            x, y, weights_overall, statistic="sum", bins=[x_bins, y_bins]
        )
        counts, _, _, _ = binned_statistic_2d(
            x, y, None, statistic="count", bins=[x_bins, y_bins]
        )

        # Calculate weighted average per bin (avoiding divide-by-zero)
        with np.errstate(divide="ignore", invalid="ignore"):
            grid_u = u_weighted_sum / weight_sum
            grid_v = v_weighted_sum / weight_sum
    else:
        # Standard unweighted mean per bin
        grid_u, _, _, _ = binned_statistic_2d(
            x, y, u_centered, statistic="mean", bins=[x_bins, y_bins]
        )
        grid_v, _, _, _ = binned_statistic_2d(
            x, y, v_centered, statistic="mean", bins=[x_bins, y_bins]
        )
        counts, _, _, _ = binned_statistic_2d(
            x, y, None, statistic="count", bins=[x_bins, y_bins]
        )

    # Compute grid cell center coordinates
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    y_centers = 0.5 * (y_bins[:-1] + y_bins[1:])
    grid_x, grid_y = np.meshgrid(x_centers, y_centers, indexing="ij")

    # Mask out empty bins or bins with fewer than min_stars_per_bin stars
    valid_mask = counts >= min_stars_per_bin
    plot_x = grid_x[valid_mask]
    plot_y = grid_y[valid_mask]
    plot_u = grid_u[valid_mask]
    plot_v = grid_v[valid_mask]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create quiver plot of binned averages
    quiver = ax.quiver(
        plot_x,
        plot_y,
        plot_u,
        plot_v,
        color="teal",
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.9,
        width=0.005,  # Thicker arrows since there are fewer vectors
    )

    # Add Quiver Key (Legend)
    ax.quiverkey(
        quiver,
        X=0.9,
        Y=0.9,
        U=key_length * conversion_factor,
        label=key_label,
        labelpos="N",
        coordinates="axes",
        color="r",
    )

    # Set axis labels and formatting
    ax.set_xlabel("RAcosDEC", fontsize = 16)
    ax.set_ylabel("DEC", fontsize = 16)
    ax.tick_params(labelsize = 14)
    ax.set_title(
        f"Gridded Star Velocity Field ({nbins_x}x{nbins_y} bins, Membership ≥ {clust_mem_min})",
        pad=20,
    )

    # Set axis ticks directly to bin boundaries
    ax.set_xticks(x_bins)
    ax.set_yticks(y_bins)

    # Set axis limits to fit the grid perfectly
    ax.set_xlim(x_bins.min(), x_bins.max())
    ax.set_ylim(y_bins.min(), y_bins.max())
    
    ax.grid(True, which="both", linestyle="--", alpha=0.6)

    # Invert X-axis for sky projection alignment
    ax.invert_xaxis()

    plt.tight_layout()
    plt.show()
    fig.savefig("binned_bulk_motion.png", dpi=300)

def plot_azimuthal_star_velocities(
    fits_file,
    radec_col="x_2010_f153m",
    dec_col="y_2010_f153m",
    pm_radec_col="fit_vx",
    pm_dec_col="fit_vy",
    pm_radec_err_col="fit_vxe",
    pm_dec_err_col="fit_vye",
    clust_mem_min=0.3,
    plate_scale=121.625,
    radec_center=335,
    dec_center=499,
    scale=None,
    key_length=0.005,
    key_label="0.005 mas/yr",
    nbins_r=5,
    nbins_theta=12,
    r_max=None,
    min_stars_per_bin=1,
    weighted=True,
):
    """Grids spatial coordinates in polar/azimuthal wedges and plots binned average star proper motions.

    Parameters:
    -----------
    nbins_r : int, optional
        Number of radial concentric rings (default is 5).
    nbins_theta : int, optional
        Number of azimuthal sector wedges across 360 degrees (default is 12, i.e. 30-degree wedges).
    r_max : float, optional
        Maximum radius for the outer bin boundary. Defaults to maximum star distance from center.
    min_stars_per_bin : int, optional
        Minimum number of stars required in a wedge to plot a vector.
    weighted : bool, optional
        If True, computes membership-weighted averages for each wedge.
    """
    # Load filtered catalog using existing function
    catalog = load_cluster_members(fits_file, clust_mem_min)

    # Extract coordinates and velocities, applying the plate_scale factor
    conversion_factor = plate_scale / 1000.0
    x = (catalog[radec_col] - radec_center) * conversion_factor
    y = (catalog[dec_col] - dec_center) * conversion_factor
    u = catalog[pm_radec_col] * conversion_factor
    v = catalog[pm_dec_col] * conversion_factor

    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )

    # Center bulk motion relative to overall weighted mean
    weights_overall = clust_membership
    u_centered = u - np.average(u, weights=weights_overall)
    v_centered = v - np.average(v, weights=weights_overall)

    # Convert Cartesian (x, y) to Polar (r, theta) coordinates
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)  # Angles in radians from -pi to +pi

    # Define polar bin edges
    r_upper = r_max if r_max is not None else r.max()
    r_bins = np.linspace(0, r_upper, nbins_r + 1)
    theta_bins = np.linspace(-np.pi, np.pi, nbins_theta + 1)

    # Perform 2D binned statistics in (r, theta) space
    if weighted:
        u_weighted_sum, _, _, _ = binned_statistic_2d(
            r,
            theta,
            u_centered * weights_overall,
            statistic="sum",
            bins=[r_bins, theta_bins],
        )
        v_weighted_sum, _, _, _ = binned_statistic_2d(
            r,
            theta,
            v_centered * weights_overall,
            statistic="sum",
            bins=[r_bins, theta_bins],
        )
        weight_sum, _, _, _ = binned_statistic_2d(
            r, theta, weights_overall, statistic="sum", bins=[r_bins, theta_bins]
        )
        counts, _, _, _ = binned_statistic_2d(
            r, theta, None, statistic="count", bins=[r_bins, theta_bins]
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            grid_u = u_weighted_sum / weight_sum
            grid_v = v_weighted_sum / weight_sum
    else:
        grid_u, _, _, _ = binned_statistic_2d(
            r, theta, u_centered, statistic="mean", bins=[r_bins, theta_bins]
        )
        grid_v, _, _, _ = binned_statistic_2d(
            r, theta, v_centered, statistic="mean", bins=[r_bins, theta_bins]
        )
        counts, _, _, _ = binned_statistic_2d(
            r, theta, None, statistic="count", bins=[r_bins, theta_bins]
        )

    # Compute midpoints for each wedge in (r, theta) to place the vectors
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
    grid_r, grid_theta = np.meshgrid(r_centers, theta_centers, indexing="ij")

    # Convert wedge centerpoints back to Cartesian coordinates for plotting
    grid_x = grid_r * np.cos(grid_theta)
    grid_y = grid_r * np.sin(grid_theta)

    # Mask low-density bins
    valid_mask = counts >= min_stars_per_bin
    plot_x = grid_x[valid_mask]
    plot_y = grid_y[valid_mask]
    plot_u = grid_u[valid_mask]
    plot_v = grid_v[valid_mask]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot azimuthal grid outlines (concentric circles & radial spokes)
    for r_edge in r_bins:
        circle = plt.Circle(
            (0, 0), r_edge, color="gray", linestyle="--", fill=False, alpha=0.5
        )
        ax.add_patch(circle)

    for t_edge in theta_bins:
        ax.plot(
            [0, r_upper * np.cos(t_edge)],
            [0, r_upper * np.sin(t_edge)],
            color="gray",
            linestyle="--",
            alpha=0.5,
        )

    # Plot velocity vectors at wedge centers
    quiver = ax.quiver(
        plot_x,
        plot_y,
        plot_u,
        plot_v,
        color="teal",
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.9,
        width=0.005,
    )

    # Add Quiver Key
    ax.quiverkey(
        quiver,
        X=0.9,
        Y=0.9,
        U=key_length * conversion_factor,
        label=key_label,
        labelpos="N",
        coordinates="axes",
        color="r",
    )

    # Center marker
    ax.plot(0, 0, "r+", markersize=10, label="Center")

    # Set boundaries and formatting
    ax.set_aspect("equal")
    ax.set_xlim(-r_upper * 1.05, r_upper * 1.05)
    ax.set_ylim(-r_upper * 1.05, r_upper * 1.05)

    ax.set_xlabel("RAcosDEC", fontsize = 16)
    ax.set_ylabel("DEC", fontsize = 16)
    ax.tick_params(labelsize = 14)
    ax.set_title(
        f"Azimuthal Binned Vector Field ({nbins_r} rings, {nbins_theta} sectors); Rui+2019 (0,0)",
        pad=20,
    )

    # Invert X-axis for sky projection alignment
    ax.invert_xaxis()

    plt.tight_layout()
    plt.show()
    fig.savefig("azimuthal_binned_motion.png", dpi=300)


def plot_azimuthal_star_velocities_sigma_old(
    fits_file,
    radec_col="x_2010_f153m",
    dec_col="y_2010_f153m",
    radec_err_col="xe_2010_f153m",
    dec_err_col="ye_2010_f153m",
    pm_radec_col="fit_vx",
    pm_dec_col="fit_vy",
    pm_radec_err_col="fit_vxe",
    pm_dec_err_col="fit_vye",
    clust_mem_min=0.3,
    plate_scale=121.625,
    radec_center=None,
    dec_center=None,
    scale=None,
    key_length=3.0,  # Default magnitude for key label in units of sigma
    key_label=r"$3\sigma$",
    nbins_r=5,
    nbins_theta=12,
    r_max=None,
    min_stars_per_bin=3,
    weighted=True,
    error_type="std",  # Options: 'mean' (error on mean) or 'std' (standard deviation)
):
    """Grids spatial coordinates in polar/azimuthal wedges and plots binned average 

    star proper motions expressed in units of sigma (average / uncertainty).

    Parameters:
    -----------
    error_type : str, optional
        Uncertainty definition used for normalizing vector magnitude:
        - 'mean': Standard error of the mean (default).
        - 'std': Standard deviation of stellar velocities within each bin.
    """
    # Load filtered catalog using existing function
    catalog = load_cluster_members(fits_file, clust_mem_min)

    if (radec_center is None) or (dec_center is None):
        radec_center = np.average(
            catalog[radec_col], weights=1 / (catalog[radec_err_col] ** 2)
        )
        dec_center = np.average(
            catalog[dec_col], weights=1 / (catalog[dec_err_col] ** 2)
        )

    # Extract coordinates and velocities, applying the plate_scale factor
    conversion_factor = plate_scale / 1000.0
    x = (catalog[radec_col] - radec_center) * conversion_factor
    y = (catalog[dec_col] - dec_center) * conversion_factor
    u = catalog[pm_radec_col] * conversion_factor
    v = catalog[pm_dec_col] * conversion_factor

    # Extract errors (scaled appropriately)
    u_err = catalog[pm_radec_err_col] * conversion_factor
    v_err = catalog[pm_dec_err_col] * conversion_factor

    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )

    # Center bulk motion relative to overall weighted mean
    weights_overall = clust_membership if weighted else np.ones_like(u)
    u_centered = u - np.average(u, weights=weights_overall)
    v_centered = v - np.average(v, weights=weights_overall)

    # Convert Cartesian (x, y) to Polar (r, theta) coordinates
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)  # Angles in radians from -pi to +pi

    # Define polar bin edges
    r_upper = r_max if r_max is not None else r.max()
    r_bins = np.linspace(0, r_upper, nbins_r + 1)
    theta_bins = np.linspace(-np.pi, np.pi, nbins_theta + 1)

    # Perform 2D binned statistics in (r, theta) space
    counts, _, _, _ = binned_statistic_2d(
        r, theta, None, statistic="count", bins=[r_bins, theta_bins]
    )

    if weighted:
        u_weighted_sum, _, _, _ = binned_statistic_2d(
            r, theta, u_centered * weights_overall, statistic="sum", bins=[r_bins, theta_bins]
        )
        v_weighted_sum, _, _, _ = binned_statistic_2d(
            r, theta, v_centered * weights_overall, statistic="sum", bins=[r_bins, theta_bins]
        )
        weight_sum, _, _, _ = binned_statistic_2d(
            r, theta, weights_overall, statistic="sum", bins=[r_bins, theta_bins]
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            grid_u_mean = u_weighted_sum / weight_sum
            grid_v_mean = v_weighted_sum / weight_sum

            if error_type == "std":
                # Weighted variance calculation in 2D bins: Var = <v^2>_w - (<v>_w)^2
                u_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r, theta, (u_centered ** 2) * weights_overall, statistic="sum", bins=[r_bins, theta_bins]
                )
                v_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r, theta, (v_centered ** 2) * weights_overall, statistic="sum", bins=[r_bins, theta_bins]
                )
                u_var = (u_sq_w_sum / weight_sum) - (grid_u_mean ** 2)
                v_var = (v_sq_w_sum / weight_sum) - (grid_v_mean ** 2)

                grid_u_sigma_err = np.sqrt(np.maximum(0.0, u_var))
                grid_v_sigma_err = np.sqrt(np.maximum(0.0, v_var))
            else:
                # Standard error of the weighted mean
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r, theta, (weights_overall * u_err) ** 2, statistic="sum", bins=[r_bins, theta_bins]
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r, theta, (weights_overall * v_err) ** 2, statistic="sum", bins=[r_bins, theta_bins]
                )
                grid_u_sigma_err = np.sqrt(u_err_sq_sum) / weight_sum
                grid_v_sigma_err = np.sqrt(v_err_sq_sum) / weight_sum
    else:
        grid_u_mean, _, _, _ = binned_statistic_2d(
            r, theta, u_centered, statistic="mean", bins=[r_bins, theta_bins]
        )
        grid_v_mean, _, _, _ = binned_statistic_2d(
            r, theta, v_centered, statistic="mean", bins=[r_bins, theta_bins]
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            if error_type == "std":
                # Unweighted standard deviation in 2D bins
                grid_u_sigma_err, _, _, _ = binned_statistic_2d(
                    r, theta, u_centered, statistic="std", bins=[r_bins, theta_bins]
                )
                grid_v_sigma_err, _, _, _ = binned_statistic_2d(
                    r, theta, v_centered, statistic="std", bins=[r_bins, theta_bins]
                )
            else:
                # Standard error of the unweighted mean
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r, theta, u_err ** 2, statistic="sum", bins=[r_bins, theta_bins]
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r, theta, v_err ** 2, statistic="sum", bins=[r_bins, theta_bins]
                )
                grid_u_sigma_err = np.sqrt(u_err_sq_sum) / counts
                grid_v_sigma_err = np.sqrt(v_err_sq_sum) / counts

    # Convert velocity components to units of sigma
    with np.errstate(divide="ignore", invalid="ignore"):
        grid_u = grid_u_mean / grid_u_sigma_err
        grid_v = grid_v_mean / grid_v_sigma_err

    # Compute midpoints for each wedge in (r, theta) to place the vectors
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
    grid_r, grid_theta = np.meshgrid(r_centers, theta_centers, indexing="ij")

    # Convert wedge centerpoints back to Cartesian coordinates for plotting
    grid_x = grid_r * np.cos(grid_theta)
    grid_y = grid_r * np.sin(grid_theta)

    # Mask low-density bins and invalid divisions
    valid_mask = (counts >= min_stars_per_bin) & np.isfinite(grid_u) & np.isfinite(grid_v)
    plot_x = grid_x[valid_mask]
    plot_y = grid_y[valid_mask]
    plot_u = grid_u[valid_mask]
    plot_v = grid_v[valid_mask]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot azimuthal grid outlines (concentric circles & radial spokes)
    for r_edge in r_bins:
        circle = plt.Circle(
            (0, 0), r_edge, color="gray", linestyle="--", fill=False, alpha=0.5
        )
        ax.add_patch(circle)

    for t_edge in theta_bins:
        ax.plot(
            [0, r_upper * np.cos(t_edge)],
            [0, r_upper * np.sin(t_edge)],
            color="gray",
            linestyle="--",
            alpha=0.5,
        )

    # Plot velocity vectors at wedge centers
    quiver = ax.quiver(
        plot_x,
        plot_y,
        plot_u,
        plot_v,
        color="k",
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.9,
        width=0.005,
    )

    # Add Quiver Key (U is now directly in units of sigma)
    ax.quiverkey(
        quiver,
        X=0.9,
        Y=0.9,
        U=key_length,
        label=key_label,
        labelpos="N",
        coordinates="axes",
        color="r",
    )

    # Center marker
    ax.plot(0, 0, "r+", markersize=10, label="Center")

    # Set boundaries and formatting
    ax.set_aspect("equal")
    ax.set_xlim(-r_upper * 1.05, r_upper * 1.05)
    ax.set_ylim(-r_upper * 1.05, r_upper * 1.05)

    ax.set_xlabel("RAcosDEC", fontsize=16)
    ax.set_ylabel("DEC", fontsize=16)
    ax.tick_params(labelsize=14)

    # Invert X-axis for sky projection alignment
    ax.invert_xaxis()

    plt.tight_layout()
    plt.show()
    fig.savefig("azimuthal_binned_motion_sigma.png", dpi=300)





