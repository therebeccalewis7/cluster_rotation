import matplotlib.pyplot as plt
import numpy as np
from cluster_rotation_analysis import clust_rot_utils as cru
import matplotlib.colors as colors
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import binned_statistic_2d
from astropy.table import Table



# def plot_individual_stars_pm_grid(
#     fits_file,
#     radec_col="x_2010_f153m",
#     dec_col="y_2010_f153m",
#     radec_err_col="xe_2010_f153m",
#     dec_err_col="ye_2010_f153m",
#     pm_radec_col="fit_vx",
#     pm_dec_col="fit_vy",
#     pm_radec_err_col="fit_vxe",
#     pm_dec_err_col="fit_vye",
#     clust_mem_min=0.3,
#     plate_scale=121.625,
#     radec_center=None,
#     dec_center=None,
#     cmap="coolwarm",
#     markersize=15,
#     alpha=1.0,
# ):
#     """
#     Plots individual cluster member stars in a 2x2 grid of spatial scatter plots
#     colored by X and Y proper motions in mas/yr (top row) and in units of sigma (bottom row),
#     with colorbars centered at zero.
#     """
#     # Load cluster catalog
#     catalog = cru.load_cluster_members(fits_file, clust_mem_min)

#     # Determine center coordinates if not provided
#     if (radec_center is None) or (dec_center is None):
#         radec_center = np.average(
#             catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
#         )
#         dec_center = np.average(
#             catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
#         )

#     # Extract coordinates and velocities with plate scale conversion (mas/yr)
#     conversion_factor = plate_scale
#     x = (catalog[radec_col] - radec_center) * conversion_factor
#     y = (catalog[dec_col] - dec_center) * conversion_factor

#     u = catalog[pm_radec_col] * conversion_factor
#     v = catalog[pm_dec_col] * conversion_factor
#     u_err = catalog[pm_radec_err_col] * conversion_factor
#     v_err = catalog[pm_dec_err_col] * conversion_factor

#     # Center bulk motion relative to overall weighted mean
#     clust_membership = (
#         catalog["Membership"].data * catalog["Membership_color"].data
#     )
#     weights = clust_membership

#     u_centered = u - np.average(u, weights=weights)
#     v_centered = v - np.average(v, weights=weights)

#     # Calculate velocities in units of sigma
#     u_sigma = u_centered / u_err
#     v_sigma = v_centered / v_err

#     # Make figure
#     fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

#     # Configure panels data and formatting
#     panels = [
#         (
#             axes[0, 0],
#             u_centered,
#             r"$\mu_x$ [mas/yr]",
#             "X Proper Motion [mas/yr]",
#         ),
#         (
#             axes[0, 1],
#             v_centered,
#             r"$\mu_y$ [mas/yr]",
#             "Y Proper Motion [mas/yr]",
#         ),
#         (
#             axes[1, 0],
#             u_sigma,
#             r"$\mu_x / \sigma_x$",
#             "X Proper Motion [$\sigma$]",
#         ),
#         (
#             axes[1, 1],
#             v_sigma,
#             r"$\mu_y / \sigma_y$",
#             "Y Proper Motion [$\sigma$]",
#         ),
#     ]

#     for ax, color_data, cbar_label, title in panels:
#         # Compute symmetric limits centered at 0
#         max_abs_val = np.nanpercentile(np.abs(color_data), 95)
#         norm = TwoSlopeNorm(
#             vmin=-max_abs_val, vcenter=0.0, vmax=max_abs_val
#         )

#         # Scatter plot for individual stars
#         sc = ax.scatter(
#             x / 1000.0,
#             y / 1000.0,
#             c=color_data,
#             cmap=cmap,
#             norm=norm,
#             s=markersize,
#             alpha=alpha,
#             marker="o",
#             edgecolors="none",
#         )

#         # Plot center marker
#         ax.plot(0, 0, "r+", markersize=10, markeredgewidth=1.5, label="Center")

#         # Colorbar
#         cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
#         cbar.set_label(cbar_label, fontsize=12)
#         cbar.ax.tick_params(labelsize = 12)

#         ax.set_title(title, fontsize=16, pad=10)
#         ax.set_aspect("equal")
#         ax.grid(True, linestyle=":", alpha=0.5)

#         # Invert X-axis for sky projection alignment
#         ax.invert_xaxis()

#         ax.tick_params(labelsize=16)

#     # Set outer axis labels
#     for ax in axes[1, :]:
#         ax.set_xlabel("RAcosDEC [arcsec]", fontsize=16)
#     for ax in axes[:, 0]:
#         ax.set_ylabel("DEC [arcsec]", fontsize=16)

#     plt.tight_layout()
#     plt.show()
#     fig.savefig("star_pm_scatter_grid.png", dpi=300)

#     return fig, axes

# def plot_pm_density_grid(
#     fits_file,
#     radec_col="x_2010_f153m",
#     dec_col="y_2010_f153m",
#     radec_err_col="xe_2010_f153m",
#     dec_err_col="ye_2010_f153m",
#     pm_radec_col="fit_vx",
#     pm_dec_col="fit_vy",
#     pm_radec_err_col="fit_vxe",
#     pm_dec_err_col="fit_vye",
#     clust_mem_min=0.3,
#     plate_scale=121.625,
#     radec_center=335,
#     dec_center=499,
#     cmap="coolwarm",
#     nbins=10,
#     min_stars_per_bin=3,
#     weighted=True,
# ):
#     """Plots a 2x2 grid of spatial density maps showing binned average X and Y proper motions

#     in mas/yr (top row) and in units of sigma (bottom row), with colorbars centered at zero.
#     """
#     # Load cluster catalog
#     catalog = cru.load_cluster_members(fits_file, clust_mem_min)

#     # Determine center coordinates if not provided
#     if (radec_center is None) or (dec_center is None):
#         radec_center = np.average(
#             catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
#         )
#         dec_center = np.average(
#             catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
#         )

#     # Extract coordinates and velocities with plate scale conversion (mas)
#     conversion_factor = plate_scale
#     x = (catalog[radec_col] - radec_center) * conversion_factor / 1000.0
#     y = (catalog[dec_col] - dec_center) * conversion_factor / 1000.0

#     u = catalog[pm_radec_col] * conversion_factor
#     v = catalog[pm_dec_col] * conversion_factor
#     u_err = catalog[pm_radec_err_col] * conversion_factor
#     v_err = catalog[pm_dec_err_col] * conversion_factor

#     # Center bulk motion relative to overall weighted mean
#     clust_membership = (
#         catalog["Membership"].data * catalog["Membership_color"].data
#     )
#     weights = clust_membership if weighted else np.ones_like(u)

#     u_centered = u - np.average(u, weights=weights)
#     v_centered = v - np.average(v, weights=weights)

#     # Calculate velocities in units of sigma
#     u_sigma = u_centered / u_err
#     v_sigma = v_centered / v_err

#     # Set up 2D spatial bin grid edges
#     x_bins = np.linspace(x.min(), x.max(), nbins + 1)
#     y_bins = np.linspace(y.min(), y.max(), nbins + 1)

#     # Compute star counts per 2D spatial bin
#     counts, _, _, _ = binned_statistic_2d(
#         x, y, None, statistic="count", bins=[x_bins, y_bins]
#     )

#     # Helper to compute weighted or unweighted binned 2D mean map
#     def compute_binned_map(values):
#         if weighted:
#             w_sum, _, _, _ = binned_statistic_2d(
#                 x, y, values * weights, statistic="sum", bins=[x_bins, y_bins]
#             )
#             weight_sum, _, _, _ = binned_statistic_2d(
#                 x, y, weights, statistic="sum", bins=[x_bins, y_bins]
#             )
#             with np.errstate(divide="ignore", invalid="ignore"):
#                 grid_map = w_sum / weight_sum
#         else:
#             grid_map, _, _, _ = binned_statistic_2d(
#                 x, y, values, statistic="mean", bins=[x_bins, y_bins]
#             )

#         # Mask bins with fewer than min_stars_per_bin
#         grid_map[counts < min_stars_per_bin] = np.nan
#         return grid_map

#     # Generate 2D grid maps for each component
#     grid_u = compute_binned_map(u_centered)
#     grid_v = compute_binned_map(v_centered)
#     grid_u_sigma = compute_binned_map(u_sigma)
#     grid_v_sigma = compute_binned_map(v_sigma)

#     # Make figure
#     fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

#     # Configure panels data and formatting
#     panels = [
#         (
#             axes[0, 0],
#             grid_u,
#             r"$\mu_x$ [mas/yr]",
#             "X Proper Motion [mas/yr]",
#         ),
#         (
#             axes[0, 1],
#             grid_v,
#             r"$\mu_y$ [mas/yr]",
#             "Y Proper Motion [mas/yr]",
#         ),
#         (
#             axes[1, 0],
#             grid_u_sigma,
#             r"$\mu_x / \sigma_x$",
#             "X Proper Motion [$\sigma$]",
#         ),
#         (
#             axes[1, 1],
#             grid_v_sigma,
#             r"$\mu_y / \sigma_y$",
#             "Y Proper Motion [$\sigma$]",
#         ),
#     ]

#     for ax, grid_data, cbar_label, title in panels:
#         # Compute symmetric limits centered at 0
#         max_abs_val = np.nanpercentile(np.abs(grid_data), 95)
#         norm = TwoSlopeNorm(
#             vmin=-max_abs_val, vcenter=0.0, vmax=max_abs_val
#         )

#         # Plot 2D density grid map
#         im = ax.pcolormesh(
#             x_bins,
#             y_bins,
#             grid_data.T,  # Transpose to align (x, y) orientation
#             cmap=cmap,
#             norm=norm,
#             shading="flat",
#         )

#         # Plot center marker
#         ax.plot(0, 0, "r+", markersize=10, markeredgewidth=1.5, label="Center")

#         # Colorbar
#         cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
#         cbar.set_label(cbar_label, fontsize=12)
#         cbar.ax.tick_params(labelsize=12)

#         ax.set_title(title, fontsize=16, pad=10)
#         ax.set_aspect("equal")
#         ax.grid(True, linestyle=":", alpha=0.5)

#         # Invert X-axis for sky projection alignment
#         ax.invert_xaxis()

#         ax.tick_params(labelsize=16)

#     # Set outer axis labels
#     for ax in axes[1, :]:
#         ax.set_xlabel("RAcosDEC [arcsec]", fontsize=16)
#     for ax in axes[:, 0]:
#         ax.set_ylabel("DEC [arcsec]", fontsize=16)

#     plt.tight_layout()
#     plt.show()
#     fig.savefig("star_pm_density_grid.png", dpi=300)

#     return fig, axes


# def transform_pixel_to_sky(
#     fits_file, 
#     x_col="x_2010_f153m", 
#     y_col="y_2010_f153m", 
#     x_err_col="xe_2010_f153m",
#     y_err_col="ye_2010_f153m",
#     pm_x_col="fit_vx", 
#     pm_y_col="fit_vy",
#     pm_x_err_col="fit_vxe", 
#     pm_y_err_col="fit_vye",
#     pa_deg=0.0, 
#     pixel_scale=121.625, 
#     scale_divisor=1000.0,
#     pixel_east_left=False
# ):
#     """
#     Reads a FITS file into an Astropy Table, rotates positions, proper motions, 
#     and their uncertainties by a position angle (degrees East of North), and 
#     scales them to on-sky coordinates where East is to the left (-X) and North is up (+Y).
    
#     Parameters:
#     -----------
#     fits_file : str
#         Path to the FITS file.
#     x_col, y_col : str
#         Column names for the x and y pixel coordinates.
#     x_err_col, y_err_col : str, optional
#         Column names for positional uncertainties (set to None if not needed).
#     pm_x_col, pm_y_col : str, optional
#         Column names for the x and y proper motions (set to None if not needed).
#     pm_x_err_col, pm_y_err_col : str, optional
#         Column names for the x and y proper motion uncertainties.
#     pa_deg : float
#         Position angle in degrees East of North.
#     pixel_scale : float
#         The scale factor to convert pixels to physical units (e.g., mas/pixel).
#     scale_divisor : float
#         Factor to divide the pixel scale by (default 1000.0 to convert mas to arcsec).
#     pixel_east_left : bool
#         Set to False (default) if raw detector pixels have +X to the right (West).
#         Set to True if raw pixel X-axis is already inverted so +X points East (left).
        
#     Returns:
#     --------
#     catalog : astropy.table.Table
#         The updated catalog with transformed sky positions, proper motions, 
#         and propagated uncertainties.
#     """
#     catalog = Table.read(fits_file, format='fits')
    
#     theta = np.deg2rad(pa_deg)
#     cos_t = np.cos(theta)
#     sin_t = np.sin(theta)
#     conversion_factor = pixel_scale / scale_divisor
#     sign_flip = 1.0 if pixel_east_left else -1.0

#     # 1. Transform Positions
#     if x_col in catalog.colnames and y_col in catalog.colnames:
#         x_east = sign_flip * catalog[x_col]
#         y_north = catalog[y_col]
        
#         x_rot_east = x_east * cos_t + y_north * sin_t
#         y_rot_north = -x_east * sin_t + y_north * cos_t
        
#         catalog['radec_sky'] = -x_rot_east * conversion_factor
#         catalog['dec_sky'] = y_rot_north * conversion_factor

#     # 2. Propagate Position Uncertainties (Quadrature)
#     if x_err_col and y_err_col and (x_err_col in catalog.colnames) and (y_err_col in catalog.colnames):
#         x_err = catalog[x_err_col]
#         y_err = catalog[y_err_col]
        
#         catalog['radec_sky_err'] = np.sqrt((x_err * cos_t)**2 + (y_err * sin_t)**2) * conversion_factor
#         catalog['dec_sky_err'] = np.sqrt((x_err * sin_t)**2 + (y_err * cos_t)**2) * conversion_factor

#     # 3. Transform Proper Motions
#     if pm_x_col and pm_y_col and (pm_x_col in catalog.colnames) and (pm_y_col in catalog.colnames):
#         u_east = sign_flip * catalog[pm_x_col]
#         v_north = catalog[pm_y_col]
        
#         u_rot_east = u_east * cos_t + v_north * sin_t
#         v_rot_north = -u_east * sin_t + v_north * cos_t
        
#         catalog['pm_radec_sky'] = -u_rot_east * conversion_factor
#         catalog['pm_dec_sky'] = v_rot_north * conversion_factor

#     # 4. Propagate Proper Motion Uncertainties (Quadrature)
#     if pm_x_err_col and pm_y_err_col and (pm_x_err_col in catalog.colnames) and (pm_y_err_col in catalog.colnames):
#         u_err = catalog[pm_x_err_col]
#         v_err = catalog[pm_y_err_col]
        
#         catalog['pm_radec_sky_err'] = np.sqrt((u_err * cos_t)**2 + (v_err * sin_t)**2) * conversion_factor
#         catalog['pm_dec_sky_err'] = np.sqrt((u_err * sin_t)**2 + (v_err * cos_t)**2) * conversion_factor
        
#     return catalog


# def transform_pixel_to_sky(
#     catalog_or_file,
#     x_col="x_2010_f153m",
#     y_col="y_2010_f153m",
#     x_err_col="xe_2010_f153m",
#     y_err_col="ye_2010_f153m",
#     pm_x_col="fit_vx",
#     pm_y_col="fit_vy",
#     pm_x_err_col="fit_vxe",
#     pm_y_err_col="fit_vye",
#     pa_deg=0.0,
#     pixel_scale=121.625,
#     scale_divisor=1000.0,
#     pixel_east_left=True,
#     clust_mem_min = 0.3
# ):
#     """
#     Accepts a FITS file path or Astropy Table, rotates positions, proper motions,
#     and their uncertainties by a position angle (degrees East of North), and
#     scales positions to arcsec and proper motions to mas/yr with East to the left (-X).
#     """
#     if isinstance(catalog_or_file, Table):
#         catalog = catalog_or_file.copy()
#     else:
#         catalog = cru.load_cluster_members(catalog_or_file, clust_mem_min = clust_mem_min)

#     theta = np.deg2rad(pa_deg)
#     cos_t = np.cos(theta)
#     sin_t = np.sin(theta)
#     pos_factor = pixel_scale / scale_divisor  # pixels -> arcsec
#     pm_factor = pixel_scale                   # pixels/yr -> mas/yr
#     sign_flip = 1.0 if pixel_east_left else -1.0

#     # 1. Transform Positions
#     if x_col in catalog.colnames and y_col in catalog.colnames:
#         x_east = sign_flip * catalog[x_col]
#         y_north = catalog[y_col]

#         x_rot_east = x_east * cos_t + y_north * sin_t
#         y_rot_north = -x_east * sin_t + y_north * cos_t

#         catalog["radec_sky"] = -x_rot_east * pos_factor
#         catalog["dec_sky"] = y_rot_north * pos_factor

#     # 2. Propagate Position Uncertainties (Quadrature)
#     if (
#         x_err_col
#         and y_err_col
#         and (x_err_col in catalog.colnames)
#         and (y_err_col in catalog.colnames)
#     ):
#         x_err = catalog[x_err_col]
#         y_err = catalog[y_err_col]

#         catalog["radec_sky_err"] = (
#             np.sqrt((x_err * cos_t) ** 2 + (y_err * sin_t) ** 2) * pos_factor
#         )
#         catalog["dec_sky_err"] = (
#             np.sqrt((x_err * sin_t) ** 2 + (y_err * cos_t) ** 2) * pos_factor
#         )

#     # 3. Transform Proper Motions
#     if (
#         pm_x_col
#         and pm_y_col
#         and (pm_x_col in catalog.colnames)
#         and (pm_y_col in catalog.colnames)
#     ):
#         u_east = sign_flip * catalog[pm_x_col]
#         v_north = catalog[pm_y_col]

#         u_rot_east = u_east * cos_t + v_north * sin_t
#         v_rot_north = -u_east * sin_t + v_north * cos_t

#         catalog["pm_radec_sky"] = -u_rot_east * pm_factor
#         catalog["pm_dec_sky"] = v_rot_north * pm_factor

#     # 4. Propagate Proper Motion Uncertainties (Quadrature)
#     if (
#         pm_x_err_col
#         and pm_y_err_col
#         and (pm_x_err_col in catalog.colnames)
#         and (pm_y_err_col in catalog.colnames)
#     ):
#         u_err = catalog[pm_x_err_col]
#         v_err = catalog[pm_y_err_col]

#         catalog["pm_radec_sky_err"] = (
#             np.sqrt((u_err * cos_t) ** 2 + (v_err * sin_t) ** 2) * pm_factor
#         )
#         catalog["pm_dec_sky_err"] = (
#             np.sqrt((u_err * sin_t) ** 2 + (v_err * cos_t) ** 2) * pm_factor
#         )

#     return catalog

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, Galactic
from astropy.table import Table
from cluster_rotation_analysis import clust_rot_utils as cru


def transform_pixel_to_sky(
    catalog_or_file,
    x_col="x_2010_f153m",
    y_col="y_2010_f153m",
    x_err_col="xe_2010_f153m",
    y_err_col="ye_2010_f153m",
    pm_x_col="fit_vx",
    pm_y_col="fit_vy",
    pm_x_err_col="fit_vxe",
    pm_y_err_col="fit_vye",
    pa_deg=0.0,
    pixel_scale=121.625,
    scale_divisor=1000.0,
    pixel_east_left=True,
    clust_mem_min=0.3,
    to_galactic=False,
    ra_zero_deg=None,
    dec_zero_deg=None,
):
    """
    Accepts a FITS file path or Astropy Table, rotates positions, proper motions,
    and their uncertainties by a position angle (degrees East of North), and
    scales positions to arcsec and proper motions to mas/yr with East to the left (-X).

    Optionally converts the sky-projected coordinates and proper motions into
    Galactic coordinates (l, b, pm_l_cosb, pm_b) given the ICRS (RA, Dec) in degrees
    corresponding to pixel (0, 0).
    """
    if isinstance(catalog_or_file, Table):
        catalog = catalog_or_file.copy()
    else:
        catalog = cru.load_cluster_members(catalog_or_file, clust_mem_min=clust_mem_min)

    theta = np.deg2rad(pa_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    pos_factor = pixel_scale / scale_divisor  # pixels -> arcsec
    pm_factor = pixel_scale                   # pixels/yr -> mas/yr
    sign_flip = 1.0 if pixel_east_left else -1.0

    # 1. Transform Positions
    if x_col in catalog.colnames and y_col in catalog.colnames:
        x_east = sign_flip * catalog[x_col]
        y_north = catalog[y_col]

        x_rot_east = x_east * cos_t + y_north * sin_t
        y_rot_north = -x_east * sin_t + y_north * cos_t

        catalog["radec_sky"] = -x_rot_east * pos_factor
        catalog["dec_sky"] = y_rot_north * pos_factor

    # 2. Propagate Position Uncertainties (Quadrature)
    if (
        x_err_col
        and y_err_col
        and (x_err_col in catalog.colnames)
        and (y_err_col in catalog.colnames)
    ):
        x_err = catalog[x_err_col]
        y_err = catalog[y_err_col]

        catalog["radec_sky_err"] = (
            np.sqrt((x_err * cos_t) ** 2 + (y_err * sin_t) ** 2) * pos_factor
        )
        catalog["dec_sky_err"] = (
            np.sqrt((x_err * sin_t) ** 2 + (y_err * cos_t) ** 2) * pos_factor
        )

    # 3. Transform Proper Motions
    if (
        pm_x_col
        and pm_y_col
        and (pm_x_col in catalog.colnames)
        and (pm_y_col in catalog.colnames)
    ):
        u_east = sign_flip * catalog[pm_x_col]
        v_north = catalog[pm_y_col]

        u_rot_east = u_east * cos_t + v_north * sin_t
        v_rot_north = -u_east * sin_t + v_north * cos_t

        catalog["pm_radec_sky"] = -u_rot_east * pm_factor
        catalog["pm_dec_sky"] = v_rot_north * pm_factor

    # 4. Propagate Proper Motion Uncertainties (Quadrature)
    if (
        pm_x_err_col
        and pm_y_err_col
        and (pm_x_err_col in catalog.colnames)
        and (pm_y_err_col in catalog.colnames)
    ):
        u_err = catalog[pm_x_err_col]
        v_err = catalog[pm_y_err_col]

        catalog["pm_radec_sky_err"] = (
            np.sqrt((u_err * cos_t) ** 2 + (v_err * sin_t) ** 2) * pm_factor
        )
        catalog["pm_dec_sky_err"] = (
            np.sqrt((u_err * sin_t) ** 2 + (v_err * cos_t) ** 2) * pm_factor
        )

    # 5. Optional: Convert Sky Offsets & Proper Motions to Galactic Coordinates
    if to_galactic:
        if ra_zero_deg is None or dec_zero_deg is None:
            raise ValueError(
                "ra_zero_deg and dec_zero_deg (the ICRS coordinates of pixel (0, 0) in degrees) "
                "must be provided when to_galactic=True."
            )

        # Tangent plane origin at pixel (0, 0)
        origin = SkyCoord(ra=ra_zero_deg * u.deg, dec=dec_zero_deg * u.deg, frame="icrs")
        offset_frame = SkyOffsetFrame(origin=origin)

        # Note: SkyOffsetFrame lon increases East (+x_rot_east = -radec_sky)
        lon_offset = (-np.asarray(catalog["radec_sky"])) * u.arcsec
        lat_offset = np.asarray(catalog["dec_sky"]) * u.arcsec

        has_pm = "pm_radec_sky" in catalog.colnames and "pm_dec_sky" in catalog.colnames

        if has_pm:
            pm_lon = (-np.asarray(catalog["pm_radec_sky"])) * (u.mas / u.yr)
            pm_lat = np.asarray(catalog["pm_dec_sky"]) * (u.mas / u.yr)

            sky_coord = SkyCoord(
                lon=lon_offset,
                lat=lat_offset,
                pm_lon_coslat=pm_lon,
                pm_lat=pm_lat,
                frame=offset_frame,
            )
        else:
            sky_coord = SkyCoord(lon=lon_offset, lat=lat_offset, frame=offset_frame)

        # Transform from tangent-plane sky offsets -> ICRS -> Galactic
        icrs_coord = sky_coord.transform_to("icrs")
        gal_coord = icrs_coord.transform_to(Galactic())

        # Store absolute ICRS and Galactic positions (degrees)
        catalog["ra_icrs"] = icrs_coord.ra.deg
        catalog["dec_icrs"] = icrs_coord.dec.deg
        catalog["l_gal"] = gal_coord.l.deg
        catalog["b_gal"] = gal_coord.b.deg

        # Also store Galactic tangent-plane offsets relative to pixel (0, 0) in arcsec
        gal_origin = origin.transform_to(Galactic())
        gal_offset = gal_coord.transform_to(SkyOffsetFrame(origin=gal_origin))
        catalog["l_gal_offset"] = gal_offset.lon.to(u.arcsec).value
        catalog["b_gal_offset"] = gal_offset.lat.to(u.arcsec).value

        # Compute parallactic/Galactic rotation angle phi between Equatorial and Galactic axes
        # at each star's position to transform proper motions and propagate uncertainties
        icrs_north = SkyCoord(
            ra=icrs_coord.ra,
            dec=np.clip(icrs_coord.dec + 1.0 * u.arcsec, -90 * u.deg, 90 * u.deg),
            frame="icrs",
        ).transform_to(Galactic())
        phi = gal_coord.position_angle(icrs_north).rad

        cos_p = np.cos(phi)
        sin_p = np.sin(phi)

        if has_pm:
            catalog["pm_l_gal"] = gal_coord.pm_l_cosb.to(u.mas / u.yr).value
            catalog["pm_b_gal"] = gal_coord.pm_b.to(u.mas / u.yr).value

        if "radec_sky_err" in catalog.colnames and "dec_sky_err" in catalog.colnames:
            ra_err = np.asarray(catalog["radec_sky_err"])
            dec_err = np.asarray(catalog["dec_sky_err"])
            catalog["l_gal_err"] = np.sqrt((ra_err * cos_p) ** 2 + (dec_err * sin_p) ** 2)
            catalog["b_gal_err"] = np.sqrt((ra_err * sin_p) ** 2 + (dec_err * cos_p) ** 2)

        if "pm_radec_sky_err" in catalog.colnames and "pm_dec_sky_err" in catalog.colnames:
            pm_ra_err = np.asarray(catalog["pm_radec_sky_err"])
            pm_dec_err = np.asarray(catalog["pm_dec_sky_err"])
            catalog["pm_l_gal_err"] = np.sqrt((pm_ra_err * cos_p) ** 2 + (pm_dec_err * sin_p) ** 2)
            catalog["pm_b_gal_err"] = np.sqrt((pm_ra_err * sin_p) ** 2 + (pm_dec_err * cos_p) ** 2)

    return catalog


# def plot_pm_density_grid(
#     fits_file,
#     radec_col="x_2010_f153m",
#     dec_col="y_2010_f153m",
#     radec_err_col="xe_2010_f153m",
#     dec_err_col="ye_2010_f153m",
#     pm_radec_col="fit_vx",
#     pm_dec_col="fit_vy",
#     pm_radec_err_col="fit_vxe",
#     pm_dec_err_col="fit_vye",
#     clust_mem_min=0.3,
#     plate_scale=121.625,
#     radec_center=335,
#     dec_center=499,
#     cmap="coolwarm",
#     nbins=10,
#     min_stars_per_bin=3,
#     weighted=True,
#     transform_to_sky=True,
#     pa_deg=-42.3,
#     scale_divisor=1.0,
#     pixel_east_left=True,
# ):
#     """Plots a 2x2 grid of spatial density maps showing binned average X and Y proper motions
#     in mas/yr (top row) and in units of sigma (bottom row), with colorbars centered at zero.
#     Optionally applies a position angle rotation and pixel-to-sky coordinate transformation.
#     """
#     # Load cluster catalog
#     catalog = cru.load_cluster_members(fits_file, clust_mem_min)

#     if transform_to_sky:
#         catalog = transform_pixel_to_sky(
#             catalog,
#             x_col=radec_col,
#             y_col=dec_col,
#             x_err_col=radec_err_col,
#             y_err_col=dec_err_col,
#             pm_x_col=pm_radec_col,
#             pm_y_col=pm_dec_col,
#             pm_x_err_col=pm_radec_err_col,
#             pm_y_err_col=pm_dec_err_col,
#             pa_deg=pa_deg,
#             pixel_scale=plate_scale,
#             scale_divisor=scale_divisor,
#             pixel_east_left=pixel_east_left,
#         )

#         # Transform provided pixel centers into the same rotated sky frame
#         if (radec_center is not None) and (dec_center is not None):
#             theta = np.deg2rad(pa_deg)
#             pos_factor = plate_scale / scale_divisor
#             sign_flip = 1.0 if pixel_east_left else -1.0

#             cx_east = sign_flip * radec_center
#             cy_north = dec_center
#             radec_center = -(cx_east * np.cos(theta) + cy_north * np.sin(theta)) * pos_factor
#             dec_center = (-cx_east * np.sin(theta) + cy_north * np.cos(theta)) * pos_factor
#         else:
#             radec_center = np.average(
#                 catalog["radec_sky"], weights=1.0 / (catalog["radec_sky_err"] ** 2)
#             )
#             dec_center = np.average(
#                 catalog["dec_sky"], weights=1.0 / (catalog["dec_sky_err"] ** 2)
#             )

#         # Extract already-scaled sky coordinates [arcsec] and proper motions [mas/yr]
#         x = catalog["radec_sky"] - radec_center
#         y = catalog["dec_sky"] - dec_center
#         u = catalog["pm_radec_sky"]
#         v = catalog["pm_dec_sky"]
#         u_err = catalog["pm_radec_sky_err"]
#         v_err = catalog["pm_dec_sky_err"]

#     else:
#         # Determine center coordinates in pixel units if not provided
#         if (radec_center is None) or (dec_center is None):
#             radec_center = np.average(
#                 catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
#             )
#             dec_center = np.average(
#                 catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
#             )

#         # Extract coordinates [arcsec] and velocities [mas/yr] with plate scale conversion
#         conversion_factor = plate_scale
#         x = (catalog[radec_col] - radec_center) * conversion_factor / scale_divisor
#         y = (catalog[dec_col] - dec_center) * conversion_factor / scale_divisor

#         u = catalog[pm_radec_col] * conversion_factor
#         v = catalog[pm_dec_col] * conversion_factor
#         u_err = catalog[pm_radec_err_col] * conversion_factor
#         v_err = catalog[pm_dec_err_col] * conversion_factor

#     # Center bulk motion relative to overall weighted mean
#     clust_membership = (
#         catalog["Membership"].data * catalog["Membership_color"].data
#     )
#     weights = clust_membership if weighted else np.ones_like(u)

#     u_centered = u - np.average(u, weights=weights)
#     v_centered = v - np.average(v, weights=weights)

#     # Calculate velocities in units of sigma
#     u_sigma = u_centered / u_err
#     v_sigma = v_centered / v_err

#     # Set up 2D spatial bin grid edges
#     x_bins = np.linspace(x.min(), x.max(), nbins + 1)
#     y_bins = np.linspace(y.min(), y.max(), nbins + 1)

#     # Compute star counts per 2D spatial bin
#     counts, _, _, _ = binned_statistic_2d(
#         x, y, None, statistic="count", bins=[x_bins, y_bins]
#     )

#     # Helper to compute weighted or unweighted binned 2D mean map
#     def compute_binned_map(values):
#         if weighted:
#             w_sum, _, _, _ = binned_statistic_2d(
#                 x, y, values * weights, statistic="sum", bins=[x_bins, y_bins]
#             )
#             weight_sum, _, _, _ = binned_statistic_2d(
#                 x, y, weights, statistic="sum", bins=[x_bins, y_bins]
#             )
#             with np.errstate(divide="ignore", invalid="ignore"):
#                 grid_map = w_sum / weight_sum
#         else:
#             grid_map, _, _, _ = binned_statistic_2d(
#                 x, y, values, statistic="mean", bins=[x_bins, y_bins]
#             )

#         # Mask bins with fewer than min_stars_per_bin
#         grid_map[counts < min_stars_per_bin] = np.nan
#         return grid_map

#     # Generate 2D grid maps for each component
#     grid_u = compute_binned_map(u_centered)
#     grid_v = compute_binned_map(v_centered)
#     grid_u_sigma = compute_binned_map(u_sigma)
#     grid_v_sigma = compute_binned_map(v_sigma)

#     # Make figure
#     fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

#     # Configure panels data and formatting
#     panels = [
#         (
#             axes[0, 0],
#             grid_u,
#             r"$\mu_x$ [mas/yr]",
#             "X Proper Motion [mas/yr]",
#         ),
#         (
#             axes[0, 1],
#             grid_v,
#             r"$\mu_y$ [mas/yr]",
#             "Y Proper Motion [mas/yr]",
#         ),
#         (
#             axes[1, 0],
#             grid_u_sigma,
#             r"$\mu_x / \sigma_x$",
#             r"X Proper Motion [$\sigma$]",
#         ),
#         (
#             axes[1, 1],
#             grid_v_sigma,
#             r"$\mu_y / \sigma_y$",
#             r"Y Proper Motion [$\sigma$]",
#         ),
#     ]

#     for ax, grid_data, cbar_label, title in panels:
#         # Compute symmetric limits centered at 0
#         max_abs_val = np.nanpercentile(np.abs(grid_data), 95)
#         norm = TwoSlopeNorm(
#             vmin=-max_abs_val, vcenter=0.0, vmax=max_abs_val
#         )

#         # Plot 2D density grid map
#         im = ax.pcolormesh(
#             x_bins,
#             y_bins,
#             grid_data.T,  # Transpose to align (x, y) orientation
#             cmap=cmap,
#             norm=norm,
#             shading="flat",
#         )

#         # Plot center marker
#         ax.plot(0, 0, "r+", markersize=10, markeredgewidth=1.5, label="Center")

#         # Colorbar
#         cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
#         cbar.set_label(cbar_label, fontsize=12)
#         cbar.ax.tick_params(labelsize=12)

#         ax.set_title(title, fontsize=16, pad=10)
#         ax.set_aspect("equal")
#         ax.grid(True, linestyle=":", alpha=0.5)

#         # Invert X-axis only when raw X increases to the right;
#         # transform_pixel_to_sky already negates radec_sky so East is on the left
#         if not transform_to_sky:
#             ax.invert_xaxis()

#         ax.tick_params(labelsize=16)

#     # Set outer axis labels
#     for ax in axes[1, :]:
#         ax.set_xlabel("RAcosDEC [arcsec]", fontsize=16)
#     for ax in axes[:, 0]:
#         ax.set_ylabel("DEC [arcsec]", fontsize=16)

#     plt.tight_layout()
#     plt.show()

#     return fig, axes

def plot_pm_density_grid(
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
    radec_center=335,
    dec_center=499,
    cmap="coolwarm",
    nbins=10,
    min_stars_per_bin=3,
    weighted=True,
    transform_to_sky=True,
    pa_deg=-42.3,
    scale_divisor=1000.0,
    pixel_east_left=True,
    to_galactic=False,
    ra_zero_deg=266.56304167,
    dec_zero_deg=-28.82630556,
    filename_out = 'star_pm_density_grid'
):
    """Plots a 2x2 grid of spatial density maps showing binned average X and Y proper motions
    in mas/yr (top row) and in units of sigma (bottom row), with colorbars centered at zero.
    Optionally applies a position angle rotation, pixel-to-sky transformation, and Galactic conversion.
    """
    # Load cluster catalog
    catalog = cru.load_cluster_members(fits_file, clust_mem_min)

    if transform_to_sky:
        catalog = transform_pixel_to_sky(
            catalog,
            x_col=radec_col,
            y_col=dec_col,
            x_err_col=radec_err_col,
            y_err_col=dec_err_col,
            pm_x_col=pm_radec_col,
            pm_y_col=pm_dec_col,
            pm_x_err_col=pm_radec_err_col,
            pm_y_err_col=pm_dec_err_col,
            pa_deg=pa_deg,
            pixel_scale=plate_scale,
            scale_divisor=scale_divisor,
            pixel_east_left=pixel_east_left,
            to_galactic=to_galactic,
            ra_zero_deg=ra_zero_deg,
            dec_zero_deg=dec_zero_deg,
        )

        # Transform provided pixel centers into the same rotated sky/Galactic frame
        if (radec_center is not None) and (dec_center is not None):
            theta = np.deg2rad(pa_deg)
            pos_factor = plate_scale / scale_divisor
            sign_flip = 1.0 if pixel_east_left else -1.0

            cx_east = sign_flip * radec_center
            cy_north = dec_center
            center_x_sky = -(cx_east * np.cos(theta) + cy_north * np.sin(theta)) * pos_factor
            center_y_sky = (-cx_east * np.sin(theta) + cy_north * np.cos(theta)) * pos_factor

            if to_galactic:
                origin = SkyCoord(ra=ra_zero_deg * u.deg, dec=dec_zero_deg * u.deg, frame="icrs")
                center_sky = SkyCoord(
                    lon=-center_x_sky * u.arcsec,
                    lat=center_y_sky * u.arcsec,
                    frame=SkyOffsetFrame(origin=origin),
                )
                gal_origin = origin.transform_to(Galactic())
                center_gal_offset = center_sky.transform_to(Galactic()).transform_to(
                    SkyOffsetFrame(origin=gal_origin)
                )
                radec_center = center_gal_offset.lon.to(u.arcsec).value
                dec_center = center_gal_offset.lat.to(u.arcsec).value
            else:
                radec_center = center_x_sky
                dec_center = center_y_sky
        else:
            if to_galactic:
                radec_center = np.average(
                    catalog["l_gal_offset"], weights=1.0 / (catalog["l_gal_err"] ** 2)
                )
                dec_center = np.average(
                    catalog["b_gal_offset"], weights=1.0 / (catalog["b_gal_err"] ** 2)
                )
            else:
                radec_center = np.average(
                    catalog["radec_sky"], weights=1.0 / (catalog["radec_sky_err"] ** 2)
                )
                dec_center = np.average(
                    catalog["dec_sky"], weights=1.0 / (catalog["dec_sky_err"] ** 2)
                )

        # Extract coordinates and proper motions from the transformed columns
        if to_galactic:
            x = catalog["l_gal_offset"] - radec_center
            y = catalog["b_gal_offset"] - dec_center
            u_pm = catalog["pm_l_gal"]
            v_pm = catalog["pm_b_gal"]
            u_err = catalog["pm_l_gal_err"]
            v_err = catalog["pm_b_gal_err"]
        else:
            x = catalog["radec_sky"] - radec_center
            y = catalog["dec_sky"] - dec_center
            u_pm = catalog["pm_radec_sky"]
            v_pm = catalog["pm_dec_sky"]
            u_err = catalog["pm_radec_sky_err"]
            v_err = catalog["pm_dec_sky_err"]

    else:
        # Determine center coordinates in pixel units if not provided
        if (radec_center is None) or (dec_center is None):
            radec_center = np.average(
                catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
            )
            dec_center = np.average(
                catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
            )

        # Extract coordinates and velocities with plate scale conversion
        conversion_factor = plate_scale
        x = (catalog[radec_col] - radec_center) * conversion_factor / scale_divisor
        y = (catalog[dec_col] - dec_center) * conversion_factor / scale_divisor

        u_pm = catalog[pm_radec_col] * conversion_factor
        v_pm = catalog[pm_dec_col] * conversion_factor
        u_err = catalog[pm_radec_err_col] * conversion_factor
        v_err = catalog[pm_dec_err_col] * conversion_factor

    # Center bulk motion relative to overall weighted mean
    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )
    weights = clust_membership if weighted else np.ones_like(u_pm)

    u_centered = u_pm - np.average(u_pm, weights=weights)
    v_centered = v_pm - np.average(v_pm, weights=weights)

    # Calculate velocities in units of sigma
    u_sigma = u_centered / u_err
    v_sigma = v_centered / v_err

    # Set up 2D spatial bin grid edges
    x_bins = np.linspace(x.min(), x.max(), nbins + 1)
    y_bins = np.linspace(y.min(), y.max(), nbins + 1)

    # Compute star counts per 2D spatial bin
    counts, _, _, _ = binned_statistic_2d(
        x, y, None, statistic="count", bins=[x_bins, y_bins]
    )

    # Helper to compute weighted or unweighted binned 2D mean map
    def compute_binned_map(values):
        if weighted:
            w_sum, _, _, _ = binned_statistic_2d(
                x, y, values * weights, statistic="sum", bins=[x_bins, y_bins]
            )
            weight_sum, _, _, _ = binned_statistic_2d(
                x, y, weights, statistic="sum", bins=[x_bins, y_bins]
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                grid_map = w_sum / weight_sum
        else:
            grid_map, _, _, _ = binned_statistic_2d(
                x, y, values, statistic="mean", bins=[x_bins, y_bins]
            )

        # Mask bins with fewer than min_stars_per_bin
        grid_map[counts < min_stars_per_bin] = np.nan
        return grid_map

    # Generate 2D grid maps for each component
    grid_u = compute_binned_map(u_centered)
    grid_v = compute_binned_map(v_centered)
    grid_u_sigma = compute_binned_map(u_sigma)
    grid_v_sigma = compute_binned_map(v_sigma)

    # Make figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

    # Set frame-appropriate labels
    if transform_to_sky and to_galactic:
        panels = [
            (axes[0, 0], grid_u, r"$\mu_{l*}$ [mas/yr]", "Galactic Longitude PM [mas/yr]"),
            (axes[0, 1], grid_v, r"$\mu_b$ [mas/yr]", "Galactic Latitude PM [mas/yr]"),
            (axes[1, 0], grid_u_sigma, r"$\mu_{l*} / \sigma_{l*}$", r"Galactic Longitude PM [$\sigma$]"),
            (axes[1, 1], grid_v_sigma, r"$\mu_b / \sigma_b$", r"Galactic Latitude PM [$\sigma$]"),
        ]
        xlabel, ylabel = r"$\Delta l \cos b$ [arcsec]", r"$\Delta b$ [arcsec]"
    elif transform_to_sky and not to_galactic:
        panels = [
            (axes[0, 0], grid_u, r"$\mu_{\alpha^{*}}$ [mas/yr]", "RAcosDEC Proper Motion [mas/yr]"),
            (axes[0, 1], grid_v, r"$\mu_{\delta}$ [mas/yr]", "DEC Proper Motion [mas/yr]"),
            (axes[1, 0], grid_u_sigma, r"$\mu_{\alpha^{*}} / \sigma_{\alpha^{*}}$", r"RAcosDEC Proper Motion [$\sigma$]"),
            (axes[1, 1], grid_v_sigma, r"$\mu_{\delta} / \sigma_{\delta}$", r"DEC Proper Motion [$\sigma$]"),
        ]
        unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
        xlabel, ylabel = f"RAcosDEC [{unit_str}]", f"DEC [{unit_str}]"
    else:
        panels = [
            (axes[0, 0], grid_u, r"$\mu_x$ [mas/yr]", "X Proper Motion [mas/yr]"),
            (axes[0, 1], grid_v, r"$\mu_y$ [mas/yr]", "Y Proper Motion [mas/yr]"),
            (axes[1, 0], grid_u_sigma, r"$\mu_x / \sigma_x$", r"X Proper Motion [$\sigma$]"),
            (axes[1, 1], grid_v_sigma, r"$\mu_y / \sigma_y$", r"Y Proper Motion [$\sigma$]"),
        ]
        unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
        xlabel, ylabel = f"X [{unit_str}]", f"Y [{unit_str}]"

    for ax, grid_data, cbar_label, title in panels:
        # Compute symmetric limits centered at 0
        max_abs_val = np.nanpercentile(np.abs(grid_data), 95)
        norm = TwoSlopeNorm(
            vmin=-max_abs_val, vcenter=0.0, vmax=max_abs_val
        )

        # Plot 2D density grid map
        im = ax.pcolormesh(
            x_bins,
            y_bins,
            grid_data.T,  # Transpose to align (x, y) orientation
            cmap=cmap,
            norm=norm,
            shading="flat",
        )

        # Plot center marker
        ax.plot(0, 0, "r+", markersize=10, markeredgewidth=1.5, label="Center")

        # Colorbar
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, fontsize=12)
        cbar.ax.tick_params(labelsize=12)

        ax.set_title(title, fontsize=16, pad=10)
        ax.set_aspect("equal")
        ax.grid(True, linestyle=":", alpha=0.5)

        # Invert X-axis when:
        # 1. transform_to_sky=False (raw detector X increases right)
        # 2. to_galactic=True (l_gal_offset increases positively toward +l, which sits to the left on sky)
        if (not transform_to_sky) or to_galactic:
            ax.invert_xaxis()

        ax.tick_params(labelsize=16)

    # Set outer axis labels
    for ax in axes[1, :]:
        ax.set_xlabel(xlabel, fontsize=16)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=16)

    plt.tight_layout()
    #plt.show()

    fig.savefig(filename_out + '.png', dpi = 300)

    return fig, axes


def plot_individual_stars_pm_grid(
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
    radec_center=335,
    dec_center=499,
    cmap="coolwarm",
    markersize=15,
    alpha=1.0,
    weighted=True,
    transform_to_sky=True,
    pa_deg=-42.3,
    scale_divisor=1000.0,
    pixel_east_left=True,
    to_galactic=False,
    ra_zero_deg=266.56304167,
    dec_zero_deg=-28.82630556,
    filename_out="star_pm_scatter_grid",
):
    """Plots individual cluster member stars in a 2x2 grid of spatial scatter plots

    colored by X and Y proper motions in mas/yr (top row) and in units of sigma
    (bottom row),
    with colorbars centered at zero. Optionally applies a position angle
    rotation,
    pixel-to-sky transformation, and Galactic conversion.
    """
    # Load cluster catalog
    catalog = cru.load_cluster_members(fits_file, clust_mem_min)

    if transform_to_sky:
        catalog = transform_pixel_to_sky(
            catalog,
            x_col=radec_col,
            y_col=dec_col,
            x_err_col=radec_err_col,
            y_err_col=dec_err_col,
            pm_x_col=pm_radec_col,
            pm_y_col=pm_dec_col,
            pm_x_err_col=pm_radec_err_col,
            pm_y_err_col=pm_dec_err_col,
            pa_deg=pa_deg,
            pixel_scale=plate_scale,
            scale_divisor=scale_divisor,
            pixel_east_left=pixel_east_left,
            to_galactic=to_galactic,
            ra_zero_deg=ra_zero_deg,
            dec_zero_deg=dec_zero_deg,
        )

        # Transform provided pixel centers into the same rotated sky/Galactic frame
        if (radec_center is not None) and (dec_center is not None):
            theta = np.deg2rad(pa_deg)
            pos_factor = plate_scale / scale_divisor
            sign_flip = 1.0 if pixel_east_left else -1.0

            cx_east = sign_flip * radec_center
            cy_north = dec_center
            center_x_sky = (
                -(cx_east * np.cos(theta) + cy_north * np.sin(theta))
                * pos_factor
            )
            center_y_sky = (
                -cx_east * np.sin(theta) + cy_north * np.cos(theta)
            ) * pos_factor

            if to_galactic:
                origin = SkyCoord(
                    ra=ra_zero_deg * u.deg, dec=dec_zero_deg * u.deg, frame="icrs"
                )
                center_sky = SkyCoord(
                    lon=-center_x_sky * u.arcsec,
                    lat=center_y_sky * u.arcsec,
                    frame=SkyOffsetFrame(origin=origin),
                )
                gal_origin = origin.transform_to(Galactic())
                center_gal_offset = center_sky.transform_to(
                    Galactic()
                ).transform_to(SkyOffsetFrame(origin=gal_origin))
                radec_center = center_gal_offset.lon.to(u.arcsec).value
                dec_center = center_gal_offset.lat.to(u.arcsec).value
            else:
                radec_center = center_x_sky
                dec_center = center_y_sky
        else:
            if to_galactic:
                radec_center = np.average(
                    catalog["l_gal_offset"],
                    weights=1.0 / (catalog["l_gal_err"] ** 2),
                )
                dec_center = np.average(
                    catalog["b_gal_offset"],
                    weights=1.0 / (catalog["b_gal_err"] ** 2),
                )
            else:
                radec_center = np.average(
                    catalog["radec_sky"],
                    weights=1.0 / (catalog["radec_sky_err"] ** 2),
                )
                dec_center = np.average(
                    catalog["dec_sky"],
                    weights=1.0 / (catalog["dec_sky_err"] ** 2),
                )

        # Extract coordinates and proper motions from the transformed columns
        if to_galactic:
            x = catalog["l_gal_offset"] - radec_center
            y = catalog["b_gal_offset"] - dec_center
            u_pm = catalog["pm_l_gal"]
            v_pm = catalog["pm_b_gal"]
            u_err = catalog["pm_l_gal_err"]
            v_err = catalog["pm_b_gal_err"]
        else:
            x = catalog["radec_sky"] - radec_center
            y = catalog["dec_sky"] - dec_center
            u_pm = catalog["pm_radec_sky"]
            v_pm = catalog["pm_dec_sky"]
            u_err = catalog["pm_radec_sky_err"]
            v_err = catalog["pm_dec_sky_err"]

    else:
        # Determine center coordinates in pixel units if not provided
        if (radec_center is None) or (dec_center is None):
            radec_center = np.average(
                catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
            )
            dec_center = np.average(
                catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
            )

        # Extract coordinates and velocities with plate scale conversion
        conversion_factor = plate_scale
        x = (
            (catalog[radec_col] - radec_center)
            * conversion_factor
            / scale_divisor
        )
        y = (catalog[dec_col] - dec_center) * conversion_factor / scale_divisor

        u_pm = catalog[pm_radec_col] * conversion_factor
        v_pm = catalog[pm_dec_col] * conversion_factor
        u_err = catalog[pm_radec_err_col] * conversion_factor
        v_err = catalog[pm_dec_err_col] * conversion_factor

    # Center bulk motion relative to overall weighted mean
    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )
    weights = clust_membership if weighted else np.ones_like(u_pm)

    u_centered = u_pm - np.average(u_pm, weights=weights)
    v_centered = v_pm - np.average(v_pm, weights=weights)

    # Calculate velocities in units of sigma
    u_sigma = u_centered / u_err
    v_sigma = v_centered / v_err

    # Make figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

    # Set frame-appropriate labels
    if transform_to_sky and to_galactic:
        panels = [
            (
                axes[0, 0],
                u_centered,
                r"$\mu_{l*}$ [mas/yr]",
                "Galactic Longitude PM [mas/yr]",
            ),
            (
                axes[0, 1],
                v_centered,
                r"$\mu_b$ [mas/yr]",
                "Galactic Latitude PM [mas/yr]",
            ),
            (
                axes[1, 0],
                u_sigma,
                r"$\mu_{l*} / \sigma_{l*}$",
                r"Galactic Longitude PM [$\sigma$]",
            ),
            (
                axes[1, 1],
                v_sigma,
                r"$\mu_b / \sigma_b$",
                r"Galactic Latitude PM [$\sigma$]",
            ),
        ]
        xlabel, ylabel = r"$\Delta l \cos b$ [arcsec]", r"$\Delta b$ [arcsec]"
    elif transform_to_sky and not to_galactic:
        panels = [
            (
                axes[0, 0],
                u_centered,
                r"$\mu_{\alpha^{*}}$ [mas/yr]",
                "RAcosDEC Proper Motion [mas/yr]",
            ),
            (
                axes[0, 1],
                v_centered,
                r"$\mu_{\delta}$ [mas/yr]",
                "DEC Proper Motion [mas/yr]",
            ),
            (
                axes[1, 0],
                u_sigma,
                r"$\mu_{\alpha^{*}} / \sigma_{\alpha^{*}}$",
                r"RAcosDEC Proper Motion [$\sigma$]",
            ),
            (
                axes[1, 1],
                v_sigma,
                r"$\mu_{\delta} / \sigma_{\delta}$",
                r"DEC Proper Motion [$\sigma$]",
            ),
        ]
        unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
        xlabel, ylabel = f"RAcosDEC [{unit_str}]", f"DEC [{unit_str}]"
    else:
        panels = [
            (axes[0, 0], u_centered, r"$\mu_x$ [mas/yr]", "X Proper Motion [mas/yr]"),
            (axes[0, 1], v_centered, r"$\mu_y$ [mas/yr]", "Y Proper Motion [mas/yr]"),
            (
                axes[1, 0],
                u_sigma,
                r"$\mu_x / \sigma_x$",
                r"X Proper Motion [$\sigma$]",
            ),
            (
                axes[1, 1],
                v_sigma,
                r"$\mu_y / \sigma_y$",
                r"Y Proper Motion [$\sigma$]",
            ),
        ]
        unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
        xlabel, ylabel = f"X [{unit_str}]", f"Y [{unit_str}]"

    for ax, color_data, cbar_label, title in panels:
        # Compute symmetric limits centered at 0
        max_abs_val = np.nanpercentile(np.abs(color_data), 95)
        norm = TwoSlopeNorm(
            vmin=-max_abs_val, vcenter=0.0, vmax=max_abs_val
        )

        # Scatter plot for individual stars (x and y are already scaled by scale_divisor)
        sc = ax.scatter(
            x,
            y,
            c=color_data,
            cmap=cmap,
            norm=norm,
            s=markersize,
            alpha=alpha,
            marker="o",
            edgecolors="none",
        )

        # Plot center marker
        ax.plot(0, 0, "r+", markersize=10, markeredgewidth=1.5, label="Center")

        # Colorbar
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, fontsize=12)
        cbar.ax.tick_params(labelsize=12)

        ax.set_title(title, fontsize=16, pad=10)
        ax.set_aspect("equal")
        ax.grid(True, linestyle=":", alpha=0.5)

        # Invert X-axis when:
        # 1. transform_to_sky=False (raw detector X increases right)
        # 2. to_galactic=True (l_gal_offset increases positively toward +l, which sits to the left on sky)
        if (not transform_to_sky) or to_galactic:
            ax.invert_xaxis()

        ax.tick_params(labelsize=16)

    # Set outer axis labels
    for ax in axes[1, :]:
        ax.set_xlabel(xlabel, fontsize=16)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=16)

    plt.tight_layout()
    fig.savefig(filename_out + ".png", dpi=300)

    return fig, axes