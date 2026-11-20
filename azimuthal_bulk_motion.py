# import matplotlib.patches as mpatches
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.stats import binned_statistic_2d
# from cluster_rotation_analysis import scatter_proper_motion as spm
# from cluster_rotation_analysis import clust_rot_utils as cru


# def plot_azimuthal_star_velocities_sigma(
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
#     scale=None,
#     key_length=1.0,
#     key_label=r"$1\sigma$",
#     nbins_r=5,
#     nbins_theta=12,
#     r_max=None,
#     min_stars_per_bin=3,
#     weighted=True,
#     error_type="std",
#     plot_type="tangential",
#     cmap="viridis",
#     highlight_threshold=1.0,  
#     highlight_color="yellow",  
#     highlight_alpha=0.3, 
# ):
#     """Grids spatial coordinates in polar/azimuthal wedges, plots binned average
#     star proper motions expressed in units of sigma, and highlights sectors
#     where the vector magnitude exceeds a specified sigma threshold.
#     """
#     # Load filtered catalog using existing function
#     catalog = cru.load_cluster_members(fits_file, clust_mem_min)

#     if (radec_center is None) or (dec_center is None):
#         radec_center = np.average(
#             catalog[radec_col], weights=1 / (catalog[radec_err_col] ** 2)
#         )
#         dec_center = np.average(
#             catalog[dec_col], weights=1 / (catalog[dec_err_col] ** 2)
#         )

#     # Extract coordinates and velocities, applying the plate_scale factor
#     conversion_factor = plate_scale / 1000.0
#     x = (catalog[radec_col] - radec_center) * conversion_factor
#     y = (catalog[dec_col] - dec_center) * conversion_factor
#     u = catalog[pm_radec_col] * conversion_factor
#     v = catalog[pm_dec_col] * conversion_factor

#     # Extract errors (scaled appropriately)
#     u_err = catalog[pm_radec_err_col] * conversion_factor
#     v_err = catalog[pm_dec_err_col] * conversion_factor

#     clust_membership = (
#         catalog["Membership"].data * catalog["Membership_color"].data
#     )

#     # Center bulk motion relative to overall weighted mean
#     weights_overall = clust_membership if weighted else np.ones_like(u)
#     u_centered = u - np.average(u, weights=weights_overall)
#     v_centered = v - np.average(v, weights=weights_overall)

#     # Convert Cartesian (x, y) to Polar (r, theta) coordinates
#     r = np.hypot(x, y)
#     theta = np.arctan2(y, x)  # Angles in radians from -pi to +pi

#     # If requested, convert centered velocities to tangential velocity vectors
#     if plot_type == "tangential":
#         v_tan = -u_centered * np.sin(theta) + v_centered * np.cos(theta)
#         v_tan_err = np.sqrt(
#             (np.sin(theta) * u_err) ** 2 + (np.cos(theta) * v_err) ** 2
#         )

#         u_centered = -v_tan * np.sin(theta)
#         v_centered = v_tan * np.cos(theta)
#         u_err = np.abs(v_tan_err * np.sin(theta))
#         v_err = np.abs(v_tan_err * np.cos(theta))

#     # Define polar bin edges
#     r_upper = r_max if r_max is not None else r.max()
#     r_bins = np.linspace(0, r_upper, nbins_r + 1)
#     theta_bins = np.linspace(-np.pi, np.pi, nbins_theta + 1)

#     # Perform 2D binned statistics in (r, theta) space
#     counts, _, _, _ = binned_statistic_2d(
#         r, theta, None, statistic="count", bins=[r_bins, theta_bins]
#     )

#     if weighted:
#         u_weighted_sum, _, _, _ = binned_statistic_2d(
#             r,
#             theta,
#             u_centered * weights_overall,
#             statistic="sum",
#             bins=[r_bins, theta_bins],
#         )
#         v_weighted_sum, _, _, _ = binned_statistic_2d(
#             r,
#             theta,
#             v_centered * weights_overall,
#             statistic="sum",
#             bins=[r_bins, theta_bins],
#         )
#         weight_sum, _, _, _ = binned_statistic_2d(
#             r,
#             theta,
#             weights_overall,
#             statistic="sum",
#             bins=[r_bins, theta_bins],
#         )

#         with np.errstate(divide="ignore", invalid="ignore"):
#             grid_u_mean = u_weighted_sum / weight_sum
#             grid_v_mean = v_weighted_sum / weight_sum
#             print(grid_u_mean, grid_v_mean)
#             if error_type == "std":
#                 u_sq_w_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     (u_centered**2) * weights_overall,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 v_sq_w_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     (v_centered**2) * weights_overall,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 u_var = (u_sq_w_sum / weight_sum) - (grid_u_mean**2)
#                 v_var = (v_sq_w_sum / weight_sum) - (grid_v_mean**2)

#                 grid_u_sigma_err = np.sqrt(np.maximum(0.0, u_var))
#                 grid_v_sigma_err = np.sqrt(np.maximum(0.0, v_var))
#             else:
#                 u_err_sq_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     (weights_overall * u_err) ** 2,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 v_err_sq_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     (weights_overall * v_err) ** 2,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 grid_u_sigma_err = np.sqrt(u_err_sq_sum) / weight_sum
#                 grid_v_sigma_err = np.sqrt(v_err_sq_sum) / weight_sum
#     else:
#         grid_u_mean, _, _, _ = binned_statistic_2d(
#             r, theta, u_centered, statistic="mean", bins=[r_bins, theta_bins]
#         )
#         grid_v_mean, _, _, _ = binned_statistic_2d(
#             r, theta, v_centered, statistic="mean", bins=[r_bins, theta_bins]
#         )

#         with np.errstate(divide="ignore", invalid="ignore"):
#             if error_type == "std":
#                 grid_u_sigma_err, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     u_centered,
#                     statistic="std",
#                     bins=[r_bins, theta_bins],
#                 )
#                 grid_v_sigma_err, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     v_centered,
#                     statistic="std",
#                     bins=[r_bins, theta_bins],
#                 )
#             else:
#                 u_err_sq_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     u_err**2,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 v_err_sq_sum, _, _, _ = binned_statistic_2d(
#                     r,
#                     theta,
#                     v_err**2,
#                     statistic="sum",
#                     bins=[r_bins, theta_bins],
#                 )
#                 grid_u_sigma_err = np.sqrt(u_err_sq_sum) / counts
#                 grid_v_sigma_err = np.sqrt(v_err_sq_sum) / counts

#     # Convert velocity components to units of sigma
#     with np.errstate(divide="ignore", invalid="ignore"):
#         grid_u = grid_u_mean / grid_u_sigma_err
#         grid_v = grid_v_mean / grid_v_sigma_err

#     # Compute total magnitude in units of sigma for each 2D bin
#     grid_sigma_mag = np.hypot(grid_u, grid_v)

#     # Compute midpoints for each wedge in (r, theta) to place the vectors
#     r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
#     theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
#     grid_r, grid_theta = np.meshgrid(r_centers, theta_centers, indexing="ij")

#     # Convert wedge centerpoints back to Cartesian coordinates for plotting
#     grid_x = grid_r * np.cos(grid_theta)
#     grid_y = grid_r * np.sin(grid_theta)

#     # Mask low-density bins and invalid divisions
#     valid_mask = (
#         (counts >= min_stars_per_bin)
#         & np.isfinite(grid_u)
#         & np.isfinite(grid_v)
#     )
#     plot_x = grid_x[valid_mask]
#     plot_y = grid_y[valid_mask]
#     plot_u = grid_u[valid_mask]
#     plot_v = grid_v[valid_mask]
#     plot_counts = counts[valid_mask]

#     # Initialize plot
#     fig, ax = plt.subplots(figsize=(8.5, 8))

#     # Highlight sectors where vector significance exceeds highlight_threshold (e.g., > 1 sigma)
#     highlight_added = False
#     for i in range(nbins_r):
#         for j in range(nbins_theta):
#             if valid_mask[i, j] and grid_sigma_mag[i, j] > highlight_threshold:
#                 r1, r2 = r_bins[i], r_bins[i + 1]
#                 t1 = np.degrees(theta_bins[j])
#                 t2 = np.degrees(theta_bins[j + 1])

#                 wedge = mpatches.Wedge(
#                     center=(0, 0),
#                     r=r2,
#                     theta1=t1,
#                     theta2=t2,
#                     width=(r2 - r1),
#                     color=highlight_color,
#                     alpha=highlight_alpha,
#                     zorder=1,
#                 )
#                 ax.add_patch(wedge)
#                 highlight_added = True

#     # Plot azimuthal grid outlines
#     for r_edge in r_bins:
#         circle = plt.Circle(
#             (0, 0),
#             r_edge,
#             color="gray",
#             linestyle="--",
#             fill=False,
#             alpha=0.5,
#             zorder=2,
#         )
#         ax.add_patch(circle)

#     for t_edge in theta_bins:
#         ax.plot(
#             [0, r_upper * np.cos(t_edge)],
#             [0, r_upper * np.sin(t_edge)],
#             color="gray",
#             linestyle="--",
#             alpha=0.5,
#             zorder=2,
#         )

#     # Plot velocity vectors color-coded by number of stars (counts)
#     quiver = ax.quiver(
#         plot_x,
#         plot_y,
#         plot_u,
#         plot_v,
#         plot_counts,
#         cmap=cmap,
#         angles="xy",
#         scale_units="xy",
#         scale=scale,
#         alpha=0.9,
#         width=0.005,
#         zorder=3,
#     )

#     # Add Colorbar for Star Counts
#     cbar = fig.colorbar(quiver, ax=ax, fraction=0.046, pad=0.04)
#     cbar.set_label("Number of stars per bin ($N$)", fontsize=14)
#     cbar.ax.tick_params(labelsize=12)

#     # Add Quiver Key
#     ax.quiverkey(
#         quiver,
#         X=0.9,
#         Y=0.9,
#         U=key_length,
#         label=key_label,
#         labelpos="N",
#         coordinates="axes",
#         color="r",
#     )

#     # Center marker
#     ax.plot(0, 0, "r+", markersize=10, label="Center", zorder=4)

#     # Add legend entry for highlighted sectors if any are drawn
#     if highlight_added:
#         highlight_patch = mpatches.Patch(
#             color=highlight_color,
#             alpha=highlight_alpha,
#             label=f"$>{highlight_threshold:.0f}\\sigma$ Motion",
#         )
#         handles, labels = ax.get_legend_handles_labels()
#         handles.append(highlight_patch)
#         ax.legend(handles=handles, loc="upper left", fontsize=11)

#     # Set boundaries and formatting
#     ax.set_aspect("equal")
#     ax.set_xlim(-r_upper * 1.05, r_upper * 1.05)
#     ax.set_ylim(-r_upper * 1.05, r_upper * 1.05)

#     ax.set_xlabel("RAcosDEC [arcsec]", fontsize=16)
#     ax.set_ylabel("DEC [arcsec]", fontsize=16)
#     ax.tick_params(labelsize=14)

#     # Invert X-axis for sky projection alignment
#     ax.invert_xaxis()

#     plt.tight_layout()
#     plt.show()
#     fig.savefig("azimuthal_binned_motion_sigma.png", dpi=300)


import astropy.units as u
from astropy.coordinates import Galactic, SkyCoord, SkyOffsetFrame
from cluster_rotation_analysis import clust_rot_utils as cru
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from cluster_rotation_analysis.scatter_proper_motion import transform_pixel_to_sky as spm
from scipy.stats import binned_statistic_2d


def plot_azimuthal_star_velocities_sigma(
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
    scale=None,
    key_length=1.0,
    key_label=r"$1\sigma$",
    nbins_r=5,
    nbins_theta=12,
    r_max=None,
    min_stars_per_bin=3,
    weighted=True,
    error_type="std",
    plot_type="tangential",
    cmap="viridis",
    highlight_threshold=1.0,
    highlight_color="yellow",
    highlight_alpha=0.3,
    transform_to_sky=True,
    pa_deg=-42.3,
    scale_divisor=1000.0,
    pixel_east_left=True,
    to_galactic=False,
    ra_zero_deg=266.56304167,
    dec_zero_deg=-28.82630556,
    filename_out="azimuthal_binned_motion_sigma",
):
    """Grids spatial coordinates in polar/azimuthal wedges, plots binned average

    star proper motions expressed in units of sigma, and highlights sectors
    where the vector magnitude exceeds a specified sigma threshold.
    Optionally applies a position angle rotation, pixel-to-sky transformation,
    and Galactic conversion.
    """
    # Load filtered catalog using existing function
    catalog = cru.load_cluster_members(fits_file, clust_mem_min)

    if transform_to_sky:
        # Note: if spm is imported as the module (e.g., `import scatter_proper_motion as spm`),
        # call `spm.transform_pixel_to_sky(...)` instead.
        catalog = spm(
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
            theta_pa = np.deg2rad(pa_deg)
            pos_factor = plate_scale / scale_divisor
            sign_flip = 1.0 if pixel_east_left else -1.0

            cx_east = sign_flip * radec_center
            cy_north = dec_center
            center_x_sky = (
                -(cx_east * np.cos(theta_pa) + cy_north * np.sin(theta_pa))
                * pos_factor
            )
            center_y_sky = (
                -cx_east * np.sin(theta_pa) + cy_north * np.cos(theta_pa)
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
        if (radec_center is None) or (dec_center is None):
            radec_center = np.average(
                catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
            )
            dec_center = np.average(
                catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
            )

        # Extract coordinates and velocities, applying the plate_scale factor
        x = (catalog[radec_col] - radec_center) * plate_scale / scale_divisor
        y = (catalog[dec_col] - dec_center) * plate_scale / scale_divisor
        u_pm = catalog[pm_radec_col] * plate_scale
        v_pm = catalog[pm_dec_col] * plate_scale
        u_err = catalog[pm_radec_err_col] * plate_scale
        v_err = catalog[pm_dec_err_col] * plate_scale

    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )

    # Center bulk motion relative to overall weighted mean
    weights_overall = clust_membership if weighted else np.ones_like(u_pm)
    u_centered = u_pm - np.average(u_pm, weights=weights_overall)
    v_centered = v_pm - np.average(v_pm, weights=weights_overall)

    # Convert Cartesian (x, y) to Polar (r, theta) coordinates
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)  # Angles in radians from -pi to +pi

    # If requested, convert centered velocities to tangential velocity vectors
    if plot_type == "tangential":
        v_tan = -u_centered * np.sin(theta) + v_centered * np.cos(theta)
        v_tan_err = np.sqrt(
            (np.sin(theta) * u_err) ** 2 + (np.cos(theta) * v_err) ** 2
        )

        u_centered = -v_tan * np.sin(theta)
        v_centered = v_tan * np.cos(theta)
        u_err = np.abs(v_tan_err * np.sin(theta))
        v_err = np.abs(v_tan_err * np.cos(theta))

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
            r,
            theta,
            weights_overall,
            statistic="sum",
            bins=[r_bins, theta_bins],
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            grid_u_mean = u_weighted_sum / weight_sum
            grid_v_mean = v_weighted_sum / weight_sum

            if error_type == "std":
                u_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (u_centered**2) * weights_overall,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (v_centered**2) * weights_overall,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                u_var = (u_sq_w_sum / weight_sum) - (grid_u_mean**2)
                v_var = (v_sq_w_sum / weight_sum) - (grid_v_mean**2)

                grid_u_sigma_err = np.sqrt(np.maximum(0.0, u_var))
                grid_v_sigma_err = np.sqrt(np.maximum(0.0, v_var))
            else:
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (weights_overall * u_err) ** 2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (weights_overall * v_err) ** 2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
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
                grid_u_sigma_err, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    u_centered,
                    statistic="std",
                    bins=[r_bins, theta_bins],
                )
                grid_v_sigma_err, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    v_centered,
                    statistic="std",
                    bins=[r_bins, theta_bins],
                )
            else:
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    u_err**2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    v_err**2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                grid_u_sigma_err = np.sqrt(u_err_sq_sum) / counts
                grid_v_sigma_err = np.sqrt(v_err_sq_sum) / counts

    # Convert velocity components to units of sigma
    with np.errstate(divide="ignore", invalid="ignore"):
        grid_u = grid_u_mean / grid_u_sigma_err
        grid_v = grid_v_mean / grid_v_sigma_err

    # Compute total magnitude in units of sigma for each 2D bin
    grid_sigma_mag = np.hypot(grid_u, grid_v)

    # Compute midpoints for each wedge in (r, theta) to place the vectors
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
    grid_r, grid_theta = np.meshgrid(r_centers, theta_centers, indexing="ij")

    # Convert wedge centerpoints back to Cartesian coordinates for plotting
    grid_x = grid_r * np.cos(grid_theta)
    grid_y = grid_r * np.sin(grid_theta)

    # Mask low-density bins and invalid divisions
    valid_mask = (
        (counts >= min_stars_per_bin)
        & np.isfinite(grid_u)
        & np.isfinite(grid_v)
    )
    plot_x = grid_x[valid_mask]
    plot_y = grid_y[valid_mask]
    plot_u = grid_u[valid_mask]
    plot_v = grid_v[valid_mask]
    plot_counts = counts[valid_mask]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8.5, 8))

    # Highlight sectors where vector significance exceeds highlight_threshold (e.g., > 1 sigma)
    highlight_added = False
    for i in range(nbins_r):
        for j in range(nbins_theta):
            if valid_mask[i, j] and grid_sigma_mag[i, j] > highlight_threshold:
                r1, r2 = r_bins[i], r_bins[i + 1]
                t1 = np.degrees(theta_bins[j])
                t2 = np.degrees(theta_bins[j + 1])

                wedge = mpatches.Wedge(
                    center=(0, 0),
                    r=r2,
                    theta1=t1,
                    theta2=t2,
                    width=(r2 - r1),
                    color=highlight_color,
                    alpha=highlight_alpha,
                    zorder=1,
                )
                ax.add_patch(wedge)
                highlight_added = True

    # Plot azimuthal grid outlines
    for r_edge in r_bins:
        circle = plt.Circle(
            (0, 0),
            r_edge,
            color="gray",
            linestyle="--",
            fill=False,
            alpha=0.5,
            zorder=2,
        )
        ax.add_patch(circle)

    for t_edge in theta_bins:
        ax.plot(
            [0, r_upper * np.cos(t_edge)],
            [0, r_upper * np.sin(t_edge)],
            color="gray",
            linestyle="--",
            alpha=0.5,
            zorder=2,
        )

    # Plot velocity vectors color-coded by number of stars (counts)
    quiver = ax.quiver(
        plot_x,
        plot_y,
        plot_u,
        plot_v,
        plot_counts,
        cmap=cmap,
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.9,
        width=0.005,
        zorder=3,
    )

    # Add Colorbar for Star Counts
    cbar = fig.colorbar(quiver, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Number of stars per bin ($N$)", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    # Add Quiver Key
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
    ax.plot(0, 0, "r+", markersize=10, label="Center", zorder=4)

    # Add legend entry for highlighted sectors if any are drawn
    if highlight_added:
        highlight_patch = mpatches.Patch(
            color=highlight_color,
            alpha=highlight_alpha,
            label=f"$>{highlight_threshold:.0f}\\sigma$ Motion",
        )
        handles, labels = ax.get_legend_handles_labels()
        handles.append(highlight_patch)
        ax.legend(handles=handles, loc="upper left", fontsize=11)

    # Set boundaries and formatting
    ax.set_aspect("equal")
    ax.set_xlim(-r_upper * 1.05, r_upper * 1.05)
    ax.set_ylim(-r_upper * 1.05, r_upper * 1.05)

    # Set frame-appropriate axis labels
    unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
    if transform_to_sky and to_galactic:
        xlabel, ylabel = r"$\Delta l \cos b$ [arcsec]", r"$\Delta b$ [arcsec]"
    elif transform_to_sky and not to_galactic:
        xlabel, ylabel = f"RAcosDEC [{unit_str}]", f"DEC [{unit_str}]"
    else:
        xlabel, ylabel = f"X [{unit_str}]", f"Y [{unit_str}]"

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(labelsize=14)

    # Invert X-axis when:
    # 1. transform_to_sky=False (raw detector X increases right)
    # 2. to_galactic=True (l_gal_offset increases positively toward +l, which sits to the left on sky)
    if (not transform_to_sky) or to_galactic:
        ax.invert_xaxis()

    plt.tight_layout()
    fig.savefig(filename_out + ".png", dpi=300)

    return fig, ax

import os

def plot_azimuthal_star_velocities_sigma_bins(
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
    scale=None,
    key_length=1.0,
    key_label=r"$1\sigma$",
    nbins_r=5,
    nbins_theta=12,
    r_max=None,
    min_stars_per_bin=3,
    weighted=True,
    error_type="std",
    plot_type="tangential",
    cmap="viridis",
    highlight_threshold=1.0,
    highlight_color="yellow",
    highlight_alpha=0.3,
    transform_to_sky=True,
    pa_deg=-42.3,
    scale_divisor=1000.0,
    pixel_east_left=True,
    to_galactic=False,
    ra_zero_deg=266.56304167,
    dec_zero_deg=-28.82630556,
    filename_out="azimuthal_binned_motion_sigma",
    save_bin_histograms=True,
    hist_dir="bin_histograms",
):
    """Grids spatial coordinates in polar/azimuthal wedges, plots binned average

    star proper motions expressed in units of sigma, and highlights sectors
    where the vector magnitude exceeds a specified sigma threshold.
    Optionally applies a position angle rotation, pixel-to-sky transformation,
    Galactic conversion, and saves per-bin velocity distribution histograms.
    """
    # Load filtered catalog using existing function
    catalog = cru.load_cluster_members(fits_file, clust_mem_min)

    if transform_to_sky:
        catalog = spm(
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
            theta_pa = np.deg2rad(pa_deg)
            pos_factor = plate_scale / scale_divisor
            sign_flip = 1.0 if pixel_east_left else -1.0

            cx_east = sign_flip * radec_center
            cy_north = dec_center
            center_x_sky = (
                -(cx_east * np.cos(theta_pa) + cy_north * np.sin(theta_pa))
                * pos_factor
            )
            center_y_sky = (
                -cx_east * np.sin(theta_pa) + cy_north * np.cos(theta_pa)
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
            x = np.asarray(catalog["l_gal_offset"] - radec_center)
            y = np.asarray(catalog["b_gal_offset"] - dec_center)
            u_pm = np.asarray(catalog["pm_l_gal"])
            v_pm = np.asarray(catalog["pm_b_gal"])
            u_err = np.asarray(catalog["pm_l_gal_err"])
            v_err = np.asarray(catalog["pm_b_gal_err"])
        else:
            x = np.asarray(catalog["radec_sky"] - radec_center)
            y = np.asarray(catalog["dec_sky"] - dec_center)
            u_pm = np.asarray(catalog["pm_radec_sky"])
            v_pm = np.asarray(catalog["pm_dec_sky"])
            u_err = np.asarray(catalog["pm_radec_sky_err"])
            v_err = np.asarray(catalog["pm_dec_sky_err"])

    else:
        if (radec_center is None) or (dec_center is None):
            radec_center = np.average(
                catalog[radec_col], weights=1.0 / (catalog[radec_err_col] ** 2)
            )
            dec_center = np.average(
                catalog[dec_col], weights=1.0 / (catalog[dec_err_col] ** 2)
            )

        # Extract coordinates and velocities, applying the plate_scale factor
        x = np.asarray(
            (catalog[radec_col] - radec_center) * plate_scale / scale_divisor
        )
        y = np.asarray(
            (catalog[dec_col] - dec_center) * plate_scale / scale_divisor
        )
        u_pm = np.asarray(catalog[pm_radec_col] * plate_scale)
        v_pm = np.asarray(catalog[pm_dec_col] * plate_scale)
        u_err = np.asarray(catalog[pm_radec_err_col] * plate_scale)
        v_err = np.asarray(catalog[pm_dec_err_col] * plate_scale)

    clust_membership = (
        catalog["Membership"].data * catalog["Membership_color"].data
    )

    # Center bulk motion relative to overall weighted mean
    weights_overall = clust_membership if weighted else np.ones_like(u_pm)
    u_centered = u_pm - np.average(u_pm, weights=weights_overall)
    v_centered = v_pm - np.average(v_pm, weights=weights_overall)

    # Convert Cartesian (x, y) to Polar (r, theta) coordinates
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)  # Angles in radians from -pi to +pi

    # Keep raw centered velocities for optional histograms before tangential projection
    u_raw_centered = u_centered.copy()
    v_raw_centered = v_centered.copy()
    v_tan = -u_centered * np.sin(theta) + v_centered * np.cos(theta)

    # If requested, convert centered velocities to tangential velocity vectors
    if plot_type == "tangential":
        v_tan_err = np.sqrt(
            (np.sin(theta) * u_err) ** 2 + (np.cos(theta) * v_err) ** 2
        )

        u_centered = -v_tan * np.sin(theta)
        v_centered = v_tan * np.cos(theta)
        u_err = np.abs(v_tan_err * np.sin(theta))
        v_err = np.abs(v_tan_err * np.cos(theta))

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
            r,
            theta,
            weights_overall,
            statistic="sum",
            bins=[r_bins, theta_bins],
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            grid_u_mean = u_weighted_sum / weight_sum
            grid_v_mean = v_weighted_sum / weight_sum

            if error_type == "std":
                u_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (u_centered**2) * weights_overall,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_sq_w_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (v_centered**2) * weights_overall,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                u_var = (u_sq_w_sum / weight_sum) - (grid_u_mean**2)
                v_var = (v_sq_w_sum / weight_sum) - (grid_v_mean**2)

                grid_u_sigma_err = np.sqrt(np.maximum(0.0, u_var))
                grid_v_sigma_err = np.sqrt(np.maximum(0.0, v_var))
            else:
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (weights_overall * u_err) ** 2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    (weights_overall * v_err) ** 2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
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
                grid_u_sigma_err, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    u_centered,
                    statistic="std",
                    bins=[r_bins, theta_bins],
                )
                grid_v_sigma_err, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    v_centered,
                    statistic="std",
                    bins=[r_bins, theta_bins],
                )
            else:
                u_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    u_err**2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                v_err_sq_sum, _, _, _ = binned_statistic_2d(
                    r,
                    theta,
                    v_err**2,
                    statistic="sum",
                    bins=[r_bins, theta_bins],
                )
                grid_u_sigma_err = np.sqrt(u_err_sq_sum) / counts
                grid_v_sigma_err = np.sqrt(v_err_sq_sum) / counts

    # Convert velocity components to units of sigma
    with np.errstate(divide="ignore", invalid="ignore"):
        grid_u = grid_u_mean / grid_u_sigma_err
        grid_v = grid_v_mean / grid_v_sigma_err

    # Compute total magnitude in units of sigma for each 2D bin
    grid_sigma_mag = np.hypot(grid_u, grid_v)

    # Compute midpoints for each wedge in (r, theta) to place the vectors
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
    grid_r, grid_theta = np.meshgrid(r_centers, theta_centers, indexing="ij")

    # Convert wedge centerpoints back to Cartesian coordinates for plotting
    grid_x = grid_r * np.cos(grid_theta)
    grid_y = grid_r * np.sin(grid_theta)

    # Set frame-appropriate axis labels
    unit_str = "mas" if scale_divisor == 1.0 else "arcsec"
    if transform_to_sky and to_galactic:
        xlabel, ylabel = r"$\Delta l \cos b$ [arcsec]", r"$\Delta b$ [arcsec]"
        x_sym, y_sym = r"\Delta l \cos b", r"\Delta b"
    elif transform_to_sky and not to_galactic:
        xlabel, ylabel = f"RAcosDEC [{unit_str}]", f"DEC [{unit_str}]"
        x_sym, y_sym = r"\mathrm{RA}\cos\mathrm{DEC}", r"\mathrm{DEC}"
    else:
        xlabel, ylabel = f"X [{unit_str}]", f"Y [{unit_str}]"
        x_sym, y_sym = "X", "Y"

    # Optionally save histograms of star values inside each polar bin
    if save_bin_histograms:
        os.makedirs(hist_dir, exist_ok=True)
        xaxis_inverted = (not transform_to_sky) or to_galactic

        # Order theta_bins clockwise as displayed on the plot, starting from Top (+Y)
        # When X-axis is NOT inverted (East left, -X left): +Y is pi/2, clockwise goes toward +X (0 -> -pi/2 -> +/-pi -> pi/2)
        # When X-axis IS inverted (East left, +X left): +Y is pi/2, clockwise goes toward -X (+/-pi -> -pi/2 -> 0 -> pi/2)
        if xaxis_inverted:
            cw_angles = (theta_centers - np.pi / 2.0) % (2.0 * np.pi)
        else:
            cw_angles = (np.pi / 2.0 - theta_centers) % (2.0 * np.pi)
        cw_order = np.argsort(cw_angles)

        for i in range(nbins_r):
            r_lo, r_hi = r_bins[i], r_bins[i + 1]
            in_r = (r >= r_lo) & (r <= r_hi if i == nbins_r - 1 else r < r_hi)

            for cw_idx, j in enumerate(cw_order):
                t_lo, t_hi = theta_bins[j], theta_bins[j + 1]
                in_t = (theta >= t_lo) & (
                    theta <= t_hi if j == nbins_theta - 1 else theta < t_hi
                )
                in_bin = in_r & in_t
                n_stars = int(np.sum(in_bin))

                if n_stars == 0:
                    continue

                cell_x = grid_x[i, j]
                cell_y = grid_y[i, j]

                fig_h, ax_h = plt.subplots(figsize=(10, 5))

                if plot_type == "tangential":
                    vals = v_tan[in_bin]
                    ax_h.hist(
                        vals,
                        bins="sturges",
                        density=False,
                        histtype="step",
                        linewidth=2,
                        color="black",
                    )
                    ax_h.set_xlabel(r"$v_{\mathrm{tan}}$ [mas/yr]", fontsize=13)
                else:
                    vals_u = u_raw_centered[in_bin]
                    vals_v = v_raw_centered[in_bin]
                    ax_h.hist(
                        vals_u,
                        bins="sturges",
                        density=False,
                        histtype="step",
                        linewidth=2,
                        color="black",
                        label="X PM",
                    )
                    ax_h.hist(
                        vals_v,
                        bins="sturges",
                        density=False,
                        histtype="step",
                        linewidth=2,
                        color="black",
                        linestyle="--",
                        label="Y PM",
                    )
                    ax_h.set_xlabel("Proper Motion [mas/yr]", fontsize=13)
                    ax_h.legend(fontsize=11)

                ax_h.set_ylabel("Number of Stars", fontsize=13)
                ax_h.set_title(
                    f"Cell Center: ({x_sym}={cell_x:.2f}, {y_sym}={cell_y:.2f}) [{unit_str}]; N = {n_stars}",
                    fontsize=13,
                )
                ax_h.grid(True, linestyle=":", alpha=0.5)
                fig_h.tight_layout()

                hist_path = os.path.join(
                    hist_dir, f"r_bin_{i}_theta_bin_{cw_idx}.png"
                )
                fig_h.savefig(hist_path, dpi=200)
                plt.close(fig_h)

    # Mask low-density bins and invalid divisions
    valid_mask = (
        (counts >= min_stars_per_bin)
        & np.isfinite(grid_u)
        & np.isfinite(grid_v)
    )
    plot_x = grid_x[valid_mask]
    plot_y = grid_y[valid_mask]
    plot_u = grid_u[valid_mask]
    plot_v = grid_v[valid_mask]
    plot_counts = counts[valid_mask]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8.5, 8))

    # Highlight sectors where vector significance exceeds highlight_threshold (e.g., > 1 sigma)
    highlight_added = False
    for i in range(nbins_r):
        for j in range(nbins_theta):
            if valid_mask[i, j] and grid_sigma_mag[i, j] > highlight_threshold:
                r1, r2 = r_bins[i], r_bins[i + 1]
                t1 = np.degrees(theta_bins[j])
                t2 = np.degrees(theta_bins[j + 1])

                wedge = mpatches.Wedge(
                    center=(0, 0),
                    r=r2,
                    theta1=t1,
                    theta2=t2,
                    width=(r2 - r1),
                    color=highlight_color,
                    alpha=highlight_alpha,
                    zorder=1,
                )
                ax.add_patch(wedge)
                highlight_added = True

    # Plot azimuthal grid outlines
    for r_edge in r_bins:
        circle = plt.Circle(
            (0, 0),
            r_edge,
            color="gray",
            linestyle="--",
            fill=False,
            alpha=0.5,
            zorder=2,
        )
        ax.add_patch(circle)

    for t_edge in theta_bins:
        ax.plot(
            [0, r_upper * np.cos(t_edge)],
            [0, r_upper * np.sin(t_edge)],
            color="gray",
            linestyle="--",
            alpha=0.5,
            zorder=2,
        )

    # Plot velocity vectors color-coded by number of stars (counts)
    quiver = ax.quiver(
        plot_x,
        plot_y,
        plot_u,
        plot_v,
        plot_counts,
        cmap=cmap,
        angles="xy",
        scale_units="xy",
        scale=scale,
        alpha=0.9,
        width=0.005,
        zorder=3,
    )

    # Add Colorbar for Star Counts
    cbar = fig.colorbar(quiver, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Number of stars per bin ($N$)", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    # Add Quiver Key
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
    ax.plot(0, 0, "r+", markersize=10, label="Center", zorder=4)

    # Add legend entry for highlighted sectors if any are drawn
    if highlight_added:
        highlight_patch = mpatches.Patch(
            color=highlight_color,
            alpha=highlight_alpha,
            label=f"$>{highlight_threshold:.0f}\\sigma$ Motion",
        )
        handles, labels = ax.get_legend_handles_labels()
        handles.append(highlight_patch)
        ax.legend(handles=handles, loc="upper left", fontsize=11)

    # Set boundaries and formatting
    ax.set_aspect("equal")
    ax.set_xlim(-r_upper * 1.05, r_upper * 1.05)
    ax.set_ylim(-r_upper * 1.05, r_upper * 1.05)

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(labelsize=14)

    # Invert X-axis when:
    # 1. transform_to_sky=False (raw detector X increases right)
    # 2. to_galactic=True (l_gal_offset increases positively toward +l, which sits to the left on sky)
    if (not transform_to_sky) or to_galactic:
        ax.invert_xaxis()

    plt.tight_layout()
    fig.savefig(filename_out + ".png", dpi=300)

    return fig, ax