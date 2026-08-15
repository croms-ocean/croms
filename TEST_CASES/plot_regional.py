#!/usr/bin/env python3
"""
CROCO Regional Visualization Tool

This script creates visualizations from CROCO model NetCDF outputs including:
- Horizontal maps of various variables (SSH, SST, SSS, velocities, and stresses)
- Vertical sections for temperature, salinity, and velocity components

Usage examples:
  python plot_regional.py --makepdf --file croco_his.nc
  python plot_regional.py --makepng --no-show --output-dir ./plots
  python plot_regional.py --file croco_his.nc --vars temp,salt --sections x
"""

import os
import math
import argparse
import logging
from typing import Dict, List, Tuple, Any
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import croco_utils as cr


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Plot data from a NetCDF file and optionally save as PDF or PNG.",
        epilog="Example usage:\n"
        + "  python plot_regional.py --makepdf --file croco_his.nc\n"
        + "  python plot_regional.py --makepng --no-show",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        type=str,
        default="CROCO_FILES/croco_his.nc",
        help="Path to the NetCDF file (default: CROCO_FILES/croco_his.nc)",
    )
    parser.add_argument(
        "--makepdf", action="store_true", help="Generate a PDF of the plots"
    )
    parser.add_argument(
        "--makepng", action="store_true", help="Generate a PNG of the plots"
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display the plots on the screen (default: display)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Directory to save the output files (default: current directory)",
    )
    parser.add_argument(
        "--ind-sec",
        type=int,
        default=17,
        help="Column/row index for vertical sections (default: 17)",
    )
    parser.add_argument(
        "--vars",
        type=str,
        default="temp,salt,u,v",
        help="Comma-separated list of variables to plot in vertical sections (default: temp,salt,u,v)",
    )
    parser.add_argument(
        "--sections",
        type=str,
        default="x,y",
        help="Comma-separated list of section directions to plot (x for zonal, y for meridional, default: x,y)",
    )
    return parser.parse_args()


def setup_output_directory(output_dir: str) -> str:
    """
    Create the output directory if it doesn't exist.

    Args:
        output_dir (str): Path to the output directory

    Returns:
        str: Path to the created output directory
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Output directory set to: {output_dir}")
        return output_dir
    except Exception as e:
        logger.error(f"Failed to create output directory: {e}")
        raise


def load_netcdf_data(filepath: str) -> Tuple[Dataset, Dict[str, np.ndarray]]:
    """
    Load data from a NetCDF file safely.

    Args:
        filepath (str): Path to the NetCDF file

    Returns:
        Tuple[Dataset, Dict]: NetCDF dataset and dictionary with grid data

    Raises:
        FileNotFoundError: If the file doesn't exist
        IOError: If there's an error reading the file
    """
    if not os.path.exists(filepath):
        logger.error(f"File not found: {filepath}")
        raise FileNotFoundError(f"File '{filepath}' not found.")

    try:
        nc = Dataset(filepath)
        logger.info(f"Successfully loaded NetCDF file: {filepath}")

        # Read lat, lon, mask for different grid points
        grid_data = {}
        for grid_type in ["r", "u", "v"]:
            lat, lon, mask = cr.read_latlonmask(filepath, grid_type)
            grid_data[f"lat_{grid_type}"] = lat
            grid_data[f"lon_{grid_type}"] = lon
            grid_data[f"mask_{grid_type}"] = mask

        # Get dimensions
        try:
            temp = nc.variables["temp"][:]
            T = np.shape(temp)[0]  # Number of time steps
            N = np.shape(temp)[1]  # Number of vertical levels
            grid_data["time_steps"] = T
            grid_data["vertical_levels"] = N
            logger.info(f"Model has {T} time steps and {N} vertical levels")
        except Exception as e:
            logger.warning(f"Could not determine dimensions from temperature: {e}")
            grid_data["time_steps"] = 1
            grid_data["vertical_levels"] = 1

        return nc, grid_data
    except Exception as e:
        logger.error(f"Error loading NetCDF file: {e}")
        raise IOError(f"Error reading file '{filepath}': {e}")


def get_horizontal_data(
    nc_file: str, tndx: int = -1, vlev: int = -1
) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, str, str]]:
    """
    Get horizontal slices of various variables from the NetCDF file.

    Args:
        nc_file (str): Path to the NetCDF file
        tndx (int): Time index (default: -1, last time step)
        vlev (int): Vertical level index (default: -1, surface level)

    Returns:
        Dict: Dictionary of variables with their coordinates and metadata
    """
    variables = {}

    try:
        # Get grid information
        lat_r, lon_r, _ = cr.read_latlonmask(nc_file, "r")
        lat_u, lon_u, _ = cr.read_latlonmask(nc_file, "u")
        lat_v, lon_v, _ = cr.read_latlonmask(nc_file, "v")

        # Define variables to extract with their metadata
        var_specs = [
            ("zeta", 0, "r", "SSH", "Sea Surface Height", "m"),
            ("temp", vlev, "r", "SST", "Sea Surface Temperature", "°C"),
            ("salt", vlev, "r", "SSS", "Sea Surface Salinity", "PSU"),
            ("ubar", 0, "u", "U-Bar", "Barotropic U", "m/s"),
            ("vbar", 0, "v", "V-Bar", "Barotropic V", "m/s"),
            ("u", vlev, "u", "U-Surf", "Surface U", "m/s"),
            ("v", vlev, "v", "V-Surf", "Surface V", "m/s"),
            ("sustr", 0, "u", "U-Wstress", "U Wind Stress", "N/m^{2}"),
            ("svstr", 0, "v", "V-Wstress", "V Wind Stress", "N/m^{2}"),
        ]

        # Get data for each variable
        for var_name, level, grid_type, key, title, unit in var_specs:
            try:
                # Select appropriate coordinates based on grid type
                if grid_type == "r":
                    lon, lat = lon_r, lat_r
                elif grid_type == "u":
                    lon, lat = lon_u, lat_u
                else:  # grid_type == "v"
                    lon, lat = lon_v, lat_v

                # Get the horizontal slice
                data = cr.get_hslice(nc_file, nc_file, var_name, tndx, level, grid_type)
                variables[key] = (lon, lat, data, title, unit)
                logger.debug(f"Successfully extracted {var_name} data")
            except Exception as e:
                logger.warning(f"Failed to extract {var_name}: {e}")

        return variables
    except Exception as e:
        logger.error(f"Error getting horizontal data: {e}")
        raise


def create_horizontal_plots(
    variables: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, str, str]],
    ind_sec: int,
) -> plt.Figure:
    """
    Create enhanced horizontal plots for the given variables.

    Args:
        variables (Dict): Dictionary of variables with their coordinates and metadata
        ind_sec (int): Index for marking section lines

    Returns:
        plt.Figure: Matplotlib figure with the created plots
    """
    # Calculate layout
    n_variables = len(variables)
    ncols = 3  # Number of columns
    nrows = (n_variables + ncols - 1) // ncols  # Calculate required number of rows

    # Set up the figure with better aesthetics
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    # Create figure and axes with better size ratio
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 4.5 * nrows),
        subplot_kw={"projection": ccrs.PlateCarree()},
        constrained_layout=True,  # Better spacing
    )
    axes = axes.flatten()  # Flatten for easier indexing

    # Define better colormaps for different variable types
    cmaps = {
        "SSH": "viridis",
        "SST": "RdYlBu_r",
        "SSS": "cmo.haline",
        "U-Bar": "RdBu_r",
        "V-Bar": "RdBu_r",
        "U-Surf": "RdBu_r",
        "V-Surf": "RdBu_r",
        "U-Wstress": "BrBG",
        "V-Wstress": "BrBG",
    }

    # Try to import cmocean colormaps (better for oceanographic data)
    try:
        import cmocean.cm as cmo

        cmaps.update(
            {
                "SSH": cmo.balance,
                "SST": cmo.thermal,
                "SSS": cmo.haline,
                "U-Bar": cmo.delta,
                "V-Bar": cmo.delta,
                "U-Surf": cmo.delta,
                "V-Surf": cmo.delta,
                "U-Wstress": cmo.curl,
                "V-Wstress": cmo.curl,
            }
        )
    except ImportError:
        logger.info("cmocean not available, using default colormaps")
        # If cmocean fails, fall back to standard colormaps
        cmaps = {
            "SSH": "viridis",
            "SST": "RdYlBu_r",
            "SSS": "YlGnBu",
            "U-Bar": "RdBu_r",
            "V-Bar": "RdBu_r",
            "U-Surf": "RdBu_r",
            "V-Surf": "RdBu_r",
            "U-Wstress": "BrBG",
            "V-Wstress": "BrBG",
        }

    # Plot each variable
    for idx, (key, (lon, lat, data, title, units)) in enumerate(variables.items()):
        ax = axes[idx]

        # Better title with clear spacing and font style
        ax.set_title(rf"{title} ($\mathrm{{{units}}}$)", fontweight="bold", pad=10)

        # Use robust min/max values (ignoring outliers) for better color scaling
        valid_data = data[~np.isnan(data)]
        if len(valid_data) > 0:
            vmin = np.percentile(valid_data, 2)
            vmax = np.percentile(valid_data, 98)

            # Center diverging colormaps at zero for velocity and stress fields
            if key in ["U-Bar", "V-Bar", "U-Surf", "V-Surf", "U-Wstress", "V-Wstress"]:
                abs_max = max(abs(vmin), abs(vmax))
                vmin, vmax = -abs_max, abs_max
        else:
            vmin, vmax = np.nanmin(data), np.nanmax(data)

        # Get appropriate colormap
        cmap = cmaps.get(key, "viridis")

        # Create the color mesh plot with improved parameters
        mesh = ax.contourf(
            lon,
            lat,
            data,
            levels=20,
            cmap=cmap,
            norm=Normalize(vmin=vmin, vmax=vmax),
            extend="both",  # To extend color for off-range values
        )

        # Add coastlines and land mask if available
        try:
            ax.coastlines(resolution="50m", color="gray", linewidth=0.75)
            ax.add_feature(cartopy.feature.LAND, facecolor="lightgray", alpha=0.3)
        except Exception:
            logger.debug("Could not add coastlines/land features")

        # Add section lines to the maps with improved visibility
        try:
            # Zonal section line (along x)
            ax.plot(
                lon[ind_sec, :],
                lat[ind_sec, :],
                color="black",
                linestyle="--",
                linewidth=1.5,
                alpha=0.8,
                zorder=10,  # Place above other elements
            )

            # Meridional section line (along y)
            ax.plot(
                lon[:, ind_sec],
                lat[:, ind_sec],
                color="black",
                linestyle="--",
                linewidth=1.5,
                alpha=0.8,
                zorder=10,
            )

            # Add annotations for the section lines with improved styling
            # Get midpoints for placing text
            x_mid = len(lon[ind_sec, :]) // 2
            y_mid = len(lon[:, ind_sec]) // 2

            # Add text annotations with section indices
            ax.text(
                lon[ind_sec, x_mid],
                lat[ind_sec, x_mid],
                f"i={ind_sec}",
                fontsize=9,
                fontweight="bold",
                bbox=dict(
                    facecolor="white",
                    alpha=0.8,
                    edgecolor="black",
                    linewidth=0.5,
                    pad=2,
                ),
                ha="center",
                va="bottom",
                zorder=11,
            )

            ax.text(
                lon[y_mid, ind_sec],
                lat[y_mid, ind_sec],
                f"j={ind_sec}",
                fontsize=9,
                fontweight="bold",
                bbox=dict(
                    facecolor="white",
                    alpha=0.8,
                    edgecolor="black",
                    linewidth=0.5,
                    pad=2,
                ),
                ha="right",
                va="center",
                zorder=11,
            )

        except IndexError:
            logger.warning(f"Section index {ind_sec} out of bounds for variable {key}")

        # Add better grid lines
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.5,
            color="gray",
            alpha=0.5,
            linestyle=":",
        )
        gl.top_labels = False
        gl.right_labels = False

        # Format labels to gain space
        gl.xformatter = cticker.LongitudeFormatter(degree_symbol="°")
        gl.yformatter = cticker.LatitudeFormatter(degree_symbol="°")

        # Adjust label size
        gl.xlabel_style = {"size": 8, "color": "black"}
        gl.ylabel_style = {"size": 8, "color": "black"}

        # Add a nicer color bar
        cbar = plt.colorbar(
            mesh,
            ax=ax,
            orientation="vertical",
            fraction=0.046,
            pad=0.04,
            aspect=30,
            shrink=0.7,
        )
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(units, fontsize=9)

        # Add a subtle border around plots
        ax.spines["geo"].set_edgecolor("gray")
        ax.spines["geo"].set_linewidth(0.5)

    # Remove unused axes
    for idx in range(n_variables, len(axes)):
        fig.delaxes(axes[idx])

    # Add a legend to the first subplot
    handles = [
        plt.Line2D(
            [0],
            [0],
            color="black",
            linestyle="--",
            linewidth=1.5,
            label="Section lines (i,j)",
        ),
    ]
    axes[0].legend(handles=handles, loc="lower left", fontsize=9, framealpha=0.7)

    # Add overall title showing section indices
    plt.suptitle(
        f"Horizontal Maps with Section Lines at Index {ind_sec}",
        fontsize=14,
        y=1.02,
        fontweight="bold",
    )

    return fig


def get_vertical_section_data(
    nc_file: str,
    var_list: List[str],
    direction_list: List[str],
    tindex: int = -1,
    ind_sec: int = 17,
) -> Dict[str, Dict[str, Any]]:
    """
    Get vertical section data for multiple variables and directions.

    Args:
        nc_file (str): Path to the NetCDF file
        var_list (List[str]): List of variables to extract
        direction_list (List[str]): List of directions ('x' or 'y')
        tindex (int): Time index
        ind_sec (int): Section index

    Returns:
        Dict: Dictionary with section data for each variable and direction
    """
    section_data = {}

    for var in var_list:
        section_data[var] = {}
        for direction in direction_list:
            try:
                section = cr.get_vertical_section(
                    nc_file, var, tindex, direction, ind_sec, nc_file
                )
                section_data[var][direction] = section
                logger.debug(f"Extracted {var} section in {direction} direction")
            except Exception as e:
                logger.warning(
                    f"Failed to extract {var} section in {direction} direction: {e}"
                )
                section_data[var][direction] = None

    return section_data


def create_vertical_section_plots(
    section_data: Dict[str, Dict[str, Any]], ind_sec: int
) -> plt.Figure:
    """
    Create vertical section plots for the given data.

    Args:
        section_data (Dict): Dictionary with section data
        ind_sec (int): Section index for the title

    Returns:
        plt.Figure: Matplotlib figure with the created plots
    """
    # Define titles and units for variables
    titles = {
        "temp": "Temperature",
        "salt": "Salinity",
        "u": "Zonal Velocity",
        "v": "Meridional Velocity",
    }

    units = {
        "temp": "°C",
        "salt": "PSU",
        "u": "m/s",
        "v": "m/s",
    }

    # Count total number of sections
    n_sections = sum(
        sum(1 for direction in directions.keys() if directions[direction] is not None)
        for var, directions in section_data.items()
    )

    if n_sections == 0:
        logger.warning("No valid section data to plot")
        return None

    # Calculate layout
    n_cols = min(2, n_sections)  # Maximum 2 columns
    n_rows = math.ceil(n_sections / n_cols)

    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 3))

    # Convert to array if only one subplot
    if n_sections == 1:
        axes = np.array([axes])

    # Flatten the axes array for easier indexing
    axes = axes.flatten()

    # Initialize axis index
    axis_index = 0

    # Define colormaps for different variables
    colormaps = {
        "temp": "RdYlBu_r",
        "salt": "YlGnBu",
        "u": "RdBu_r",
        "v": "RdBu_r",
    }

    # Try to import cmocean colormaps (better for oceanographic data)
    try:
        import cmocean.cm as cmo

        colormaps.update(
            {
                "temp": cmo.thermal,
                "salt": cmo.haline,
                "u": cmo.delta,
                "v": cmo.delta,
            }
        )
    except ImportError:
        logger.info("cmocean not available, using default colormaps")
        # If cmocean fails, fall back to standard colormaps
        colormaps = {
            "temp": "RdYlBu_r",
            "salt": "YlGnBu",
            "u": "RdBu_r",
            "v": "RdBu_r",
        }

    # Plot each variable and direction
    for var, directions in section_data.items():
        for direction, section in directions.items():
            if section is None:
                continue

            # Get data and mask zeros
            data = section["variable"]
            m_data = np.ma.masked_equal(data, 0)

            # Configure the axes
            ax = axes[axis_index]

            # Get colormap
            cmap = colormaps.get(var, "viridis")

            # Set color limits based on data type
            if var in ["u", "v"]:
                # Symmetric colormap for velocity
                abs_max = np.max(np.abs(m_data))
                vmin, vmax = -abs_max, abs_max
            else:
                # Use robust statistics for other variables
                valid_data = m_data.compressed()
                if len(valid_data) > 0:
                    vmin = np.percentile(valid_data, 2)
                    vmax = np.percentile(valid_data, 98)
                else:
                    vmin, vmax = np.nanmin(m_data), np.nanmax(m_data)

            # Plot the data
            im = ax.pcolormesh(
                section["distance"],
                section["depth"],
                m_data,
                cmap=cmap,
                shading="auto",
                vmin=vmin,
                vmax=vmax,
            )

            # Add grid
            ax.grid(True, linestyle="--", alpha=0.7)

            # Add bathymetry/topography line
            topo_profile = section["topo"]
            distance_profile = section["distance"]
            ax.plot(
                distance_profile,
                topo_profile,
                color="black",
                linestyle="-",
                linewidth=1.5,
            )

            # Set title and labels
            direction_name = "Zonal" if direction == "x" else "Meridional"
            var_title = titles.get(var, var.capitalize())
            var_unit = units.get(var, "")

            # Add the index to the title
            index_label = f"i={ind_sec}" if direction == "y" else f"j={ind_sec}"
            ax.set_title(
                f"{var_title} ({var_unit}) - {direction_name} Section at {index_label}"
            )

            if direction == "x":
                ax.set_xlabel("Distance along longitude (km)")
            elif direction == "y":
                ax.set_xlabel("Distance along latitude (km)")

            ax.set_ylabel("Depth (m)")

            # Add section index as text annotation on the plot
            # Position at top-right corner
            ax.text(
                0.95,
                0.95,
                index_label,
                transform=ax.transAxes,
                horizontalalignment="right",
                verticalalignment="top",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=3),
                fontsize=10,
                weight="bold",
            )

            # Add a color bar
            cbar = fig.colorbar(im, ax=ax, orientation="vertical")
            cbar.set_label(var_unit)

            # Update axis index
            axis_index += 1

    # Remove unused axes
    for j in range(n_sections, len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout
    plt.tight_layout()

    # Add overall title
    plt.suptitle(f"Vertical Sections at Index {ind_sec}", fontsize=16, y=1.02)

    return fig


def save_figure(
    fig: plt.Figure,
    output_dir: str,
    basename: str,
    save_pdf: bool = False,
    save_png: bool = False,
) -> None:
    """
    Save a figure to disk in PDF and/or PNG format.

    Args:
        fig (plt.Figure): Figure to save
        output_dir (str): Directory to save to
        basename (str): Base filename
        save_pdf (bool): Whether to save as PDF
        save_png (bool): Whether to save as PNG
    """
    if save_pdf:
        pdf_path = os.path.join(output_dir, f"{basename}.pdf")
        try:
            fig.savefig(pdf_path, transparent=True, bbox_inches="tight")
            logger.info(f"PDF file '{pdf_path}' has been created.")
        except Exception as e:
            logger.error(f"Failed to save PDF: {e}")

    if save_png:
        png_path = os.path.join(output_dir, f"{basename}.png")
        try:
            fig.savefig(png_path, dpi=600, bbox_inches="tight")
            logger.info(f"PNG file '{png_path}' has been created.")
        except Exception as e:
            logger.error(f"Failed to save PNG: {e}")


def main():
    """Main function to run the CROCO visualization tool."""
    # Parse arguments
    args = parse_arguments()

    # Set up output directory
    output_dir = setup_output_directory(args.output_dir)

    try:
        # Load NetCDF file
        logger.info(f"Loading NetCDF file: {args.file}")
        nc, grid_data = load_netcdf_data(args.file)

        # Get vertical levels
        N = grid_data["vertical_levels"]

        # Process horizontal plots
        if True:  # Always process horizontal plots (can be changed to a parameter)
            logger.info("Creating horizontal plots")
            variables = get_horizontal_data(
                args.file,
                tndx=-1,  # Last time step
                vlev=N - 1,  # Surface level (in Python indexing)
            )

            horiz_fig = create_horizontal_plots(variables, args.ind_sec)

            # Save horizontal plots
            save_figure(
                horiz_fig,
                output_dir,
                "REGIONAL_maps",
                save_pdf=args.makepdf,
                save_png=args.makepng,
            )

        # Process vertical section plots
        if True:  # Always process vertical plots (can be changed to a parameter)
            logger.info("Creating vertical section plots")

            # Parse variable and direction lists
            var_list = [v.strip() for v in args.vars.split(",")]
            direction_list = [d.strip() for d in args.sections.split(",")]

            # Get section data
            section_data = get_vertical_section_data(
                args.file,
                var_list,
                direction_list,
                tindex=-1,  # Last time step
                ind_sec=args.ind_sec,
            )

            # Create section plots
            section_fig = create_vertical_section_plots(section_data, args.ind_sec)

            if section_fig:
                # Save section plots
                save_figure(
                    section_fig,
                    output_dir,
                    "REGIONAL_sections",
                    save_pdf=args.makepdf,
                    save_png=args.makepng,
                )

        # Show plots if requested
        if not args.no_show:
            logger.info("Showing plots")
            plt.show()
        else:
            logger.info("Plot display suppressed")

    except Exception as e:
        logger.error(f"Error processing data: {e}", exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
