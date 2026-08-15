#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot and analyze RIVER test case output from a CROCO NetCDF file.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="river_his.nc",
    help="Path to the NetCDF file (default: river_his.nc)",
)
parser.add_argument(
    "--vname", type=str, default="salt", help="Variable name to plot (default: salt)"
)
parser.add_argument(
    "--step",
    type=int,
    default=2,
    help="Step size for vector field quivers (default: 2)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate PDF files for plots"
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate PNG files for plots"
)
parser.add_argument(
    "--makemovie", action="store_true", help="Generate an MP4 movie from the plots"
)
parser.add_argument("--no-show", action="store_true", help="Suppress displaying plots")
parser.add_argument(
    "--output-dir",
    type=str,
    default="./",
    help="Directory to save output files (default: current directory)",
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Parameters
vname = args.vname
step = args.step
makepdf = args.makepdf
makepng = args.makepng
makemovie = args.makemovie
no_show = args.no_show
nc_file = args.file

# Read NetCDF file
try:
    nc = Dataset(nc_file, "r")
except FileNotFoundError:
    print(f"Error: File '{nc_file}' not found.")
    exit(1)

# Extract variables
tis = nc.variables["scrum_time"][:] / (24 * 3600)  # Convert time to days
h = nc.variables["h"][:]
x = nc.variables["x_rho"][:] / 1000  # Convert to km
y = nc.variables["y_rho"][:] / 1000  # Convert to km
N = len(nc.dimensions["s_rho"])
mask = nc.variables["mask_rho"][:]
mask[mask == 0] = np.nan
hmax = np.nanmax(h)

# Prepare for movie creation
if makemovie:
    from matplotlib.animation import FFMpegWriter

    writer = FFMpegWriter(fps=5, metadata={"title": "RIVER Animation"})
    movie_path = os.path.join(args.output_dir, "river_movie.mp4")
    fig, ax = plt.subplots()

# Determine time steps to process
time_steps = range(len(tis)) if makemovie else [len(tis) - 1]

# Loop over time steps
for tndx in time_steps:
    time = tis[tndx]
    s = np.squeeze(nc.variables[vname][tndx, -1, :, :])
    u = np.squeeze(nc.variables["u"][tndx, -1, :, :])
    v = np.squeeze(nc.variables["v"][tndx, -1, :, :])

    # Plot variable with mask
    plt.figure(figsize=(8, 6))
    plt.pcolor(x, y, mask * s, shading="auto", cmap="jet")
    plt.colorbar(label=vname)
    plt.clim(18, 36)  # Color scale for salt
    plt.contour(x, y, h, levels=np.arange(0, 301, 25), colors="k")

    # Quiver plot for velocity vectors
    u_rho = cr.u2rho_2d(u)  # Convert u to rho-points
    v_rho = cr.v2rho_2d(v)  # Convert v to rho-points
    plt.quiver(
        x[::step, ::step],
        y[::step, ::step],
        u_rho[::step, ::step],
        v_rho[::step, ::step],
        scale=10,
        color="k",
    )

    # Format plot
    plt.title(f"RIVER: {vname} - Day {time:.2f}")
    plt.xlabel("X [km]")
    plt.ylabel("Y [km]")
    plt.axis("image")
    plt.axis([0, 40, 0, 80])

    # Save individual frame
    if makepng:
        png_path = os.path.join(args.output_dir, f"river_{vname}_{tndx:03d}.png")
        plt.savefig(png_path, dpi=300, transparent=True)
    if makepdf:
        pdf_path = os.path.join(args.output_dir, f"river_{vname}_{tndx:03d}.pdf")
        plt.savefig(pdf_path, dpi=300, transparent=True)

    # Add frame to movie
    if makemovie:
        plt.savefig("temp_frame.png", dpi=300)  # Save temporary frame
        img = plt.imread("temp_frame.png")
        ax.clear()
        ax.imshow(img)
        writer.grab_frame()

    # Show or close plot
    if no_show:
        plt.close()
    else:
        plt.show()

# Close NetCDF files
nc.close()

# Finalize movie
if makemovie:
    writer.finish()
    print(f"Movie saved to {movie_path}")
