#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SEAMOUNT test case.",
    epilog="Example usage:\n"
    + "  python plot_seamount.py --makepdf --file seamount_his.nc\n"
    + "  python plot_seamount.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="seamount_his.nc",
    help="Path to the NetCDF file (default: seamount_his.nc)",
)
parser.add_argument("--makepdf", action="store_true", help="Generate a PDF of the plot")
parser.add_argument("--makepng", action="store_true", help="Generate a PNG of the plot")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Read data from NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Parameters
tndx = 5  # Time index (adjusted for Python's 0-based indexing)
N = 0  # First vertical level (MATLAB 1 -> Python 0)
time = nc.variables["scrum_time"][tndx] / (24 * 3600)  # Convert time to days
h = nc.variables["h"][:]
x = nc.variables["x_rho"][:] / 1000  # Convert to km
y = nc.variables["y_rho"][:] / 1000  # Convert to km
u = np.squeeze(nc.variables["u"][tndx, N, :, :])
v = np.squeeze(nc.variables["v"][tndx, N, :, :])
nc.close()

# Compute bottom speed (magnitude of velocity at the bottom level)
spd = 1000 * np.sqrt(cr.u2rho_2d(u) ** 2 + cr.v2rho_2d(v) ** 2)  # Convert to mm/s

# Plot
plt.figure(figsize=(10, 8))

# Contour plot of speed
contour = plt.contourf(x, y, spd, levels=np.arange(0, 5.5, 0.5), cmap="viridis")
plt.colorbar(contour, label="Speed [mm/s]")
plt.clim(0, 3)

# Overlay bathymetry contours
plt.contour(x, y, h, levels=np.arange(1000, 4500, 500), colors="k")

# Plot formatting
plt.axis("image")
plt.axis([0, 500, 0, 500])
plt.xlabel("X [km]")
plt.ylabel("Y [km]")
plt.title(f"SEAMOUNT - Bottom Speed [mm/s] - Day = {time:.1f}", fontsize=14)
plt.grid()

# Save to files if requested
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

if args.makepdf:
    pdf_path = os.path.join(output_dir, "seamount_plot.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "seamount_plot.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Display plots unless suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (--no-show used).")
