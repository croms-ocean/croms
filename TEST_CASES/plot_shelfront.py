#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SHELFRONT test case.",
    epilog="Example usage:\n"
    + "  python plot_shelfront.py --makepdf --file shelfront_his.nc\n"
    + "  python plot_shelfront.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="shelfront_his.nc",
    help="Path to the NetCDF file (default: shelfront_his.nc)",
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
tndx = 10  # Adjusted for Python's 0-based indexing (MATLAB's tndx = 11)
time = nc.variables["scrum_time"][tndx] / (24 * 3600)  # Convert time to days
h = nc.variables["h"][:]
y = np.squeeze(
    nc.variables["y_rho"][:, 1]
)  # Second column (MATLAB: `(:,2)` -> Python: `[:,1]`)
zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
t = np.squeeze(
    nc.variables["temp"][tndx, :, :, 1]
)  # Second column (MATLAB: `(:,:,2)` -> Python: `[:,:,1]`)
u = np.squeeze(nc.variables["u"][tndx, :, :, 1])
N, M = t.shape
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
vtransform = np.squeeze(nc.variables.get("Vtransform", 1))  # Default to 1 if not found
nc.close()

# Compute depths using zlevs
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtransform)
zr = zr[:, :, 0]  # First column (MATLAB: `(:,:,1)` -> Python: `[:,:,0]`)
yr = np.reshape(y, (1, M))
yr = np.tile(yr, (N, 1)) / 1000  # Convert to kilometers

# Plot
plt.figure(figsize=(10, 8))

# Contour plot of temperature
contour = plt.contourf(
    yr, zr, t, levels=np.arange(12, 18.5, 0.5), cmap="viridis", linestyles="none"
)
plt.colorbar(contour, label="Temperature [°C]")
plt.clim(12, 18)

# Overlay velocity contours
C1 = plt.contour(yr, zr, 100 * u, levels=np.arange(-5, 10, 2), colors="k")
plt.clabel(C1)

# Plot formatting
plt.title(f"SHELFRONT - Temp [°C] - Day = {time:.1f}")
plt.xlabel("Y [km]")
plt.ylabel("Z [m]")
plt.grid()

# Save to files if requested
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

if args.makepdf:
    pdf_path = os.path.join(output_dir, "shelfront_plot.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "shelfront_plot.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Display plot unless suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (--no-show used).")
