#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr


# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot data from a NetCDF file and optionally save as PDF or PNG.",
    epilog="Example usage:\n"
    + "  python script.py --makepdf --file basin_his.nc\n"
    + "  python script.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="basin_his.nc",
    help="Path to the NetCDF file (default: basin.nc)",
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

args = parser.parse_args()

# Initialization
tndx = -1
j = 25

# Read data from NetCDF file
try:
    nc = Dataset(args.file)
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

time = nc.variables["scrum_time"][tndx] / 86400
h = nc.variables["h"][:]
x1 = nc.variables["x_rho"][:]
y1 = nc.variables["y_rho"][:]
x = x1[j, :]
zeta = nc.variables["zeta"][tndx, :, :]
t = nc.variables["temp"][tndx, :, j, :]
R0, TCOEF = 30, 0.28
rho = R0 - TCOEF * t
N, M = t.shape
theta_s = nc.theta_s
theta_b = nc.theta_b
vtrans = np.squeeze(nc.variables.get("Vtransform", None))
hc = nc.hc
nc.close()

zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtrans)
zr = zr[:, j, :]
xr = np.tile(x / 1000, (N, 1))  # Convert to km

# First plot
plt.figure(figsize=(6, 9))
plt.subplot(2, 1, 1)
contour = plt.contourf(xr, zr, rho, np.arange(29, 30.1, 0.1), cmap="viridis")
plt.colorbar(contour)
plt.clim(29, 30)
plt.xlabel("X [km]")
plt.ylabel("Z [m]")
plt.title(f"BASIN - $\\sigma_t$ [kg/m^3] vertical section at {time:.2f} days")

# Second plot
plt.subplot(2, 1, 2)
contour = plt.contourf(
    x1[1:-1, 1:-1] / 1000,
    y1[1:-1, 1:-1] / 1000,
    100 * zeta[1:-1, 1:-1],
    np.arange(-20, 22, 2),
    cmap="RdYlBu",
)
plt.clim(-8, 8)
plt.colorbar(contour)
plt.xlabel("X [km]")
plt.ylabel("Y [km]")
plt.title(f"BASIN - sea surface elevation [cm] at {time:.2f} days")

# Determine output file paths
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(output_dir, "basin_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(output_dir, "basin_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Show plots if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --show to enable).")
