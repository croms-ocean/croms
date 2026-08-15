#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr


# Command-line options
parser = argparse.ArgumentParser(description="Generate plots for the CANYON test case.")
parser.add_argument(
    "--file", type=str, default="canyon_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate a PDF of the plots"
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate a PNG of the plots"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir",
    type=str,
    default=".",
    help="Directory to save the output files (default: current directory)",
)
args = parser.parse_args()

# Parameters
tndx = 2  # Equivalent to MATLAB's tndx=3
i = 31  # Equivalent to MATLAB's i=32 (adjusted for Python indexing)

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
y = y1[:, i]
zeta = nc.variables["zeta"][tndx, :, :]
t = nc.variables["rho"][tndx, :, :, i]
N, M = t.shape
sst = nc.variables["temp"][tndx, -1, :, :]
u = nc.variables["u"][tndx, -1, :, :]
v = nc.variables["v"][tndx, -1, :, :]
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
vtrans = np.squeeze(nc.variables.get("Vtransform", None))
nc.close()

# Adjust u, v to rho-points
ur = 0.5 * (u[:, :-1] + u[:, 1:])  # Example replacement for MATLAB's u2rho_2d
vr = 0.5 * (v[:-1, :] + v[1:, :])  # Example replacement for MATLAB's v2rho_2d

zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtrans)
zr = zr[:, :, i]
yr = np.tile(y / 1000, (N, 1))

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# First plot
plt.figure(figsize=(6, 9))
plt.subplot(2, 1, 1)
contour = plt.contourf(yr.T, zr.T, t.T, np.arange(28, 32.1, 0.1), cmap="viridis")
plt.colorbar(contour)
plt.title(f"CANYON - $\\sigma_t$ [kg/m^3] vertical section at {time:.2f} days")
plt.clim(28, 31)

# Second plot
plt.subplot(2, 1, 2)
contour = plt.contourf(
    x1[1:-1, 1:-1] / 1000,
    y1[1:-1, 1:-1] / 1000,
    100 * zeta[1:-1, 1:-1],
    np.arange(-0.5, 2.1, 0.1),
    cmap="RdYlBu",
)
plt.colorbar(contour)
plt.clim(-0.1, 1)
plt.contour(
    x1[1:-1, 1:-1] / 1000,
    y1[1:-1, 1:-1] / 1000,
    h[1:-1, 1:-1],
    colors="k",
)
plt.title(f"CANYON - sea surface elevation [cm] at {time:.2f} days")

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "canyon_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(args.output_dir, "canyon_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()
