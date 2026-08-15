#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import cm, colors
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the INTERNAL test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="internal_his.nc", help="Path to the NetCDF file"
)
parser.add_argument("--makepdf", action="store_true", help="Generate a PDF of the plot")
parser.add_argument("--makepng", action="store_true", help="Generate a PNG of the plot")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save the output file"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Parameters
j = 2  # Adjust for Python's 0-based indexing (MATLAB j=3)

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read data
time = nc.variables["scrum_time"][:] / 86400  # Convert time to days
tndx = len(time)  # Number of time steps
print(f"tndx = {tndx} - Time = {time[-1] * 24 / 12.4:.2f} M2 cycles")
h = np.squeeze(nc.variables["h"][j, :])
x = np.squeeze(nc.variables["x_rho"][j, :])
zeta = np.squeeze(nc.variables["zeta"][tndx - 1, j, :])
t = np.squeeze(nc.variables["rho"][tndx - 1, :, j, :])
t0 = np.squeeze(nc.variables["rho"][0, :, j, :])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
vtransform = np.squeeze(nc.variables.get("Vtransform", 1))  # Default to 1 if not found
nc.close()

# Dimensions
N, L = t.shape

# Compute vertical levels
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtransform)
xr = np.tile(x / 1000, (N, 1))  # Convert x to km

# Create a custom colormap
colormap = cm.jet(np.linspace(0, 1, 40))  # Generate colormap array
colormap[19:21, :] = [1, 1, 1, 1]  # Add white band around zero
custom_cmap = colors.ListedColormap(colormap)

# Plot rho anomaly
plt.figure(figsize=(10, 5))
contour = plt.contourf(xr, zr, t - t0, levels=40, cmap=custom_cmap)
plt.colorbar(contour, label="Rho Anomaly")
plt.plot(x / 1000, -h, color="k", linewidth=2, label="Bathymetry")
plt.clim(-0.01, 0.01)

# Add plot details
plt.title(f"Internal Case - Rho Anomaly at {time[-1] * 24 / 12.4:.2f} M2 cycles")
plt.xlabel("X [km]")
plt.ylabel("Depth [m]")
plt.legend()
plt.grid()
plt.tight_layout()

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "internal_rho.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved to {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "internal_rho.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG saved to {png_path}")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()
