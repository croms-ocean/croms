#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Gravitational OVERFLOW test case.",
    epilog="Example usage:\n"
    + "  python plot_overflow.py --makepdf --file over_his.nc\n"
    + "  python plot_overflow.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="over_his.nc",
    help="Path to the NetCDF file (default: over_his.nc)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate a PDF of the plots"
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate a PNG of the plots"
)
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
tndx = 6  # Time index to process (convert from MATLAB's 1-based index)
h = nc.variables["h"][:]
time = nc.variables["scrum_time"][tndx] / 86400  # Convert time to days
y = nc.variables["y_rho"][:, 1]  # Second column (MATLAB: `(:,2)` -> Python: `[:,1]`)
zeta = nc.variables["zeta"][tndx, :, :]
t0 = nc.variables["temp"][
    0, :, :, 1
]  # First time step, third column (MATLAB: `(:,:,2)` -> Python: `[:,:,1]`)
t = nc.variables["temp"][tndx, :, :, 1]
N, M = t.shape
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
vtrans = np.squeeze(nc.variables.get("Vtransform", None))
nc.close()

# Compute depths using zlevs
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtrans)
zr = zr[:, :, 0]  # First column (MATLAB: `(:,:,1)` -> Python: `[:,:,0]`)
yr = np.tile(y / 1000, (N, 1))  # Convert y to km and replicate for contour plots

# Plot
plt.figure(figsize=(10, 8))

# First subplot: Initial temperature profile
plt.subplot(2, 1, 1)
contour = plt.contourf(yr, zr, t0, levels=np.arange(-0.1, 1.1, 0.1), cmap="viridis")
plt.colorbar(contour)
plt.clim(0, 1)
plt.title("OVERFLOW - Initial Temperature Profile")
plt.xlabel("Y [km]")
plt.ylabel("Z [m]")
plt.text(10, -30, "t=0", fontsize=14)
plt.grid()

# Second subplot: Temperature profile at selected time step
plt.subplot(2, 1, 2)
contour = plt.contourf(yr, zr, t, levels=np.arange(-0.1, 1.1, 0.1), cmap="viridis")
plt.colorbar(contour)
plt.clim(0, 1)
plt.title(f"OVERFLOW - Temperature Profile at t={time:.1f} days")
plt.xlabel("Y [km]")
plt.ylabel("Z [m]")
plt.text(10, -30, f"t={time:.1f} days", fontsize=14)
plt.grid()

# Save to files if requested
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

if args.makepdf:
    pdf_path = os.path.join(output_dir, "overflow_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "overflow_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Display plots unless suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (--no-show used).")
