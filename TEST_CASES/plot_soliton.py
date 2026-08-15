#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SOLITON test case.",
    epilog="Example usage:\n"
    + "  python plot_soliton.py --makepdf --file soliton_his.nc\n"
    + "  python plot_soliton.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="soliton_his.nc",
    help="Path to the NetCDF file (default: soliton_his.nc)",
)
parser.add_argument("--makepdf", action="store_true", help="Generate a PDF of the plot")
parser.add_argument("--makepng", action="store_true", help="Generate a PNG of the plot")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Parameters
tndx = 15  # Time index (adjusted for Python's 0-based indexing)
H = 1  # Non-dimensional height (corresponding to 40 cm in physical terms)
L = 1  # Non-dimensional length (corresponding to 295 km in physical terms)
T = 1  # Non-dimensional time (corresponding to 1.71 days in physical terms)

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read variables
time = nc.variables["scrum_time"][tndx] * T
x = nc.variables["x_rho"][:] * L
y = nc.variables["y_rho"][:] * L
z1 = np.squeeze(nc.variables["zeta"][tndx, :, :]) * H
z0 = np.squeeze(nc.variables["zeta"][0, :, :]) * H
nc.close()

# Create plot
plt.figure(figsize=(10, 8))

# Contour plot of zeta at time t
cnt = np.arange(-0.05, 0.21, 0.02) * H
contour = plt.contourf(x, y, z1, levels=cnt, cmap="viridis", linestyles="none")
plt.colorbar(contour, orientation="horizontal", label="Zeta")

# Overlay contours for the initial state
plt.contour(x, y, z0, levels=cnt, colors="k")

# Format plot
plt.axis("image")
plt.title(f"SOLITON - Zeta at t={time:.2f} [non-dimensional]")
plt.xlabel("X [non-dimensional]")
plt.ylabel("Y [non-dimensional]")
plt.grid()

# Save outputs

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

if args.makepdf:
    pdf_path = os.path.join(output_dir, "soliton_plot.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "soliton_plot.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

# Compute theoretical and CROCO soliton properties
maxzi = np.max(z0)
maxzf = np.max(z1)
Xi = np.mean(x[np.where(z0 == maxzi)])
Xf = np.mean(x[np.where(z1 == maxzf)])
X0 = 48 * L
X_croco = X0 - (Xf - Xi)

print(f"Final amplitude (vs. theory): {maxzf:.3f} ({maxzi:.3f})")
print(f"Final position (vs. theory): {X_croco:.3f} ({X0:.3f})")
