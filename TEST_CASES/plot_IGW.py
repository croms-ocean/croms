#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the IGW test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="igw_his.nc", help="Path to the NetCDF file"
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

# User-defined parameters
jj = 599  # Location of validation (adjusted for Python indexing)
valid = False  # Set to True for validation

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read data
tndx = len(nc.variables["scrum_time"][:]) - 1
h = nc.variables["h"][:]
hsec = np.squeeze(nc.variables["h"][1, :])
lonu = np.squeeze(nc.variables["lon_u"][1, :])
lonr = np.squeeze(nc.variables["lon_rho"][1, :])
N = len(nc.dimensions["s_rho"])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
vtransform = np.squeeze(nc.variables.get("Vtransform", 1))
ssh = np.squeeze(nc.variables["zeta"][:, 1, :])
zeta = np.squeeze(nc.variables["zeta"][tndx, :, :])
u = np.squeeze(nc.variables["ubar"][:, 1, :])
v = np.squeeze(nc.variables["vbar"][:, 1, :])
usec = np.squeeze(nc.variables["u"][tndx, :, 1, :])
vsec = np.squeeze(nc.variables["v"][tndx, :, 1, :])
wsec = np.squeeze(nc.variables["w"][tndx, :, 1, :])
tsec = np.squeeze(nc.variables["temp"][tndx, :, 1, :])
rsec = np.squeeze(nc.variables["rho"][tndx, :, 1, :])
drsec = rsec - np.squeeze(nc.variables["rho"][0, :, 1, :])
time = nc.variables["scrum_time"][tndx] / 86400
nc.close()

# Compute depths
zeta_u = cr.rho2u_2d(zeta)
h_u = cr.rho2u_2d(h)
z = cr.zlevs(h_u, zeta_u, theta_s, theta_b, hc, N, "r", vtransform)
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", vtransform)
zsec = np.squeeze(z[:, 1, :])
xsec = np.tile(lonu, (N, 1))
zrsec = np.squeeze(zr[:, 1, :])
xrsec = np.tile(lonr, (N, 1))

# Plot internal tides section for u, w, and rho anomaly
plt.figure(figsize=(8, 10))
map_jet = plt.cm.jet(np.linspace(0, 1, 20))
map_jet[9:11, :] = [1, 1, 1, 1]  # Add white band
custom_cmap = plt.cm.colors.ListedColormap(map_jet)

# Plot U
plt.subplot(3, 1, 1)
plt.contourf(xsec, zsec, usec, levels=20, cmap=custom_cmap)
plt.colorbar(label="U [m/s]")
plt.plot(lonr, -hsec, color="k", linewidth=3)
plt.clim(-0.5, 0.5)
plt.ylabel("Z [m]")
plt.title(f"IGW - U [m/s] at {time:.1f} days")
plt.grid()

# Plot W
plt.subplot(3, 1, 2)
plt.contourf(xrsec, zrsec, wsec, levels=20, cmap=custom_cmap)
plt.colorbar(label="W [m/s]")
plt.plot(lonr, -hsec, color="k", linewidth=3)
plt.clim(-0.03, 0.03)
plt.ylabel("Z [m]")
plt.title(f"IGW - W [m/s] at {time:.1f} days")
plt.grid()

# Plot Rho Anomaly
plt.subplot(3, 1, 3)
plt.contourf(xrsec, zrsec, drsec, levels=20, cmap=custom_cmap)
plt.colorbar(label="Rho Anomaly [kg/m³]")
plt.plot(lonr, -hsec, color="k", linewidth=3)
plt.clim(-0.05, 0.05)
plt.xlabel("Longitude")
plt.ylabel("Z [m]")
plt.title(f"IGW - ρₐ [kg/m³] at {time:.1f} days")
plt.grid()

plt.tight_layout()

# Save or display plots
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "IGW.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved to {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "IGW.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG saved to {png_path}")

if not args.no_show:
    plt.show()
else:
    plt.close()
