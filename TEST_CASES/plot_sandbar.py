#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SANDBAR test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="sandbar_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--mycase",
    type=str,
    choices=["1B", "1C"],
    default="1B",
    help="LIP experiment case: '1B' (erosion) or '1C' (accretion)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Save the plots as a PDF file"
)
parser.add_argument(
    "--makepng", action="store_true", help="Save the plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Read input parameters
mycase = args.mycase
morph_fac = 18 if mycase == "1B" else 13
morph_cpl = True

# Read grid and parameters
yindex = 1  # Python indexing (MATLAB yindex = 2)
tindex = len(nc.variables["scrum_time"][:]) - 1  # Last time record

time = morph_fac * nc.variables["scrum_time"][:] / 3600  # Time in hours
tindex0 = np.argmin(np.abs(time - (4 if mycase == "1B" else 3)))

# Horizontal grid
hr = np.squeeze(nc.variables["h"][yindex, :])
xr = np.squeeze(nc.variables["x_rho"][yindex, :])
hu = 0.5 * (hr[:-1] + hr[1:])
xu = 0.5 * (xr[:-1] + xr[1:])

# New bathymetry
if morph_cpl:
    hnew = np.squeeze(nc.variables["hmorph"][tindex, yindex, :])
    h0 = np.squeeze(nc.variables["hmorph"][tindex0, yindex, :])
    h = hnew if hnew is not None else hr
else:
    h, hnew, h0 = hr, hr, hr

# Vertical grid
N = len(nc.dimensions["s_rho"])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
zeta = np.squeeze(nc.variables["zeta"][tindex, yindex, :])
Dcrit = 1.1 * nc.variables["Dcrit"][:]

zeta[h < Dcrit] -= h[h < Dcrit]
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
zw = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "w", 2)

zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
zwu = 0.5 * (zw[:, :-1] + zw[:, 1:])

# Model fields
u = np.squeeze(nc.variables["u"][tindex, :, yindex, :])
C = 2 * np.squeeze(nc.variables["sand_01"][tindex, :, yindex, :])
sup = np.squeeze(nc.variables["zeta"][tindex0, yindex, :])
sup[hr < 0] += hr[hr < 0] - Dcrit

# Align dimensions
xu_2d = np.tile(xu, (zru.shape[0], 1))  # Expand xu to match zru

# Plotting
fig, axs = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)

# Plot U velocity
cmin, cmax, nbcol = -0.5, 0.5, 20
cont = axs[0].contourf(
    xu_2d, zru, u, levels=np.linspace(cmin, cmax, nbcol), cmap="jet", extend="both"
)
axs[0].plot(xr, -hr, "k:", linewidth=3, label="Initial")
axs[0].plot(xr, -hnew, "k", linewidth=3, label="Final Model")
axs[0].set_xlim(60, 190)
axs[0].set_ylim(-2.5, 0.5)
axs[0].legend()
axs[0].set_title(f"Sandbar {mycase} - U Velocity at Time {time[tindex]:.2f} hours")

# Plot Hrms
axs[1].plot(xr, sup, "b-", label="Model")
axs[1].legend()
axs[1].set_xlim(60, 190)
axs[1].set_title("Hrms at the Sandbar")
axs[1].grid()

# Save outputs
if args.makepng:
    png_paths = [
        os.path.join(args.output_dir, f"sandbar_{mycase}_u_velocity.png"),
        os.path.join(args.output_dir, f"sandbar_{mycase}_hrms.png"),
    ]
    for ax, path in zip(axs, png_paths):
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path, bbox_inches=extent)
        print(f"PNG saved: {path}")

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, f"sandbar_{mycase}.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved: {pdf_path}")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()

# Close NetCDF file
nc.close()
