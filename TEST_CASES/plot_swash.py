#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command args
parser = argparse.ArgumentParser(
    description="Plot results from the SWASH test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="swash_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--varname",
    type=str,
    choices=["u", "w"],
    default="u",
    help="Variable to plot: 'u' (zonal velocity) or 'w' (vertical velocity)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Save the plots as a PDF file"
)
parser.add_argument("--makemovie", action="store_true", help="Generate an MP4 movie")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Be sure output directory exist
os.makedirs(args.output_dir, exist_ok=True)

# Open file NetCDF
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Grid parameters
yindex = 1  # Index Python (MATLAB yindex = 2)
tindex_last = len(nc.variables["scrum_time"][:]) - 1  # Last record
tstart = 0 if args.makemovie else tindex_last
tend = tindex_last
g = 9.81

hr = np.squeeze(nc.variables["h"][yindex, :])
xr = np.squeeze(nc.variables["x_rho"][yindex, :])
dx = xr[1] - xr[0]
xmin, xmax = np.min(xr), 90
zmin, zmax = -np.max(hr), 0.4 * np.max(hr)

N = len(nc.dimensions["s_rho"])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
Dcrit = np.squeeze(nc.variables["Dcrit"][:])

# Prepare figure
fig, ax = plt.subplots(figsize=(10, 6))

# Time loop
for tindex in range(tstart, tend + 1):
    print(f"Processing time index: {tindex}")

    time = nc.variables["scrum_time"][tindex]
    zeta = np.squeeze(nc.variables["zeta"][tindex, yindex, :])
    zr = cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "r", 2)
    zw = cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "w", 2)

    zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
    xr_u = 0.5 * (xr[:-1] + xr[1:])  # Grid shift for u
    xr_u_2d = np.tile(xr_u, (zru.shape[0], 1))
    D = hr + zeta
    D_u = 0.5 * (D[:-1] + D[1:])  # Grid shift for D
    D2d_u = np.tile(D_u, (zru.shape[0], 1))

    # Read vars
    if args.varname == "u":
        u = np.squeeze(nc.variables["u"][tindex, :, yindex, :])
        var = u
        cmin, cmax = -1.5, 1.5
    elif args.varname == "w":
        w = np.squeeze(nc.variables["w"][tindex, :, yindex, :])
        var = w
        cmin, cmax = -0.7, 0.7

    var[D2d_u < Dcrit] = np.nan  # Mask

    # Plots
    ax.clear()
    nbcol = 20
    cmap = plt.cm.jet
    cmap.set_under("white")

    cont = ax.contourf(
        xr_u_2d,
        zru,
        var,
        levels=np.linspace(cmin, cmax, nbcol),
        cmap=cmap,
        extend="both",
    )
    fig.colorbar(cont, ax=ax)
    ax.plot(xr, -hr, "k", linewidth=3, label="Bathymetry")
    ztop = zw[-1, :]
    ztop[D < Dcrit + 0.01] = np.nan
    ax.plot(xr, ztop, "r", linewidth=2, label="Surface")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(zmin, zmax)
    ax.set_title(f"SWASH: {args.varname} at {time:.2f} sec", fontsize=14)
    ax.set_xlabel("Distance (m)", fontsize=12)
    ax.set_ylabel("Depth (m)", fontsize=12)
    ax.legend()
    ax.grid()

    # Save outputs
    if args.makemovie:
        plt.savefig(os.path.join(args.output_dir, f"frame_{tindex:04d}.png"))

    if args.makepng:
        png_path = os.path.join(
            args.output_dir, f"swash_{args.varname}_t{tindex:04d}.png"
        )
        plt.savefig(png_path)
        print(f"PNG saved: {png_path}")

    if not args.no_show and not args.makemovie:
        plt.show()

# Save in PDF
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, f"swash_{args.varname}.pdf")
    plt.savefig(pdf_path)
    print(f"PDF saved: {pdf_path}")

# Close file NetCDF
nc.close()
