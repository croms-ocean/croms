#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the RIP test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="rip_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--makemovie", action="store_true", help="Generate an MP4 movie (frame-by-frame)"
)
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument(
    "--makepdf", action="store_true", help="Save the final plot as a PDF file"
)
parser.add_argument(
    "--pltvort", action="store_true", help="Plot vorticity instead of speed"
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

# Read parameters
tindex_last = len(nc.variables["scrum_time"][:]) - 1
tstart = 0 if args.makemovie else tindex_last
tend = tindex_last

h = nc.variables["h"][:]
x = nc.variables["x_rho"][:]
y = nc.variables["y_rho"][:]
pm = nc.variables["pm"][:]
pn = nc.variables["pn"][:]
N = len(nc.dimensions["s_rho"])
Dcrit = nc.variables["Dcrit"][:]

# Initialize plot
fig, ax = plt.subplots(figsize=(10, 8))

# Loop over time steps
for tindex in range(tstart, tend + 1):
    print(f"Tindex = {tindex}")

    time = nc.variables["scrum_time"][tindex] / 86400  # Convert time to days
    zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
    u = np.squeeze(nc.variables["u"][tindex, -1, :, :])
    v = np.squeeze(nc.variables["v"][tindex, -1, :, :])
    ub = np.squeeze(
        nc.variables["ubar"][tindex, :, :] + nc.variables["ust2d"][tindex, :, :]
    )
    vb = np.squeeze(
        nc.variables["vbar"][tindex, :, :] + nc.variables["vst2d"][tindex, :, :]
    )

    # Mask based on minimum depth
    mask = np.ones_like(zeta)
    mask[(h + zeta) <= Dcrit] = np.nan
    ur = cr.u2rho_2d(ub)
    vr = cr.v2rho_2d(vb)

    # Plot vorticity or speed
    if args.pltvort:
        vort = mask * cr.psi2rho(cr.vorticity(u, v, pm, pn))
    else:
        speed = mask * ur

    # Plot
    ax.clear()
    cmin, cmax = (-0.01, 0.01) if args.pltvort else (-0.25, 0.25)
    cmap = plt.cm.cool

    if args.pltvort:
        cont = ax.contourf(
            x, y, vort, levels=np.linspace(cmin, cmax, 10), cmap=cmap, extend="both"
        )
    else:
        cont = ax.contourf(
            x, y, speed, levels=np.linspace(cmin, cmax, 10), cmap=cmap, extend="both"
        )

    plt.colorbar(cont, ax=ax)
    var_I = ~(np.sqrt(ur**2 + vr**2) < 0.05)
    ax.quiver(x[var_I], y[var_I], ur[var_I], vr[var_I])
    ax.set_xlim([100, 760])
    ax.set_ylim([10, 760])
    ax.set_title(f"Rip Test Case - time = {int(time * 24)} h", fontsize=14)
    ax.set_xlabel("Cross-shore distance [m]", fontsize=14)
    ax.set_ylabel("Along-shore distance [m]", fontsize=14)
    ax.grid()
    ax.set_aspect("equal")

    # Save outputs
    if args.makemovie:
        plt.savefig(os.path.join(args.output_dir, f"frame_{tindex:04d}.png"))
    if args.makepng:
        png_path = os.path.join(args.output_dir, f"rip_plot_t{tindex:04d}.png")
        plt.savefig(png_path)
        print(f"PNG saved to {png_path}")

    # Display or suppress the plot
    if not args.no_show and not args.makemovie:
        plt.show()
    elif args.makemovie:
        plt.close()

# Save the final plot as a PDF if requested
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "rip_final_plot.pdf")
    plt.savefig(pdf_path)
    print(f"PDF saved to {pdf_path}")

# Close NetCDF file
nc.close()
