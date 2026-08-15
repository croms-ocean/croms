#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SED_TOY_CONSOLID test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="sed_toy_consolid_his.nc",
    help="Path to the NetCDF file",
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


idy, idx = 1, 1
depth = 20  # 20m water column depth
beddepth = 0.041  # 4cm sediment depth
srho = 2650  # rho sediment classes
nl = 41  # 41 sediment layers
a = 1  # Equilibrium bulk critical stress for erosion profile
slope = 2
offset = 3.4
masslayer0 = 0
days = [0.5, 2, 32, 38]
tplot = [6, 24, 384, 456]  # record to plot
ntplot = len(tplot)


h = nc.variables["h"][idy, idx]
bostr = nc.variables["bostr"][:, idy, idx]
Timeindays = nc.variables["scrum_time"][:] / 86400  # Convert to days
nt = len(Timeindays)

zbed = np.zeros((nl, 4))
taucb = np.zeros((nl, 4))
bedtaucrit = np.zeros((nl, 4))

for it0 in range(ntplot):
    itime = tplot[it0]
    bedthick = nc.variables["bed_thick"][itime, :, idy, idx]
    fracs1 = nc.variables["bed_frac_sand_01"][itime, :, idy, idx]
    fracs2 = nc.variables["bed_frac_sand_02"][itime, :, idy, idx]
    fracm1 = nc.variables["bed_frac_mud_01"][itime, :, idy, idx]
    fracm2 = nc.variables["bed_frac_mud_02"][itime, :, idy, idx]

    # Instantaneous profile of bulk critical stress for erosion
    bedtaucrit[:, it0] = nc.variables["bed_tau_crit"][itime, :, idy, idx]

    # Equilibrium bulk critical stress for erosion profile
    masslayer = srho * (
        fracs1 * bedthick + fracs2 * bedthick + fracm1 * bedthick + fracm2 * bedthick
    )
    massdepth = np.hstack(([masslayer0], np.cumsum(masslayer[1:])))
    massdepth = np.maximum(massdepth, 1e-10)
    taucb[:, it0] = a * np.exp((np.log(massdepth) - offset) / slope)

    zbed0 = np.hstack(([-beddepth], -beddepth + np.cumsum(np.flip(bedthick[1:]))))
    zbed[:, it0] = np.flip(zbed0) * 100  # Convert to cm

nc.close()

# ======================================================================
# Plot

fig, axes = plt.subplots(ncols=4, nrows=2, figsize=(10, 7))
gs = axes[1, 1].get_gridspec()
# remove the underlying Axes
for ax in axes[0, :]:
    ax.remove()
axbig = fig.add_subplot(gs[0, :])

axbig.plot(Timeindays, bostr, linewidth=1)
axbig.set_title("Sequence of depth-limited erosion, deposition, and compaction")
axbig.set_ylabel("Bottom stress (Pa)", fontsize=12)
axbig.set_xlabel("Days")
axbig.grid(True, linewidth=0.5, color="lightgray")
axbig.plot(Timeindays[tplot], bostr[tplot], "r*")

for it0 in range(ntplot):
    ax = axes[1, it0]
    ax.plot(
        bedtaucrit[:, it0], zbed[:, it0], label="Critical stress for erosion (N/m²)"
    )
    ax.plot(
        taucb[:, it0], zbed[:, it0], label="Equilibrium bulk critical stress (N/m²)"
    )
    ax.set_xlim([0, 1.5])
    ax.set_ylim([-2.5, 0])
    if it0 == 0:
        ax.set_ylabel(
            "Z(cm)",
            fontsize=12,
        )
        ax.legend(
            bbox_to_anchor=(0.5, 0),
            loc="lower center",
            bbox_transform=fig.transFigure,
            ncol=2,
        )
    for iz in range(len(zbed[:, it0])):
        ax.axhline(y=zbed[:, it0][iz], linestyle="--", color="gray")

    ax.set_xlabel(f"{days[it0]} days")
    ax.grid(True, linewidth=0.5, color="lightgray")


# Save outputs
if args.makepng:
    png_path = os.path.join(args.output_dir, "sed_toy_consolid.png")
    fig.savefig(png_path, dpi=300)
    print(f"PNG saved: {png_path}")

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "sed_toy_consolid.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved: {pdf_path}")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()
