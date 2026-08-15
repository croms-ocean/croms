#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse


# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot data from a NetCDF file and optionally save as PDF or PNG.",
    epilog="Example usage:\n"
    + "  python script.py --makepdf --file croco_abl_his.nc\n"
    + "  python script.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="croco_abl_his.nc",
    help="Path to the NetCDF file (default: croco_abl_his.nc)",
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

# Read data from NetCDF file
try:
    nc = Dataset(args.file)
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

time = nc.variables["scrum_time"][tndx] / 86400
t_abl = nc.variables["t_abl"][tndx, :, 1, :] - 273.15

plt.figure(figsize=(6, 4))
plt.subplot(1, 1, 1)
contour = plt.contourf(t_abl, np.arange(14, 22.1, 0.1), cmap="viridis")
plt.colorbar(contour)
plt.clim(14, 22)
plt.xlabel("X [-]")
plt.ylabel("Z [-]")
plt.title(f"KILPATRICK - Air temperature [°C] at {time:.2f} days")

# Determine output file paths
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(output_dir, "kilpatrick_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(output_dir, "kilpatrick_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Show plots if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --show to enable).")
