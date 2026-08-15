#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import argparse

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the ISOLITON test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="isoliton_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex", type=int, nargs="+", default=[30, 50, 70], help="Time indices to plot"
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as a PDF file")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# User defined parameters
vname = "zeta"
tindex = 10

# Caxis settings based on variable name
caxis_params = {
    "te": ("r", 1, 17, 0.1, 22, 1),
    "ze": ("r", 0, -100, 5, 100, 100),
    "ub": ("u", 0, -50, 5, 50, 100),
    "vb": ("v", 0, -50, 5, 50, 100),
}

var_type, ddd, cmin, dc, cmax, cff = caxis_params.get(
    vname[:2], ("r", 1, -100, 20, 100, 100)
)

# Load parent dataset
parent_ds = xr.open_dataset(args.file)
N = parent_ds.dims["s_rho"]
time = round(parent_ds["scrum_time"][tindex].item() / (24 * 3600))
print(f"Day: {time}")
X = 1e-3 * parent_ds["x_rho"]
Y = 1e-3 * parent_ds["y_rho"]
t1 = cff * parent_ds[vname].isel(time=tindex)
parent_ds.close()

# Load child dataset if available
childhis = f"{args.file}.1"
nestvortex = os.path.exists(childhis)

if nestvortex:
    child_ds = xr.open_dataset(childhis)
    X2 = 1e-3 * child_ds["x_rho"]
    Y2 = 1e-3 * child_ds["y_rho"]
    t = cff * child_ds[vname].isel(time=tindex)
    child_ds.close()

# Plot parent data
plt.figure()
plt.contourf(X, Y, t1, levels=np.arange(cmin, cmax + dc, dc), cmap="jet")
plt.colorbar()
plt.contour(
    X,
    Y,
    t1,
    levels=np.arange(cmin, cmax + dc, dc),
    colors="grey",
    linestyles="dashed",
    linewidths=0.2,
    colorbar=False,
)

plt.xlabel("X [km]")
plt.ylabel("Y [km]")
plt.title(f"{vname} - Day {time}")


# Plot child data if available
if nestvortex:
    plt.contour(
        X2,
        Y2,
        t,
        levels=np.arange(cmin, cmax + dc, dc),
        colors="k",
        linestyles="solid",
        linewidths=0.3,
        colorbar=False,
    )

# Save in PNG
if args.makepng:
    png_path = os.path.join(args.output_dir, "vortex.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Save in PDF
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "vortex.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

# Show or not
if not args.no_show:
    plt.show()
else:
    plt.close()

# Plot diff parent - child data if available
if nestvortex:
    # Interpolation using xarray's built-in method
    # Now create the DataArray with sorted coordinates
    parent = xr.DataArray(
        data=t1.values,
        dims=("eta_rho", "xi_rho"),
        coords={"xi_rho": X[0, :], "eta_rho": Y[:, 0]},
    )
    child = xr.DataArray(
        data=t.values,
        dims=("eta_rho", "xi_rho"),
        coords={"xi_rho": X2[0, :], "eta_rho": Y2[:, 0]},
    )

    # Interpolate parent onto the grid of X2, Y2
    parentint = parent.interp(xi_rho=X2[0, :], eta_rho=Y2[:, 0], method="linear")

    # Calculate difference
    tdiff = parentint.values - child.values
    max_diff = np.max(np.abs(tdiff))
    rel_diff = 100 * max_diff / np.max(np.abs(parentint))
    print(f"Max difference = {max_diff:.2f}")
    print(f"Relative difference = {rel_diff:.2f} %")

    plt.figure()
    # plt.imshow(np.flipud(tdiff), cmap='coolwarm', aspect='auto')
    plt.contourf(Y2, X2, tdiff, cmap="coolwarm", aspect="auto")
    plt.xlabel("X [km]")
    plt.ylabel("Y [km]")
    plt.axis([X.min(), X.max(), Y.min(), Y.max()])
    plt.colorbar()
    plt.title(f"Parent - Child (cm) : {vname} - Day {time}")
    plt.show()

    # Save in PNG
    if args.makepng:
        png_path = os.path.join(args.output_dir, "vortex_diff.png")
        plt.savefig(png_path, dpi=300)
        print(f"PNG file '{png_path}' has been created.")

    # Save in PDF
    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "vortex_diff.pdf")
        plt.savefig(pdf_path, transparent=True)
        print(f"PDF file '{pdf_path}' has been created.")

    # Show or not
    if not args.no_show:
        plt.show()
    else:
        plt.close()
