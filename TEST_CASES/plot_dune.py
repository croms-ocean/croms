#!/usr/bin/env python3

import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
import xarray
import xarray.plot


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
    data = xarray.open_dataset(args.file)
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)


lastonly = True
# plot fraction of sand in sediment for dune test case
# (migration of a dune under stationnary current)
# if lastonly is TRUE then only the last time step of
# the netcdf file is plotted; if lastonly is FALSE then
# all time step are plotted

# Check if variables exist in dataset
if (
    "hmorph" in data
    and "HSED" in data
    and "DZS" in data
    and "GRAV_sed" in data
    and "SAND_sed" in data
):
    bathy = data.hmorph
    hsed = data.HSED
    dzstot = data.DZS
    grav_sed = data.GRAV_sed
    sand_sed = data.SAND_sed
    is_mustang = True
    title = "DUNE - MUSTANG"

elif (
    "hmorph" in data
    and "bed_thick" in data
    and "bed_frac_sand_01" in data
    and "bed_frac_sand_02" in data
):
    bathy = data.hmorph
    hsed = data.bed_thick.sum(axis=1)
    dzstot = data.bed_thick
    grav_sed = data.bed_frac_sand_01
    sand_sed = data.bed_frac_sand_02
    is_mustang = False
    title = "DUNE - USGS"

else:
    print("Error: Required variables not found in the NetCDF file")
    exit(1)

hmin = 7.0
hsed0 = 3.0
dx = 2.0
if lastonly:
    itlist = [len(data.time) - 1]
else:
    itlist = range(len(data.time))
for it in itlist:
    depthsinit = hmin - bathy[0, 1, :] - hsed0
    depthst = hmin - bathy[it, 1, :] - hsed0

    if is_mustang:
        test = sand_sed[it, :, 1, :] / (grav_sed[it, :, 1, :] + sand_sed[it, :, 1, :])
        dzs = dzstot[it, :, 1, :]
    else:
        test = sand_sed[it, :, 1, :]
        dzs = dzstot[it, :, 1, :]

    nz = test.shape[0]
    d = data.xi_rho * dx
    colnorm = matplotlib.colors.Normalize(vmin=0.2, vmax=0.8)
    varcolor = matplotlib.cm.ScalarMappable(norm=colnorm, cmap="jet")  # ,cmap='YlOrBr')
    varcolor.set_array(test)
    col = varcolor.to_rgba(test)

    fig, ax = plt.subplots(figsize=(10, 4))
    if is_mustang:
        offset_dzs = hmin - bathy[it, 1, :] - hsed[it, 1, :] - hsed0
        for i in range(nz):
            if i > 0:
                offset_dzs = offset_dzs + dzs[i - 1, :]
            p1 = ax.bar(
                d,
                dzs[i, :],
                bottom=offset_dzs,
                width=dx,
                linewidth=0.5,
                edgecolor="black",
                color=col[i, :, :],
            )
        ax.set_title(
            "%s\nsand fraction after %s"
            % (title, data.indexes["time"][it] - data.indexes["time"][0])
        )
    else:
        offset_dzs = hmin - bathy[it, 1, :] - hsed0
        for i in range(nz):
            offset_dzs = offset_dzs - dzs[i, :]
            p1 = ax.bar(
                d,
                dzs[i, :],
                bottom=offset_dzs,
                width=dx,
                linewidth=0.5,
                edgecolor="black",
                color=col[i, :, :],
            )
        ax.set_title(
            "%s\nsand fraction after %s days"
            % (
                title,
                (data.indexes["time"][it] - data.indexes["time"][0]) / (3600.0 * 24.0),
            )
        )
    ax.set_xlabel("distance (m)")
    ax.set_ylabel("height (m)")
    p2 = ax.plot(d, depthsinit, color="red", linewidth=2, linestyle="dashed")
    p3 = ax.plot(d, depthst, color="grey", linewidth=1)
    fig.colorbar(varcolor, ax=ax, extend="both")


# Save outputs
if args.makepng:
    png_path = os.path.join(args.output_dir, "dune.png")
    fig.savefig(png_path, dpi=300)
    print(f"PNG saved: {png_path}")

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "dune.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved: {pdf_path}")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()
