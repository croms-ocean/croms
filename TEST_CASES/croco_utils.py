# croco_utils.py

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.colors import Normalize
import cartopy.crs as ccrs

# Add other utility functions here as needed.


################################################################
def zlevs(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    Compute the depth of rho or w points for CROCO.

        Parameters:
        h (np.ndarray): Bathymetry (depths of the ocean bottom).
        zeta (np.ndarray): Free surface elevation.
        theta_s (float): S-coordinate surface control parameter.
        theta_b (float): S-coordinate bottom control parameter.
        hc (float): Width (m) of the surface or bottom boundary layer where higher resolution is required.
        N (int): Number of vertical levels.
        type (str): 'r' for rho points, 'w' for w points.
        vtransform (int): Vertical transformation equation (1 or 2).

    Returns:
        np.ndarray: Depths of RHO- or W-points (3D matrix).
    """
    Ndim = np.ndim(h)
    if Ndim == 2:
        M, L = h.shape
    elif Ndim == 1:
        M = len(h)
        L = 1
        h = h[:, np.newaxis]
        zeta = zeta[:, np.newaxis]
    else:
        raise ValueError("zlevs: error - incorrect dimension for h")

    ds = 1.0 / float(N)

    if type == "w":
        sc = ds * (np.arange(0, N + 1) - N)
        Nmax = N + 1
    elif type == "r":
        sc = ds * (np.arange(1, N + 1) - N - 0.5)
        Nmax = N
    else:
        raise ValueError(f"Problem with type = {type}")

    if vtransform == 1:
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5
        )
    elif vtransform == 2:
        Cs = get_csf(sc, theta_s, theta_b)
    else:
        raise ValueError(f"Problem with vtransform = {vtransform}")

    z = np.zeros([Nmax, M, L])

    for k in np.arange(0, Nmax):
        if vtransform == 1:
            cff = hc * (sc[k] - Cs[k])
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / h
            z[k, :, :] = z0 + zeta * (1.0 + z0 * hinv)
        elif vtransform == 2:
            cff = hc * sc[k]
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / (h + hc)
            z[k, :, :] = z0 * h * hinv + zeta * (1.0 + z0 * hinv)

    return z.squeeze()


################################################################
def get_csf(sc, theta_s, theta_b):
    """
    function CSF = get_csf(sc,theta_s,theta_b)
    Get CS-Curves for the new s-coordinate system
    NOTE: Mathematical limits of CSF,csrf for theta_s, theta_b --> 0 match that under "else" logical branches.
    """

    if theta_s > 0.0:
        csrf = (1.0 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - 1.0)
    else:
        csrf = -(sc**2)

    if theta_b > 0.0:
        CSF = (np.exp(theta_b * csrf) - 1.0) / (1.0 - np.exp(-theta_b))
    else:
        CSF = csrf

    return CSF


################################################################
def rho2u_2d(rho_in):
    """
    Convert a 2D field at rho points to a field at u points
    """

    def _r2u(rho_in, Lp):
        u_out = rho_in[:, : Lp - 1]
        u_out += rho_in[:, 1:Lp]
        u_out *= 0.5
        return u_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2u(rho_in, Lshp)


#################################################################
def rho2u_3d(var):
    """
    Convert 3D variable from rho-grid to u-grid.
    """
    return 0.5 * (var[:, :, :-1] + var[:, :, 1:])


################################################################
def rho2v_2d(rho_in):
    """
    Convert a 2D field at rho points to a field at v points
    """

    def _r2v(rho_in, Mp):
        v_out = rho_in[: Mp - 1]
        v_out += rho_in[1:Mp]
        v_out *= 0.5
        return v_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2v(rho_in, Mshp)


################################################################
def u2rho_2d(u_in):
    """
    Convert a 2D field at u points to a field at rho points
    """

    def _uu2ur(uu_in, Mp, Lp):
        L, Lm = Lp - 1, Lp - 2
        u_out = np.zeros((Mp, Lp))
        u_out[:, 1:L] = 0.5 * (u_in[:, 0:Lm] + u_in[:, 1:L])
        u_out[:, 0] = u_out[:, 1]
        u_out[:, L] = u_out[:, Lm]
        return u_out.squeeze()

    assert u_in.ndim == 2, "u_in must be 2d"
    Mp, Lp = u_in.shape
    return _uu2ur(u_in, Mp, Lp + 1)


################################################################
def v2rho_2d(v_in):
    """
    Convert a 2D field at v points to a field at rho points
    """

    def _vv2vr(v_in, Mp, Lp):
        M, Mm = Mp - 1, Mp - 2
        v_out = np.zeros((Mp, Lp))
        v_out[1:M] = 0.5 * (v_in[:Mm] + v_in[1:M])
        v_out[0] = v_out[1]
        v_out[M] = v_out[Mm]
        return v_out.squeeze()

    assert v_in.ndim == 2, "v_in must be 2d"
    Mp, Lp = v_in.shape
    return _vv2vr(v_in, Mp + 1, Lp)


#################################################################
def rho2v_3d(var):
    """
    Convert 3D variable from rho-grid to v-grid.
    """
    return 0.5 * (var[:, :-1, :] + var[:, 1:, :])


################################################################
def rho2uvp(rfield):
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1

    vfield = 0.5 * (rfield[np.arange(0, M), :] + rfield[np.arange(1, Mp), :])
    ufield = 0.5 * (rfield[:, np.arange(0, L)] + rfield[:, np.arange(1, Lp)])
    pfield = 0.5 * (ufield[np.arange(0, M), :] + ufield[np.arange(1, Mp), :])

    return ufield, vfield, pfield


################################################################
def read_latlonmask(fname, var_type):
    """
    Read latitudes and longitudes and mask from CROCO netcdf file.

    Args:
        fname (str): CROCO netcdf file
        var_type (str): variable type ('r' for rho, 'w', 'u', or 'v').

    Returns:
        lat (np.ndarray): Latitude (matrix 2D).
        lon (np.ndarray): Longitude (matrix 2D).
        mask (np.ndarray): Mask (1 for sea, NaN for land).
    """
    # Open NetCDF file
    nc = Dataset(fname)
    lat = nc.variables["lat_rho"][:]
    lon = nc.variables["lon_rho"][:]
    mask = nc.variables.get("mask_rho", np.ones_like(lon))

    # Manage specific types
    if var_type == "u":
        lat = rho2u_2d(lat)
        lon = rho2u_2d(lon)
        mask = mask[:, :-1] * mask[:, 1:]
    elif var_type == "v":
        lat = rho2v_2d(lat)
        lon = rho2v_2d(lon)
        mask = mask[:-1, :] * mask[1:, :]

    # Convert mask : 0 -> NaN
    mask = np.where(mask == 0, np.nan, mask)

    return lat, lon, mask


###################################################
def scoordinate(theta_s, theta_b, N, hc, vtransform):
    """
    function [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)
    Define S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    """
    if vtransform == 2:
        # "NEW_S_COORD"
        ds = 1.0 / N
        sc_r = ds * (np.arange(1, N + 1) - N - 0.5)
        Cs_r = get_csf(sc_r, theta_s, theta_b)
        sc_w = ds * (np.arange(0, N + 1) - N)
        Cs_w = get_csf(sc_w, theta_s, theta_b)

    else:
        # "OLD_S_COORD"

        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        sc_w = (np.arange(0, N + 1) - N) / N
        Cs_w = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc_w) + theta_b * (
            cff2 * np.tanh(theta_s * (sc_w + 0.5)) - 0.5
        )
        sc_r = (np.arange(1, N + 1) - N - 0.5) / N
        Cs_r = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc_r) + theta_b * (
            cff2 * np.tanh(theta_s * (sc_r + 0.5)) - 0.5
        )

    return sc_r, Cs_r, sc_w, Cs_w


###################################################
def rempoints(var, npts):
    """
    Reduce a 2D matrix by removing border points

    Args:
        var (np.ndarray): input 2D matrix
        npts (list[int]): number of points to remove at each border
                          with format [west, east, south, north].

    Returns:
        np.ndarray: 2D matrix after removing border points
    """
    # Initial dimensions
    M, L = var.shape

    # Extract index from input
    west, east, south, north = npts

    # Reduce dimensions from input number of point to remove
    var = var[south : M - north, west : L - east]

    return var


#################################################################
def get_type(fname, vname, vlevin):
    """
        Determines the "type" of a CROCO variable: rho, u, or v.

        Args:
        fname (str): Name of the NetCDF file.
        vname (str): Name of the variable to examine.
        vlevin (int): Input vertical level.

    Returns:
        tuple:
            - type (str): Type of the variable ('r', 'u', 'v', 'w', or '').
            - vlevout (int): Output vertical level (unchanged or modified).
    """
    vlevout = vlevin
    var_type = "r"  # By default, type is 'r'

    try:
        # NetCDF Open NetCDF file
        with Dataset(fname, "r") as nc:
            if vname not in nc.variables:
                return "", vlevout

            var = nc.variables[vname]
            dims = var.dimensions
            ndim = len(dims)

            if ndim == 1:  # Case of 1D variable
                return "", vlevout

            i = 0
            name = dims[i]

            # Check if temporal dimension
            if not (name.endswith("time") or name.startswith("time")):
                print("Avertissement : Pas dépendant du temps.")
            else:
                i += 1

            # Check if vertical dimension
            name = dims[i]
            if name.startswith("s"):
                if name.endswith("w"):
                    return "w", vlevout
                else:
                    i += 1
            else:
                vlevout = 0

            # Check if 2D-H dimension : lat/y or lon/x
            name = dims[i]
            if not (name.startswith("e") or name.startswith("y")):
                return "", vlevout
            else:
                if name.endswith("v"):
                    return "v", vlevout
                if name.endswith("u"):
                    return "u", vlevout

            # Check second variable on 2D-H
            if i + 1 < ndim:
                name = dims[i + 1]
                if name.startswith("x"):
                    if name.endswith("u"):
                        return "u", vlevout
                    if name.endswith("v"):
                        return "v", vlevout

    except Exception as e:
        print(f"Error while reading the NetCDF file : {e}")
        return "", vlevout

    return var_type, vlevout


#################################################################
def vinterp(var, z, depth):
    """
    Interpolates a 3D variable onto a horizontal level of constant depth.

    Args:
        var (numpy.ndarray): Variable to process (3D matrix of size [N, Mp, Lp]).
        z (numpy.ndarray): Depths (m) of RHO or W points (3D matrix of size [N, Mp, Lp]).
        depth (float): Depth of the slice (scalar; meters, negative value).

    Returns:
        numpy.ndarray: Interpolated horizontal slice (2D matrix of size [Mp, Lp]).
    """

    N, Mp, Lp = z.shape

    # Find the position of the nearest vertical levels
    below = z < depth
    levs = np.sum(below, axis=0)  # Number of levels below the specified depth
    levs[levs == N] = N - 1  # Limit to the second-to-last level to avoid overflows

    # Create a mask to identify valid indices
    mask = np.where(levs > 0, 1, np.nan)

    # Generate linear positions for quick access
    imat, jmat = np.meshgrid(np.arange(Lp), np.arange(Mp), indexing="ij")
    pos = N * Mp * imat + N * jmat + levs
    pos = pos.astype(int)  # Assurer que les indices sont des entiers
    pos[np.isnan(mask)] = 0  # Éviter les erreurs sur les valeurs NaN

    # Do the interpolation
    z1 = z.ravel()[pos + 1]
    z2 = z.ravel()[pos]
    v1 = var.ravel()[pos + 1]
    v2 = var.ravel()[pos]

    #
    vnew = mask * (((v1 - v2) * depth + v2 * z1 - v1 * z2) / (z1 - z2))

    return vnew.reshape(Mp, Lp)


#################################################################
def get_depths(fname, gname, tindex, point_type):
    """
    Calculates the depths of sigma levels.

    Args:
        fname (str): Path to the history file.
        gname (str): Path to the grid file.
        tindex (int): Time index.
        point_type (str): Point type ('rho', 'u', 'v', 'w').

    Returns:
        numpy.ndarray: Depths of sigma levels (3D matrix).
    """

    # Reading grid data
    with Dataset(gname, mode="r") as nc:
        h = nc.variables["h"][:]

    # Reading variables
    with Dataset(fname, mode="r") as nc:
        zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
        hmorph = nc.variables.get("hmorph", None)
        if hmorph is not None:
            h = np.squeeze(hmorph)

        # Reading sigma grid parameter
        try:
            theta_s = nc.variables["theta_s"][:]
            theta_b = nc.variables["theta_b"][:]
            Tcline = nc.variables["Tcline"][:]
            hc = nc.variables.get("hc", None)
        except KeyError:
            theta_s = nc.theta_s
            theta_b = nc.theta_b
            Tcline = nc.Tcline
            hc = getattr(nc, "hc", None)

        if Tcline is not None:
            hmin = np.min(h)
            hc = min(hmin, Tcline)

        N = len(nc.dimensions["s_rho"])
        VertCoordType = getattr(nc, "VertCoordType", "")
        vtrans = np.squeeze(nc.variables.get("Vtransform", None))
        s_coord = 2 if VertCoordType == "NEW" or vtrans == 2 else 1

        if s_coord == 2:
            hc = Tcline

    if zeta is None:
        zeta = np.zeros_like(h)

    # Calculating depths of sigma levels
    vtype = point_type
    if point_type in ("u", "v"):
        vtype = "r"

    z = zlevs(h, zeta, theta_s, theta_b, hc, N, vtype, vtrans)

    # Swicth to U- or V- grid
    if point_type == "u":
        z = rho2u_3d(z)
    elif point_type == "v":
        z = rho2v_3d(z)

    return z


#################################################################
def get_hslice(fname, gname, vname, tindex, level, var_type):
    """
    Extracts a horizontal slice of a CROCO variable.

    Args:
        fname (str): Path to the CROCO file (average or history).
        gname (str): Path to the CROCO grid file.
        vname (str): Name of the variable.
        tindex (int): Time index.
        level (float): Vertical level for the slice:
                       - level = 0: 2D horizontal variable (e.g., zeta)
                       - level > 0: slice at the sigma level (integer >= 1)
                       - level < 0: interpolation at z = level.
        var_type (str): Type of the variable ('r', 'u', 'v', or 'w').

    Returns:
        numpy.ndarray: Horizontal slice (2D matrix).
    """

    with Dataset(fname, mode="r") as nc:
        if level == 0:
            # Variable 2D
            var = np.squeeze(nc.variables[vname][tindex, :, :])
            var[var == 0] = np.nan

            # Correction for wetting/drying
            with Dataset(gname, mode="r") as ng:
                h = ng.variables["h"][:]
            hmorph = nc.variables.get("hmorph", None)
            if hmorph is not None:
                h = np.squeeze(hmorph[tindex, :, :])
            zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
            D = zeta + h

            Dcrit = nc.variables.get("Dcrit", None)
            if Dcrit is not None:
                Dcrit = Dcrit[:] + 1e-5
            else:
                Dcrit = 0.2 + 1e-5

            if var_type == "u":
                D = rho2u_2d(D)
            elif var_type == "v":
                D = rho2v_2d(D)

            var[D <= Dcrit] = np.nan

        elif level > 0:
            # Slice at a specific sigma level
            var = np.squeeze(nc.variables[vname][tindex, level - 1, :, :])
            var[var == 0] = np.nan

        else:
            # Horizontal slice at a given depth level
            # Get the depths of sigma levels
            z = get_depths(fname, gname, tindex, var_type)

            # Read the 3D variable and interpolate
            var_sigma = np.squeeze(nc.variables[vname][tindex, :, :, :])
            var = vinterp(var_sigma, z, level)

    return var


#################################################################
def get_var(hisfile, gridfile, vname, tindex, vlevel, coef, rempts):
    """
    Retrieves a variable from a CROCO NetCDF file.

    Args:
        hisfile (str): Name of the NetCDF history file.
        gridfile (str): Name of the NetCDF grid file.
        vname (str): Name of the variable to retrieve.
        tindex (int): Time index of the variable.
        vlevel (float): Vertical level.
        coef (float): Multiplicative coefficient.
        rempts (int): Points to remove at the boundaries.

    Returns:
        tuple: (lat, lon, mask, var) where:
            lat (np.ndarray): Latitudes.
            lon (np.ndarray): Longitudes.
            mask (np.ndarray): Mask (1 = sea, NaN = land).
            var (np.ndarray): Requested variable.
    """

    print(f"-- {vname} --")

    if not gridfile:
        gridfile = hisfile
    if vname in ["h", "pm", "pn"]:
        lat, lon, mask = read_latlonmask(gridfile, "r")
        with Dataset(gridfile, "r") as nc:
            var = nc.variables[vname][:]

        mask = np.where(mask == 0, np.nan, mask)
        var = mask * var

    else:
        var_type, vlev = get_type(hisfile, vname, vlevel)
        lat, lon, mask = read_latlonmask(gridfile, var_type)
        var = coef * mask * get_hslice(hisfile, gridfile, vname, tindex, vlev, var_type)

    return lat, lon, mask, var


#################################################################
def do_plot(lon, lat, var, colmin, colmax, ncol, hisfile, gridfile, tindex, vlevel):
    """
    Do the plot
    """
    # Pseudocolor plot
    plt.figure()

    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(
        lon,
        lat,
        var,
        shading="flat",
        cmap="jet",
        norm=Normalize(vmin=colmin, vmax=colmax),
    )
    plt.colorbar(mesh, ax=ax)
    plt.show()


#################################################################
def get_vertical_section(fname, vname, tindex, direction, ind_sec, gname):
    """
    Extracts a vertical section of a 3D variable.

    Args:
        fname (str): Name of the NetCDF history file.
        vname (str): Name of the variable (e.g., "temp").
        tindex (int): Time index.
        direction (str): Section direction ("x" for zonal or "y" for meridional).
        ind_sec (int): Index of the row or column for the section.
        gname (str): Name of the NetCDF grid file.

    Returns:
        dict: Contains depths,
              distance (number of grid points),
              and interpolated values.
    """

    with Dataset(fname, "r") as nc:
        topo = np.squeeze(nc.variables["h"][:, :])
        # Déterminer le type de grille pour obtenir les profondeurs
        if vname == "w":
            point_type = "w"
        elif vname == "u":
            point_type = "u"
            topo = rho2u_2d(topo)
        elif vname == "v":
            point_type = "v"
            topo = rho2v_2d(topo)
        else:
            point_type = "r"

        # Get depths
        z = get_depths(fname, gname, tindex, point_type)
        temp = np.squeeze(nc.variables[vname][tindex, :, :, :])

    if direction == "y":
        section_data = temp[:, :, ind_sec]
        depth_section = z[:, :, ind_sec]
        distance = np.arange(section_data.shape[1])  # Distance in grid points
        topo_section = -topo[:, ind_sec]
    elif direction == "x":
        section_data = temp[:, ind_sec, :]
        depth_section = z[:, ind_sec, :]
        distance = np.arange(section_data.shape[1])  # Distance in grid points
        topo_section = -topo[ind_sec, :]
    else:
        raise ValueError("Direction invalide. Choisissez 'x' ou 'y'.")

    return {
        "topo": topo_section,
        "depth": depth_section,
        "distance": distance,
        "variable": section_data,
    }


#################################################################
def plot_vertical_section(section_data, title="Vertical Section", unit="°C"):
    """
    Plots a vertical section.

    Args:
        section_data (dict): Contains depths, distance, and values.
        title (str): Title of the plot.
        unit (str): Unit of the variable.
    """

    depth = section_data["depth"]
    topo = section_data["topo"]
    distance = section_data["distance"]
    variable = section_data["variable"]
    variable[variable == 0] = np.nan

    plt.figure(figsize=(12, 6))
    plt.pcolormesh(distance, depth, variable, cmap="viridis", shading="auto")
    plt.plot(
        distance,
        topo,
        color="black",
        linestyle="--",
        linewidth=1.5,
    )
    plt.colorbar(label=unit)
    plt.title(title)
    plt.xlabel("Distance (grille)")
    plt.ylabel("Profondeur (m)")
    plt.grid(True)
    # plt.show()
