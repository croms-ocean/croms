#!/bin/env python3
##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import shutil
import pytest

# numpy
import numpy

# netcdf
from netCDF4 import Dataset

# internal
from libs.check import compare_netcdf_files


##########################################################
def helper_create_netcdf_file(file: str) -> None:
    dataset = Dataset(file, "w", format="NETCDF4")

    # dimensions
    dataset.createDimension("time", None)
    dataset.createDimension("lat", 10)
    dataset.createDimension("lon", 10)

    # create some vars
    lats = dataset.createVariable("lat", "f4", ("lat",))
    lons = dataset.createVariable("lon", "f4", ("lon",))
    value = dataset.createVariable(
        "value",
        "f4",
        (
            "time",
            "lat",
            "lon",
        ),
    )

    # fill values
    lats[:] = numpy.arange(40.0, 50.0, 1.0)
    lons[:] = numpy.arange(-110.0, -100.0, 1.0)
    value[0, :, :] = numpy.random.uniform(
        0, 100, size=(10, 10)
    )  # unifrom random values
    value[1, :, :] = numpy.random.uniform(
        0, 100, size=(10, 10)
    )  # unifrom random values
    value[2, :, :] = numpy.random.uniform(
        0, 100, size=(10, 10)
    )  # unifrom random values

    # close
    dataset.close()


##########################################################
def help_patch_netcdf_file(file: str, variable_name: str, make_close: bool) -> None:
    # load & copy out the given variable
    dataset = Dataset(file, "r", format="NETCDF4")
    out = dataset.variables[variable_name][:]
    dataset.close()

    # change one value inside
    if make_close:
        out[1][5][5] += 0.00001
    else:
        out[1][5][5] = 3.869

    # write again (need to re-open to get it working)
    dataset = Dataset(file, "r+", format="NETCDF4")
    dataset.variables[variable_name][:] = out
    dataset.close()


##########################################################
def test_compare_ok(tmp_path):
    # create file
    helper_create_netcdf_file(f"{tmp_path}/test_check.nc")

    # make some copies
    shutil.copyfile(f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_ok.nc")

    # compare both files are same
    compare_netcdf_files(
        f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_ok.nc"
    )


##########################################################
def test_compare_not_ok(tmp_path):
    # create file
    helper_create_netcdf_file(f"{tmp_path}/test_check.nc")

    # copy & patch
    shutil.copyfile(
        f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_not_ok.nc"
    )
    help_patch_netcdf_file(f"{tmp_path}/test_check_copy_not_ok.nc", "value", False)

    # compare
    with pytest.raises(Exception, match="Non strict equality in variable 'value'"):
        compare_netcdf_files(
            f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_not_ok.nc"
        )


##########################################################
def test_compare_not_ok_close(tmp_path):
    # create file
    helper_create_netcdf_file(f"{tmp_path}/test_check.nc")

    # copy & patch
    shutil.copyfile(
        f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_not_ok.nc"
    )
    help_patch_netcdf_file(f"{tmp_path}/test_check_copy_not_ok.nc", "value", True)

    # compare
    with pytest.raises(Exception, match="Non strict equality in variable 'value'"):
        compare_netcdf_files(
            f"{tmp_path}/test_check.nc", f"{tmp_path}/test_check_copy_not_ok.nc"
        )
