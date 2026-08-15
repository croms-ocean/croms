#!/bin/env python3
##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# imports
from netCDF4 import Dataset
import math
import numpy
from .helpers import Messaging


##########################################################
class CompareErrorLogger:
    def __init__(self, max_stored: int = 10, max_total: int = 50):
        self.count = 0
        self.var_error_logs = {}
        self.var_error_count = {}
        self.max_stored = max_stored
        self.max_total = max_total

    def append(
        self,
        varname: str,
        ref,
        actual,
        coord: list,
        is_strict_compare: bool = True,
    ):
        # to help readability with right operator displayed
        if is_strict_compare:
            compare_name = "strict"
            operator = "!="
        else:
            compare_name = "close"
            operator = "!~="

        # calc diff
        diff = abs(ref - actual)

        # build message
        self.append_raw(
            varname,
            "Non %s equality in variable '%s' at %s: ref %s actual : %s %s %s (diff=%s)"
            % (
                compare_name,
                varname,
                {",".join(coord)},
                operator,
                ref,
                operator,
                actual,
                diff,
            ),
        )

    def append_raw(self, varname: str, message: str):
        # count it
        self.count += 1

        # first seen for this var
        if varname not in self.var_error_logs:
            self.var_error_logs[varname] = []
            self.var_error_count[varname] = 0

        # count & get log
        self.var_error_count[varname] += 1
        var_error_log = self.var_error_logs[varname]

        # extract some meaning on limits
        can_still_log_var = len(var_error_log) <= self.max_stored

        # if can still log
        if can_still_log_var:
            var_error_log.append(message)

    def has_error(self) -> bool:
        return self.count > 0

    def __str__(self):
        # prepare some vars
        error_var_names = ", ".join(self.var_error_logs.keys())
        total_error_count = self.count

        # prep per variable messages
        var_messages = []
        log_count = 0
        for varname, log in self.var_error_logs.items():
            var_messages.append(f"-------------------- {varname} ---------------------")
            for entry in log:
                log_count += 1
                if log_count < self.max_total:
                    var_messages.append(entry)
                elif log_count == self.max_total:
                    var_messages.append(
                        ".................... too many errors "
                        + f"(>{self.max_total}), stop logging details ............."
                    )

        # assemble details
        details = "\n".join(var_messages)

        # build full message
        message = f"Found {total_error_count} errors in : {error_var_names}\n{details}"

        # ok
        return message


##########################################################
def recurse_compare_current_dim(
    error_log: CompareErrorLogger,
    varname: str,
    shape_ref: tuple,
    shape_actual: tuple,
    cusor_ref,
    cursor_actual,
    dim_id: int,
    current_dims: list = [],
) -> None:
    """
    Recursively compare the values of each dimensions of the mesh.

    It is usefull mostly if we want to see which value is incorrect opposite to
    numpy who currently say only what is the mistake. It also offer a different
    approach to crosscheck we are right with the way of using numpy on top of
    netcdf.

    Parameters:
    -----------
    varname: str
        Name of the variable we are looping over.
    shape_ref: tuple
        Size of each dimensions in ref file.
    shape_actual: tuple
        Size of each dimensions in ref file.
    cursor_ref:
        Dimension we are looping over. At first call it is a netCDF4.Variable,
        after if is an array in principle.
    cursor_actual:
        Dimension we are looping over. At first call it is a netCDF4.Variable,
        after if is an array in principle.
    dim_id: int
        The current dimension we are looping in.
    """
    # check
    if len(shape_ref) != len(shape_actual):
        f"Shape is not same in both file for '{varname}'"

    # nothing to check
    if len(shape_ref) == 0:
        return

    # check
    assert dim_id < len(shape_ref)

    # check size of the dim
    if shape_actual[dim_id] != shape_ref[dim_id]:
        Messaging.step_error(f"Shape is not same in both file for '{varname}'")
        raise Exception(f"Shape is not same in both file for '{varname}'")

    # is last level
    if dim_id == len(shape_ref) - 1:
        # loop
        for i in range(shape_actual[dim_id]):
            if not math.isclose(cusor_ref[i], cursor_actual[i]):
                error_log.append(
                    varname,
                    cusor_ref[i],
                    cursor_actual[i],
                    current_dims + [str(i)],
                )
            if not cusor_ref[i] == cursor_actual[i]:
                error_log.append(
                    varname,
                    cusor_ref[i],
                    cursor_actual[i],
                    current_dims + [str(i)],
                )
    else:
        for i in range(shape_actual[dim_id]):
            recurse_compare_current_dim(
                error_log,
                varname,
                shape_ref,
                shape_actual,
                cusor_ref[i],
                cursor_actual[i],
                dim_id + 1,
                current_dims=current_dims + [str(i)],
            )


def compare_netcdf_variables(
    ref: Dataset, actual: Dataset, skipped=["spherical", "hc", "time_step"]
) -> None:
    """
    Compare variables between two NetCDF datasets, checking for shape consistency
    and value equality. Logs detailed errors if discrepancies are found.

    Args:
        ref: Reference NetCDF dataset.
        actual: Actual NetCDF dataset to compare with.
        skipped: List of variable names to skip during comparison (default is ["spherical", "hc"]).
    """
    for var in ref.variables.keys():
        if var in skipped:
            continue
        # Extract and check shapes
        compare_shapes(ref, actual, var)

        # Compare values: interior domain or entire array
        need_value_compare, error_logger = compare_values(ref, actual, var)

        if need_value_compare:
            # Log detailed errors with recursion if needed
            log_comparison_errors(ref, actual, var, error_logger)


def compare_shapes(ref: Dataset, actual: Dataset, var: str) -> None:
    shape_ref = ref.variables[var].shape
    shape_actual = actual.variables[var].shape
    if shape_ref != shape_actual:
        Messaging.step_error(f"Shape is not the same for variable '{var}'")
        raise Exception(f"Shape is not the same for variable '{var}'")


def compare_values(ref: Dataset, actual: Dataset, var: str) -> None:
    np_ref, np_actual = extract_comparable_arrays(ref, actual, var)

    need_value_compare = False
    # Error logger for differences
    error_logger = CompareErrorLogger(max_stored=10, max_total=50)

    if not numpy.allclose(np_ref, np_actual):
        need_value_compare = True
        error_logger.append_raw(
            var, f"Variable '{var}' not close equal via numpy.allclose()"
        )

    if (np_ref != np_actual).any():
        need_value_compare = True
        error_logger.append_raw(
            var, f"Variable '{var}' not strict equal via numpy.any()"
        )

    return need_value_compare, error_logger


def extract_comparable_arrays(ref: Dataset, actual: Dataset, var: str):
    shape_ref = ref.variables[var].shape

    # [1:-1, 1:-1] to check only interior domain
    if len(shape_ref) == 2:
        return (
            numpy.array(ref.variables[var])[1:-1, 1:-1],
            numpy.array(actual.variables[var])[1:-1, 1:-1],
        )
    elif len(shape_ref) == 3:
        return (
            numpy.array(ref.variables[var])[:, 1:-1, 1:-1],
            numpy.array(actual.variables[var])[:, 1:-1, 1:-1],
        )
    elif len(shape_ref) == 4:
        return (
            numpy.array(ref.variables[var])[:, :, 1:-1, 1:-1],
            numpy.array(actual.variables[var])[:, :, 1:-1, 1:-1],
        )
    else:
        return numpy.array(ref.variables[var]), numpy.array(actual.variables[var])


def log_comparison_errors(
    ref: Dataset, actual: Dataset, var: str, error_logger: CompareErrorLogger
) -> None:
    recurse_compare_current_dim(
        error_logger,
        var,
        ref.variables[var].shape,
        actual.variables[var].shape,
        ref.variables[var],
        actual.variables[var],
        0,
    )

    if not error_logger.has_error():
        raise Exception(
            "Ref is different from actual because need_value_compare "
            + "has been set to true. Nevertheless, error_logger did "
            + "not catch the error(s). This is a bug."
        )

    raise Exception(f"Errors detected for variable '{var}': \n{error_logger}")


##########################################################
def compare_netcdf_files(ref_file: str, actual_file: str) -> None:
    # load
    ref = Dataset(ref_file, "r")
    actual = Dataset(actual_file, "r")

    # check & forward exception if has one
    try:
        compare_netcdf_variables(ref, actual)
    except Exception as e:
        raise Exception(
            "Error while checking\n"
            + f" - refere : {ref_file}\n"
            + f" - actual : {actual_file}\n"
            + "-----------------------------------------------\n"
            + str(e)
        )
    finally:
        ref.close()
        actual.close()
