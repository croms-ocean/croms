##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import glob
import json
import shutil
import platform
import subprocess
from datetime import timedelta
import math


# internal
from .config import Config
from .helpers import (
    move_in_dir,
    apply_vars_in_str,
    Messaging,
    patch_lines,
    copy_tree_with_absolute_symlinks,
    extract_elements_from_file,
    copy_and_replace,
    delete_lines_from_file,
    parse_datetime,
)
from .hyperfine import run_hyperfine
from .check import compare_netcdf_files
from .crocosetup.jobcomp import JobcompCrocoSetup
from .rawdrawer import RawDrawer


##########################################################
class Croco:
    def __init__(
        self, config: Config, case_name: str, variant_name: str, restarted: bool
    ):
        # keep track
        self.config = config
        self.case_name = case_name
        self.variant_name = variant_name
        self.restarted = restarted
        if self.restarted:
            self.full_name = f"{case_name}-{variant_name}-rst"
        else:
            self.full_name = f"{case_name}-{variant_name}"

        # build directory name
        self.dirname = self.calc_rundir(variant_name, case_name, restarted)
        self.dirname = os.path.abspath(self.dirname)
        self.dirname_result = self.calc_resdir(variant_name, case_name, restarted)
        self.dirname_result = os.path.abspath(self.dirname_result)

        # extract infos
        self.case = config.config["cases"][case_name]
        self.variant = config.config["variants"][variant_name]

        # input file
        self.croco_inputfile = "TEST_CASES/croco.in.%s" % self.case["case"].capitalize()
        if "input_file" in self.case:
            if len(self.case["input_file"]) > 0:
                self.croco_inputfile = self.case["input_file"]

        self.croco_build = JobcompCrocoSetup(
            self.config, self.dirname, self.variant["tuning_familly"]
        )

    def calc_rundir(self, variant_name, case_name, restarted):
        if restarted:
            return f"{self.config.workdir}/{case_name}/{variant_name}-rst"
        else:
            return f"{self.config.workdir}/{case_name}/{variant_name}"

    def calc_resdir(self, variant_name, case_name, restarted):
        if restarted:
            return f"{self.config.results}/{case_name}/{variant_name}-rst"
        else:
            return f"{self.config.results}/{case_name}/{variant_name}"

    def reset(self):
        # display
        Messaging.step("Reset")

        # rm the old one
        shutil.rmtree(self.dirname, ignore_errors=True)

    def manage_dir_and_paths(self):
        dirname = self.dirname
        dirname_result = self.dirname_result
        os.makedirs(dirname, exist_ok=True)
        os.makedirs(dirname_result, exist_ok=True)
        # link data dir
        if "croco_files_path" in self.case:
            self.input_dir = os.path.join(
                self.config.data_root_path, self.case["croco_files_path"]
            )
            if not (os.path.isdir(self.config.data_root_path)):
                Messaging.step_error(
                    "Folder data_root_path not found :  %s" % self.config.data_root_path
                )
                raise Exception("Folder not found : %s" % self.config.data_root_path)

            if not (os.path.isdir(self.input_dir)):
                Messaging.step_error(
                    "Folder for input data not found :  %s" % self.input_dir
                )
                raise Exception("Folder not found : %s" % self.input_dir)
            copy_tree_with_absolute_symlinks(self.input_dir, dirname)

    def configure(self):
        # Configure
        Messaging.step("Configure")
        # extract some needed vars
        croco_source_dir = self.config.croco_source_dir

        tuning_familly = self.variant["tuning_familly"]
        if self.config.debug:
            tuning_flags = self.config.host["tuning"][tuning_familly]["debug"]
        else:
            tuning_flags = self.config.host["tuning"][tuning_familly]["standard"]
        tuning_flags_extra = self.config.host["tuning"][tuning_familly]["extra"]

        # create dir
        self.manage_dir_and_paths()

        # extract options
        configure_variant_options = self.variant["configure"]
        configure_case_option = f"--with-case={self.case['case']}"
        configure_compiler_option = ""
        if tuning_flags != "":
            configure_compiler_option = f'FFLAGS="{tuning_flags} {tuning_flags_extra}"'

        # Add cppkeys
        configure_cppkeys_options = ""
        if "cppkeys" in self.case:
            reshape = []
            for key, value in self.case["cppkeys"].items():
                if value:
                    reshape.append(f"+{key}")
                else:
                    reshape.append(f"-{key}")
            if len(reshape) > 0:
                reshape_str = ",".join(reshape)
                configure_cppkeys_options = f"--with-keys={reshape_str}"

        # build command
        command = "%s/configure %s %s %s %s" % (
            croco_source_dir,
            configure_compiler_option,
            configure_case_option,
            configure_cppkeys_options,
            configure_variant_options,
        )
        command = apply_vars_in_str(
            command, {"case": self.case, "tuning": self.config.host["tuning"]}
        )

        # jump in & configure
        self.croco_build.configure(command)

    def compile(self):
        # display
        Messaging.step("Compile")

        # extract some needed vars
        croco_build = self.croco_build
        make_jobs = f"-j{self.config.make_jobs}"

        # jump in & build
        croco_build.make(make_jobs)

    def enable_cvtk_checking(self):
        # vars
        croco_build = self.croco_build
        variant_name = self.variant_name
        cvtk_ref_variant_name = self.config.variant_ref_name

        # enable what we need
        croco_build.cppdef_h_set_key("CVTK_DEBUG", True)

        # except for variant we enable CVTK_DEBUG_READ
        if variant_name == cvtk_ref_variant_name and not self.restarted:
            croco_build.cppdef_h_set_key("CVTK_DEBUG_WRITE", True)
        else:
            croco_build.cppdef_h_set_key("CVTK_DEBUG_READ", True)

        if self.config.restart:
            croco_build.cppdef_h_set_key("CVTK_DEBUG_PERFRST", True)

    def build(self, extra_info: str = "", force_rebuild: bool = False):
        build_error = False
        build_error_details = None
        if self.config.continue_on_error:
            try:
                self.build_internal(extra_info, force_rebuild)
            except Exception as e:
                Messaging.step_error("Fail to build !")
                build_error = True
                build_error_details = e

        else:
            self.build_internal(extra_info, force_rebuild)

        if build_error:
            self.config.report.report_status(
                self.case_name,
                self.variant_name,
                self.restarted,
                "build",
                False,
                str(build_error_details),
            )
        else:
            self.config.report.report_status(
                self.case_name, self.variant_name, self.restarted, "build", True
            )
        # copy log file in result dir
        if os.path.exists(os.path.join(self.dirname, "jobcomp.log")):
            os.makedirs(self.dirname_result, exist_ok=True)
            shutil.copyfile(
                os.path.join(self.dirname, "jobcomp.log"),
                os.path.join(self.dirname_result, "jobcomp.log"),
            )

    def build_internal(self, extra_info: str = "", force_rebuild: bool = False):
        # display
        Messaging.section(f"Building CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")

        # perform steps
        if force_rebuild:
            self.reset()
        elif not force_rebuild and os.path.exists(f"{self.dirname}/croco"):
            Messaging.step("Already built...")
            return

        # effectively built
        self.configure()
        self.setup_case()
        self.setup_variant()

        # some specific handling
        Messaging.step("Special command line configs...")
        if self.config.cvtk:
            self.enable_cvtk_checking()

        if self.restarted:
            croco_build = self.croco_build
            croco_build.cppdef_h_set_key("EXACT_RESTART", True)

        # compile
        self.compile()

    def symlink_all_cvtk_ref_files(self):
        # vars
        case_name = self.case_name
        cvtk_ref_variant_name = self.config.variant_ref_name

        # refdir
        ref_rundir = self.calc_rundir(cvtk_ref_variant_name, case_name, False)
        refdir = os.path.join(ref_rundir, "check_file_*")

        # debug_infos
        extra_files = [f"{ref_rundir}/debug_infos"]

        # link all
        for file in glob.glob(refdir) + extra_files:
            if not os.path.exists(os.path.basename(file)):
                os.symlink(file, os.path.basename(file))

    def run(self, extra_info: str = ""):
        # display
        Messaging.section(f"Running CROCO - {self.full_name}{extra_info}")
        Messaging.step(f"Directory: {self.dirname}")

        # extract some vars
        command_prefix = self.variant["command_prefix"]
        environ = self.variant.get("environ", {}).copy()
        dirname = self.dirname
        runs = self.config.runs
        host_tuning = self.config.host["tuning"]
        cvtk = self.config.cvtk
        is_cvtk_ref = (
            self.variant_name == self.config.variant_ref_name and not self.restarted
        )
        restart = self.restarted

        # load tuning env & override if needed
        for key, value in host_tuning["environ"].items():
            environ[key] = value

        # apply vars
        command_prefix = apply_vars_in_str(
            command_prefix, {"case": self.case, "tuning": host_tuning}
        )

        # build end
        env_line = ""
        for var, value in environ.items():
            env_line += f'{var}="{value}" '

        # build command and run
        command = "%s ../../../scripts/correct_end.sh %s ./croco %s" % (
            env_line,
            command_prefix,
            self.croco_inputfile,
        )
        if restart:
            # execute twice one without restart and one with
            command = "%s && %s_rst" % (
                command.replace("../../../scripts/correct_end.sh", "").strip(),
                command,
            )

        with move_in_dir(dirname):
            # link ref files
            if cvtk and not is_cvtk_ref:
                self.symlink_all_cvtk_ref_files()

            # run
            if self.config.continue_on_error:
                try:
                    result = run_hyperfine(
                        command,
                        logfilename="croco.log",
                        runs=runs,
                        verbose=self.config.verbose,
                    )
                except Exception as e:
                    Messaging.step_error("Error while running...")
                    result = {
                        "results": [
                            {
                                "command": command,
                                "mean": 0,
                                "median": 0,
                                "stddev": 0,
                                "min": 0,
                                "max": 0,
                                "times": [0],
                            }
                        ],
                        "error": str(e),
                    }
            else:
                result = run_hyperfine(
                    command,
                    logfilename="croco.log",
                    runs=runs,
                    verbose=self.config.verbose,
                )

        # status
        if "error" in result:
            self.config.report.report_status(
                self.case_name,
                self.variant_name,
                self.restarted,
                "run",
                False,
                result["error"],
            )
        else:
            self.config.report.report_status(
                self.case_name, self.variant_name, self.restarted, "run", True
            )

        # calc output
        os.makedirs(self.dirname_result, exist_ok=True)
        output_filename = f"{self.dirname_result}/result-{self.full_name}.json"

        # add some infos
        result["name"] = self.full_name
        result["original_path"] = output_filename
        result["hostname"] = platform.node()
        result["selected_host_config"] = self.config.use_host_config
        result["run_config"] = {
            "host": self.config.host,
            "case": self.case,
            "variant": self.variant,
        }

        # dump results
        with open(output_filename, "w+") as fp:
            json.dump(result, fp=fp, indent="\t")

        # copy log file in result dir
        for iteration in range(runs):
            logfilename = "croco_%i.log" % (iteration + 1)
            if os.path.exists(os.path.join(self.dirname, logfilename)):
                shutil.copyfile(
                    os.path.join(self.dirname, logfilename),
                    os.path.join(self.dirname_result, logfilename),
                )

    def apply_patches(self, patches):
        # loop
        with move_in_dir(self.dirname):
            for file, changes in patches.items():
                # convert to list if needed
                if isinstance(changes, dict):
                    changes = [changes]

                # loop on all changes
                for change in changes:
                    patch_lines(file, [change])

    def apply_debug_patches(self):
        filename = self.croco_inputfile
        self.change_card_time_stepping_ntimes(filename, 6)
        self.change_card_history_nwrt(filename)
        # and for USE_CALENDAR
        self.change_card_end_date(filename, 6)
        self.change_card_output_time_steps_dthis(filename, 6)

    def apply_restart_patches(self):
        filename = self.croco_inputfile

        # for all case (write/read), put ldefhis to F
        self.change_card_history_ldefhis(filename)

        if self.restarted:
            # prepare 2 files for the restarted run
            full_filename = os.path.join(self.dirname, filename)
            filename_rst = filename + "_rst"
            full_filename_rst = os.path.join(self.dirname, filename_rst)
            shutil.copy(full_filename, full_filename_rst)
            file_nc_rst = "croco_restart.nc"

            # first run with filename
            self.change_card_time_stepping_ntimes(filename, 3)
            self.change_card_restart(filename, 3, file_nc_rst)

            # second run with filename_rst
            self.change_card_time_stepping_ntimes(filename_rst, 3)
            self.change_card_initial(filename_rst, file_nc_rst)

            # and for USE_CALENDAR
            self.change_card_end_date(filename, 3)
            self.change_card_output_time_steps_dtrst(filename, 3)
            # no need to change end_date or dtrsr for filename_rst

    def change_card_restart(self, filename, nrst, file_nc_rst):
        full_filename = os.path.join(self.dirname, filename)
        patches = {
            filename: {
                "file": filename,
                "mode": "insert-after",
                "what": "restart:",
                "insert": ["%i   0" % nrst, file_nc_rst],
                "descr": "change restart NRST=3",
            }
        }
        self.apply_patches(patches)
        delete_lines_from_file(full_filename, "restart", line_offset=3, num_lines=2)

    def change_card_initial(self, filename, file_nc_rst):
        full_filename = os.path.join(self.dirname, filename)
        patches = {
            filename: {
                "mode": "insert-after",
                "what": "initial:",
                "insert": ["     2", file_nc_rst],
                "descr": "change restart for reading step to NRREC=2 and file=%s"
                % file_nc_rst,
            }
        }
        self.apply_patches(patches)
        delete_lines_from_file(full_filename, "initial", line_offset=3, num_lines=2)

    def change_card_output_time_steps_dthis(self, filename, ntimes, min_dt=1.0):
        full_filename = os.path.join(self.dirname, filename)
        # Check end_date is a card in this case
        if len(extract_elements_from_file(full_filename, "output_time_steps")) > 0:
            TIME_LINE = extract_elements_from_file(full_filename, "time_stepping")
            OUTPUT_TIME_STEPS = extract_elements_from_file(
                full_filename, "output_time_steps"
            )
            dt = float(TIME_LINE[1])
            duration = math.ceil(max(dt * ntimes, min_dt))
            dt_his_hours = max(dt / 3600.0, duration / (ntimes * 3600))
            # in case of very small dt put put a minimum
            NEW_OUTPUT_TIME_STEPS = copy_and_replace(OUTPUT_TIME_STEPS, 0, dt_his_hours)
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " output_time_steps:",
                    "insert": " ".join(map(str, NEW_OUTPUT_TIME_STEPS)),
                    "descr": f"change output_time_steps to DT_HIS(H)={dt_his_hours}",
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "output_time_steps", line_offset=2, num_lines=1
            )

    def change_card_output_time_steps_dtrst(self, filename, ntimes, min_dt=1.0):
        full_filename = os.path.join(self.dirname, filename)
        # Check end_date is a card in this case
        if len(extract_elements_from_file(full_filename, "output_time_steps")) > 0:
            TIME_LINE = extract_elements_from_file(full_filename, "time_stepping")
            OUTPUT_TIME_STEPS = extract_elements_from_file(
                full_filename, "output_time_steps"
            )
            dt = float(TIME_LINE[1])
            duration = math.ceil(max(dt * ntimes, min_dt))
            dt_rst_hours = duration / 3600.0
            # in case of very small dt put put a minimum
            NEW_OUTPUT_TIME_STEPS = copy_and_replace(OUTPUT_TIME_STEPS, 2, dt_rst_hours)
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " output_time_steps:",
                    "insert": " ".join(map(str, NEW_OUTPUT_TIME_STEPS)),
                    "descr": f"change output_time_steps to DT_RST(H)={dt_rst_hours}",
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "output_time_steps", line_offset=2, num_lines=1
            )

    def change_card_end_date(self, filename, ntimes, min_dt=1.0):
        full_filename = os.path.join(self.dirname, filename)
        # Check end_date is a card in this case
        if len(extract_elements_from_file(full_filename, "end_date")) > 0:
            TIME_LINE = extract_elements_from_file(full_filename, "time_stepping")
            START_DATE = extract_elements_from_file(full_filename, "start_date")
            dt = float(TIME_LINE[1])
            duration = math.ceil(max(dt * ntimes, min_dt))
            # in case of very small dt put a minimum
            datetime_start = parse_datetime(START_DATE[0] + " " + START_DATE[1])
            datetime_end = datetime_start + timedelta(seconds=duration)
            end_date = datetime_end.strftime("%Y-%m-%d %H:%M:%S")
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " end_date:",
                    "insert": end_date,
                    "descr": f"change end_date to {end_date}",
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "end_date", line_offset=2, num_lines=1
            )

    def change_card_time_stepping_ntimes(self, filename, ntimes):
        full_filename = os.path.join(self.dirname, filename)
        # Check time_stepping is a card in this case
        if len(extract_elements_from_file(full_filename, "time_stepping")) > 0:
            TIME_LINE = extract_elements_from_file(full_filename, "time_stepping")
            NEW_TIME_LINE = copy_and_replace(TIME_LINE, 0, ntimes)
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " time_stepping:",
                    "insert": " ".join(map(str, NEW_TIME_LINE)),
                    "descr": "change duration to NTIMES=%i" % ntimes,
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "time_stepping", line_offset=2, num_lines=1
            )

    def change_card_history_nwrt(self, filename):
        full_filename = os.path.join(self.dirname, filename)
        # Check time_stepping is a card in this case
        if len(extract_elements_from_file(full_filename, "history")) > 0:
            HISTORY_LINE = extract_elements_from_file(full_filename, "history")
            NEW_HISTORY_LINE = copy_and_replace(HISTORY_LINE, 1, 1)
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " history:",
                    "insert": " ".join(map(str, NEW_HISTORY_LINE)),
                    "descr": "change output to NWRT=1",
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "history:", line_offset=2, num_lines=1
            )

    def change_card_history_ldefhis(self, filename):
        full_filename = os.path.join(self.dirname, filename)
        # Check time_stepping is a card in this case
        if len(extract_elements_from_file(full_filename, "history")) > 0:
            HISTORY_LINE = extract_elements_from_file(full_filename, "history")
            NEW_HISTORY_LINE = copy_and_replace(HISTORY_LINE, 0, "F")
            patches = {
                filename: {
                    "mode": "insert-after",
                    "what": " history:",
                    "insert": " ".join(map(str, NEW_HISTORY_LINE)),
                    "descr": "change output to LDEFHIS=F",
                }
            }
            self.apply_patches(patches)
            delete_lines_from_file(
                full_filename, "history:", line_offset=2, num_lines=1
            )

    def setup_case(self):
        # apply the case paches
        case_name = self.case_name

        Messaging.step(f"Apply case config : {case_name}")
        patches = self.case.get("patches", {})
        self.apply_patches(patches)

        if self.config.debug or self.config.restart:
            Messaging.step("Reduce number of time steps")
            self.apply_debug_patches()

        if self.config.restart:
            Messaging.step("Prepare files for restarted run")
            self.apply_restart_patches()

    def setup_variant(self):
        # apply the case paches
        variant_name = self.variant_name
        patches = self.variant.get("patches", {})

        # display
        Messaging.step(f"Apply variant config : {variant_name}")
        self.apply_patches(patches)

    def check_one_file_from_seq_ref(self, filename: str) -> None:
        # extract vars
        case_name = self.case_name
        dirname = self.dirname
        variant_ref_name = self.config.variant_ref_name
        Messaging.step(f"Checking {case_name} / {filename}")
        seq_dir = self.calc_rundir(variant_ref_name, case_name, False)
        seq_file = os.path.join(seq_dir, filename)
        if not os.path.exists(seq_file):
            Messaging.step_error(
                "Reference file does not exist, check you ran reference : %s "
                % variant_ref_name
            )
            raise Exception(
                "Missing '%s', are you sure you ran case '%s' first to get a reference for checks ?"
                % (seq_file, variant_ref_name)
            )

        # compare
        actual_file = f"{dirname}/{filename}"
        compare_netcdf_files(seq_file, actual_file)

    def check_one_file_from_ref_dir(self, filename: str) -> None:
        # vars
        refdir = self.config.use_ref
        case_name = self.case_name
        dirname = self.dirname

        # progress
        Messaging.step(f"Checking {case_name} / {filename} from ref={refdir}")

        # error
        ref_file = f"{refdir}/{case_name}/{filename}"
        if not os.path.exists(ref_file):
            Messaging.step_error(
                "Reference file does not exist, check your reference directory : %s "
                % refdir
            )
            raise Exception(
                "Missing '%s', are you sure you ran case %s in reference directory %s ?"
                % (ref_file, case_name, refdir)
            )

        # compare
        actual_file = f"{dirname}/{filename}"
        compare_netcdf_files(ref_file, actual_file)

    def check_one_file(self, filename: str):
        # compare
        if self.config.continue_on_error:
            try:
                self.check_one_file_internal(filename)
            except Exception as e:
                Messaging.step_error("Failed to compare !")
                self.config.report.report_status(
                    self.case_name,
                    self.variant_name,
                    self.restarted,
                    "check",
                    False,
                    str(e),
                )
                return
        else:
            self.check_one_file_internal(filename)

        # report
        self.config.report.report_status(
            self.case_name, self.variant_name, self.restarted, "check", True
        )

    def check_one_file_internal(self, filename: str):
        # vars
        refdir = self.config.use_ref

        # swtich
        if refdir:
            self.check_one_file_from_ref_dir(filename)
        else:
            self.check_one_file_from_seq_ref(filename)

    def check(self):
        if (
            not self.config.use_ref
            and self.variant_name == self.config.variant_ref_name
            and not self.restarted
        ):
            Messaging.step(f"Checking skiped for '{self.config.variant_ref_name}'")
            return
        else:
            for filename in self.case["check_outputs"]:
                self.check_one_file(filename)

    def plotphy(self):
        dirname = self.dirname
        full_name = self.full_name
        result_dir = self.dirname_result

        if "plotphy" in self.case:
            if "plotphy_script" in self.case["plotphy"]:
                self.plot_diag_script = os.path.join(
                    dirname, self.case["plotphy"]["plotphy_script"]
                )
                Messaging.step(f"Plotting {full_name}")

                filename = self.case["plotphy"]["plotphy_script_input"]
                actual_file = f"{dirname}/{filename}"

                command = [
                    self.plot_diag_script,
                    "--no-show",
                    "--makepng",
                    "--file",
                    actual_file,
                    "--output-dir",
                    result_dir,
                ]

                try:
                    subprocess.run(command, check=True)  # Exécute la commande
                    Messaging.step(
                        "Successfully executed %s with arguments --no-show --makepng"
                        % self.plot_diag_script
                    )
                    # report
                    self.config.report.report_status(
                        self.case_name,
                        self.variant_name,
                        self.restarted,
                        "plotphy",
                        True,
                    )
                except subprocess.CalledProcessError as e:
                    Messaging.step_error(
                        f"Error during execution of {self.plot_diag_script}: {e}"
                    )
                    self.config.report.report_status(
                        self.case_name,
                        self.variant_name,
                        self.restarted,
                        "plotphy",
                        False,
                        str(e),
                    )
            else:
                Messaging.step("No Plotting script provided")
        else:
            Messaging.step("Plotphy not available for this case")

    def plotraw(self, anim: bool = False):
        # vars
        dirname = self.dirname
        full_name = self.full_name
        result_dir = self.dirname_result

        # loop each
        for filename in self.case["check_outputs"]:
            actual_file = f"{dirname}/{filename}"
            dump_in_file = f"{result_dir}/plotraw/"
            if anim:
                Messaging.step(f"Drawing plotraw-anim {full_name} / {filename}")
            else:
                Messaging.step(f"Drawing plotraw {full_name} / {filename}")
            raw_drawer = RawDrawer(actual_file, self.config, self.case_name)
            if anim:
                raw_drawer.plot_all_variables_animation(dump_in_file)
            else:
                raw_drawer.plot_all_variables(dump_in_file)

    def make_ref(self):
        # extract vars
        variant_ref_name = self.config.variant_ref_name
        refdir = self.config.build_ref
        case_name = self.case_name
        dirname = self.dirname
        patches_and_keys = {
            "case": {
                "patches": self.case.get("patches", {}),
                "keys": self.case.get("keys", {}),
            },
            "variant": {
                "patches": self.variant.get("patches", {}),
                "keys": self.variant.get("keys", {}),
            },
        }

        # check
        assert refdir

        # progress
        Messaging.step(f"Storing ref {case_name} into ref={refdir}")

        # we create only for seq
        if self.variant_name != variant_ref_name:
            return

        # create dir if not exist
        refdir_case = f"{refdir}/{case_name}"
        os.makedirs(refdir_case, exist_ok=True)

        # copy the files there
        for filename in self.case["check_outputs"]:
            orig_path = f"{dirname}/{filename}"
            dest_path = f"{refdir_case}/{filename}"
            Messaging.step(f"Storing file {case_name} / {filename}")
            shutil.copyfile(orig_path, dest_path)

        # also copy the config
        self.croco_build.copy_config(refdir_case, case_name, patches_and_keys)

        # add the case config files
        case_file = self.croco_inputfile
        os.makedirs(f"{refdir}/{case_name}/TEST_CASES", exist_ok=True)
        shutil.copyfile(f"{dirname}/{case_file}", f"{refdir}/{case_name}/{case_file}")

        # copy the case file under croco.in
        if not os.path.exists(f"{refdir}/{case_name}/croco.in"):
            os.symlink(case_file, f"{refdir}/{case_name}/croco.in")

        # dump case info
        with open(f"{refdir}/{case_name}/case.json", "w+") as fp:
            json.dump(self.case, fp, indent="\t")

        # dump git info
        with open(f"{refdir}/{case_name}/git.log", "w+") as fp:
            fp.write(
                subprocess.getoutput(
                    "git log HEAD~1..HEAD 2>/dev/null || echo 'No git !'"
                )
                + "\n"
            )
