##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import shlex
import argparse
import shutil
import json

# internal
from ..helpers import Messaging, patch_lines, move_in_dir, run_shell_command
from .setup import AbstractCrocoSetup
from ..config import Config


##########################################################
class JobcompCrocoConfig:
    """
    In the old croco there is several tricks to setup OpenMP.....
    We agregate here them to keep OldCrocoSetup readable.
    """

    def __init__(self, builddir: str):
        self.builddir = builddir

    def cppdef_h_select_case(self, case_name: str):
        """
        Patch the cppdef.h to select the CASE to be used.
        """

        # progress
        Messaging.step(f"Select case : {case_name}")
        self.case = case_name

        # set patching rules
        rules = []

        # disable REGIONAL
        rules.append(
            {
                "mode": "replace",
                "what": "#define REGIONAL        /* REGIONAL Applications */\n",
                "by": "#undef REGIONAL        /* REGIONAL Applications */\n",
                "descr": "Disable default REGIONAL case",
            }
        )

        # insert the wanted one
        rules.append(
            {
                "mode": "insert-after",
                "what": "#undef REGIONAL        /* REGIONAL Applications */\n",
                "insert": f"#define {case_name}\n",
                "descr": f"Enabled wanted {case_name} case",
            }
        )

        # apply
        patch_lines(os.path.join(self.builddir, "cppdefs.h"), rules)

    def cppdef_h_set_key(self, key_name: str, status: bool):
        """
        Force enabling or disabling some cppkeys in cppdefs.h. In order to
        be applied to all cases without having to search the in case position.

        We place it just after the cases and before inclusion of cppdefs_dev.h
        which need to have everything well defined.

        Note: It is necessary to put # undef before # define as some keys can be
              already defined.
        """
        # caution here, the space in '# define' is required due to croco tricks
        # do not remove it or it behave strangely
        if status:
            status_tochange = "undef"
            status_wanted = "define"
        else:
            status_tochange = "define"
            status_wanted = "undef"

        if self.case == "REGIONAL":
            case_line_to_find = f"#if defined {self.case}"
        else:
            case_line_to_find = f"#elif defined {self.case}"

        # First, change possible existant status to match the wanted one
        # -------------------------------------------------------------------
        # build the patch rule
        rules = [
            {
                "mode": "replace",
                "what": f"# {status_tochange} {key_name}",
                "by": f"# {status_wanted} {key_name}",
                "descr": f"Set {key_name} to {status}",
            }
        ]
        # apply
        patch_lines(os.path.join(self.builddir, "cppdefs.h"), rules)
        # Then, add wanted key and status in the right case place
        # -------------------------------------------------------------------
        rules = [
            {
                "mode": "insert-after",
                "what": case_line_to_find,
                "insert": f"# {status_wanted} {key_name}",
                "descr": f"Set {key_name} to {status}",
            }
        ]
        # apply
        patch_lines(os.path.join(self.builddir, "cppdefs.h"), rules)

    def param_h_configure_openmp_split(self, splitting: str):
        """
        Configure the OpenMP splitting in param.h
        """
        # split
        dims = splitting.split("x")
        assert len(dims) == 2

        # extrat
        np_x = int(dims[0])
        np_eta = int(dims[1])
        threads = np_x * np_eta

        # build the patching rule
        rules = [
            {
                "mode": "replace",
                "what": "      parameter (NPP=4)\n",
                "by": f"      parameter (NPP={threads})\n",
                "descr": f"Set OPENMP threads = {threads}",
            }
        ]

        # apply
        patch_lines(os.path.join(self.builddir, "param.h"), rules)

        # build the patching rule
        rules = [
            {
                "mode": "replace",
                "what": "      parameter (NSUB_X=1, NSUB_E=NPP)\n",
                "by": f"      parameter (NSUB_X={np_x}, NSUB_E={np_eta})\n",
                "descr": f"Set OPENMP splitting : {splitting} : ranks={threads}, np_xi={np_x}, np_eta={np_eta}",
            }
        ]

        # apply
        patch_lines(os.path.join(self.builddir, "param.h"), rules)

    def param_h_configure_mpi_split(self, splitting: str):
        """
        Configure the MPI domain splitting and usable number of ranks.

        Parameters
        ----------
        splitting: str
            Splitting under the form of 2x8 as used for the CMake scripts.
        """
        # split
        dims = splitting.split("x")
        assert len(dims) == 2

        # extrat
        np_x = int(dims[0])
        np_eta = int(dims[1])
        ranks = np_x * np_eta

        # build rule
        rules = [
            {
                "mode": "replace",
                "what": "      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)\n",
                "by": f"      parameter (NP_XI={np_x},  NP_ETA={np_eta},  NNODES=NP_XI*NP_ETA)\n",
                "descr": f"Set MPI splitting : {splitting} : ranks={ranks}, np_xi={np_x}, np_eta={np_eta}",
            }
        ]

        # apply
        patch_lines(os.path.join(self.builddir, "param.h"), rules)


##########################################################
class JobcompCrocoSetup(AbstractCrocoSetup):
    """
    Implement configure wrapping of old croco build system to be able to also
    use bench on the master branch of croco.
    """

    def __init__(self, config: Config, builddir: str, tuning_familly: str):
        super().__init__(config, builddir)
        self.croco_config = JobcompCrocoConfig(
            builddir,
        )
        self.tuning_familly = tuning_familly
        self.fflags = ""
        self.fc = "gfortran"
        self.mpif90 = "mpif90"

    @staticmethod
    def convert_arg_for_argparse(arg_string: str) -> list:
        """
        Splits the arguments to get them as a list of args to
        be passed to argparser.

        Note: We need to reshape a bit, because the configure script
              is using `--with-xxx=value`, where argparse want with space :
              `--with-xxx value`.
        """

        # make splitting by keeping the spaces protected by "" or '' as in shell.
        basic_list = shlex.split(arg_string)

        # need to make vars (FC=gfortran) close together because argparse don't
        # like to have them mixed with options
        tmp_opts = []
        tmp_vars = []

        # pass over the args
        for entry in basic_list:
            # need to split --with-xxx=value => ['--with-xxx', 'value']
            if entry.startswith("--") and "=" in entry:
                # add a fix to avoid -key to be parse add a new argument rather than a value
                if entry.split("=", maxsplit=1)[1][0] == "-":
                    tmp_opts.append(entry.split("=", maxsplit=1)[0])
                    tmp_opts.append("," + entry.split("=", maxsplit=1)[1])
                else:
                    tmp_opts += entry.split("=", maxsplit=1)
            else:
                tmp_vars.append(entry)

        # fuse again vars + args
        return tmp_vars[1:] + tmp_opts

    def cppdef_h_set_mpi(self, status: bool):
        """
        Enable of disable MPI
        """
        self.croco_config.cppdef_h_set_key("MPI", status)

    def cppdef_h_set_openmp(self, status: bool):
        """
        Enable of disable OpenMP
        """
        self.croco_config.cppdef_h_set_key("OPENMP", status)

    def cppdef_h_enable_openacc(self, status: bool):
        """
        Enable OPENACC
        """
        self.croco_config.cppdef_h_set_key("OPENACC", status)

    def cppdef_h_enable_openacc_psyclone(self, status: bool):
        """
        Enable OPENACC_PSYCLONE
        """
        self.croco_config.cppdef_h_set_key("OPENACC_PSYCLONE", status)

    def cppdef_h_set_key(self, key: str, status: bool) -> None:
        """
        Set any cpp key
        """
        self.croco_config.cppdef_h_set_key(key, status)

    def handle_variables(self, arg_vars: list, is_mpi: bool, extra_vars: dict) -> dict:
        """
        Apply to ops required when getting some variables on the command line.

        Currently consider (exemple of values):
         - FC=gfortran
         - FFLAGS=-march=native
        """
        vars = {}
        self.fflags = ""
        for entry in arg_vars:
            # split name & value
            fields = entry.split("=", maxsplit=1)

            # extract
            var_name = fields[0]
            var_value = fields[1]

            # append extra
            if var_name in extra_vars:
                var_value += " " + " ".join(extra_vars[var_name])

            # store if we want to use latter
            vars[var_name] = var_value

            # apply vars
            if var_name == "FFLAGS":
                self.fflags = self.fflags + " " + var_value
            elif var_name == "FC":
                if is_mpi:
                    self.mpif90 = var_value
                self.fc = var_value
            else:
                raise Exception(f"Unsupported variable : {entry}")

        # ok
        return vars

    def emulate_cmd_line_parse_args(self, args: str):
        """
        Convert the original CMake + configure semantic in the old one.

        Note: This might permit to one day use this script as a configure to
              be added in the old CROCO.
        """

        # convert conventional configure agument way (with --xxx=v => -xxx v)
        args = self.convert_arg_for_argparse(args)

        # build parser
        parser = argparse.ArgumentParser(
            prog="croco_configure",
            description="A wrapper to amulate the new build API onto the old croco.",
        )

        # add options
        parser.add_argument(
            "--with-case", required=True, help="Select the case to run."
        )
        parser.add_argument(
            "--with-optim",
            required=False,
            help="Select the optimization mode to use.",
            default="seq",
        )
        parser.add_argument(
            "--with-threads",
            required=False,
            help="Select the number of threads.",
            default="4",
        )
        parser.add_argument(
            "--with-splitting",
            required=False,
            help="Select the MPI domain splitting.",
            default="1x4",
        )
        parser.add_argument(
            "--with-keys",
            required=False,
            type=str,
            help="Change cppdefs definition",
            default=None,
        )
        parser.add_argument(
            "VARS",
            nargs="*",
            type=str,
            help="Extra variable definitions, like compilers : FC=gfortran.",
        )

        # parser
        options = parser.parse_args(args)

        # ok
        return options

    def create_build_dir(self) -> None:
        """
        Create the build directory by calling `create_config.bash` so it copies
        everything.
        """
        # create directory
        with move_in_dir(self.sourcedir):
            # because this !**ù$$m script does not like getting a full path
            rel_path = os.path.relpath(self.builddir)
            run_shell_command(
                f"./create_config.bash -f -n {rel_path}",
                capture=self.config.capture,
            )

    def configure(self, args: str) -> None:
        """
        Perform the configuration.

        Parameters
        ----------
        args: str
            Arguments on the form normally passed to the
            script which will be translated in what is needed for CROCO.
        """

        Messaging.step(f"jobcomp-configure {args}")

        # get options
        options = self.emulate_cmd_line_parse_args(args)

        # create build dir
        self.create_build_dir()

        # set case
        self.croco_config.cppdef_h_select_case(options.with_case)

        # configure optimization and set flags
        self.configure_optimisation(options)

        # configure keys
        self.configure_keys(options)

        # handle extra variables
        extra_vars = {"FFLAGS": [], "LDFLAGS": []}
        self.handle_variables(options.VARS, self.use_mpi, extra_vars)

    def configure_optimisation(self, options):
        """Handle the optimisation settings."""
        self.use_openmp = False
        self.use_mpi = False
        self.use_openacc = False
        self.use_openacc_psyclone = False

        if options.with_optim == "seq":
            pass
        elif options.with_optim == "openmp":
            self.use_openmp = True
        elif options.with_optim == "mpi":
            self.use_mpi = True
        elif options.with_optim == "openacc-native":
            self.use_openacc = True
            self.use_openacc_psyclone = False
        elif options.with_optim == "openacc-psyclone":
            self.use_openacc = True
            self.use_openacc_psyclone = True
        elif options.with_optim == "poseidon":
            raise Exception(f"Variant {options.with_optim} not yet supported")
        else:
            raise Exception(
                f"Invalid optimisation mode from commande line : {options.with_optim}"
            )

        self.cppdef_h_set_mpi(self.use_mpi)
        self.cppdef_h_set_openmp(self.use_openmp)
        self.cppdef_h_enable_openacc(self.use_openacc)
        self.cppdef_h_enable_openacc_psyclone(self.use_openacc_psyclone)
        if self.use_openmp:
            self.croco_config.param_h_configure_openmp_split(options.with_splitting)
        if self.use_mpi:
            self.croco_config.param_h_configure_mpi_split(options.with_splitting)

    def configure_keys(self, options):
        """Configure additional keys if provided."""
        if options.with_keys is not None:
            # Split the string into individual key-value representations
            keys_list = options.with_keys.split(",")
            for key in keys_list:
                status = key.startswith("+")
                if key.startswith("+") or key.startswith("-"):
                    self.croco_config.cppdef_h_set_key(key[1:], status)

    def make(self, make_jobs: str):
        """
        Perform the build with jobcomp.
        """

        # move in dir & call jobcomp
        with move_in_dir(self.builddir):
            run_shell_command(
                "./jobcomp --fc %s --mpif90 %s --fflags '%s' --jobs %s"
                % (self.fc, self.mpif90, self.fflags, self.config.make_jobs),
                logfilename="jobcomp.log",
                capture=self.config.capture,
            )

    def copy_config(
        self, refdir_case: str, case_name: str, patches_and_keys: dict
    ) -> None:
        """
        Copy the required config files to pass the case in the reference dir so
        we can possibly also reproduce by hand easily if needed one day
        (outside of BENCH).

        Parameters
        ----------
        refdir_case: str
            The path in which to put the files.
        case_name: str
            Name of the case to know which file to take (possibly).
        patches_and_keys: dict
            The information from the case on what was patched.
        """

        # vars
        builddir = self.builddir

        # create subdir
        put_int = f"{refdir_case}/jobcomp-mode/"
        os.makedirs(put_int, exist_ok=True)

        # copy files from build sys
        to_copy = [
            # keep track of the exact case config being used.
            "cppdefs.h",
            "cppdefs_dev.h",
            "param.h",
        ]

        # copy them
        for file in to_copy:
            shutil.copyfile(f"{builddir}/{file}", f"{put_int}/{file}")

        # write extra
        with open(f"{put_int}/patches-and-keys.json", "w+") as fp:
            json.dump({"case": case_name, "changes": patches_and_keys}, fp)
