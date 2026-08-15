##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import copy
import shutil
import fnmatch
import argparse
import platform
from datetime import datetime
from .messaging import Messaging
from .configfile import ConfigFile
from .report import Report
from .htmlreport.htmlreport import generate_global_html


##########################################################
class Config:
    """Configure the benchmarking"""

    def __init__(self):
        self.args = None

    def parse(self) -> None:
        """Parse the program options to init the run"""
        parser = argparse.ArgumentParser(
            prog="bench-croco.py",
            description="Automate the benchmarking of CROCO by doing rewrite.",
            formatter_class=argparse.RawTextHelpFormatter,  # Preserve formatting
        )

        # get hostname from system
        selected_host = platform.node()

        # register options
        parser.add_argument(
            "-c",
            "--cases",
            help="Specific case to run. For example BASIN.\n"
            "Could also be meta-case as defined in config file. \n"
            "(default=BASIN,BASIN-LARGE)",
            default="@default",
        )
        parser.add_argument(
            "-v",
            "--variants",
            help="Variant to run defined in config file. \n"
            "For example sequential or mpi-4.\n"
            "Could also be meta-variant. (default=@cpu)",
            default="@default",
        )
        parser.add_argument(
            "-r", "--rebuild", help="Force rebuilding.", action="store_true"
        )
        parser.add_argument(
            "-C",
            "--config",
            help="""Config file to use. (default="config.jsonc")""",
            default="config.jsonc",
        )
        parser.add_argument(
            "-s",
            "--steps",
            help="List of steps to run. \n"
            "Available steps are :\n"
            "- build : code compilation\n"
            "- run : run simulation\n"
            "- check : perform comparison with a reference (see --compare-to option)\n"
            "  By default, comparison is made over netcdf results file. \n"
            "  A more precise option is available through --cvtk.\n"
            "- plotperf : plot runtime of each variant\n"
            "- trackperf : track runtime of each case and variant in one file\n"
            "- plotphy : physical plot using python script specified in case config in plot_diag_script.\n"
            "- plotraw : plot raw map of variables at start/middle/end of simulation.\n"
            "- animraw : same plot as plotraw with an animation over the simulation.\n"
            "(default=build,run,check,plotperf)",
            dest="modes",
            default="build,run,check,plotperf",
        )

        parser.add_argument(
            "--compare-to",
            help="Set which variant is used as reference (default='sequential').",
            type=str,
            default="sequential",
        )
        parser.add_argument(
            "--cvtk",
            help="Enable usage of CVTK_DEBUG for comparison with a reference",
            action="store_true",
        )
        parser.add_argument(
            "-R",
            "--runs",
            help="Number or runs to perform. \n"
            "Several runs could be used to make performance tests. (default=1)",
            default="1",
        )
        parser.add_argument(
            "-l",
            "--load-perf",
            help="Load previous runtime tracking",
            type=str,
            default="",
        )
        parser.add_argument(
            "-j",
            "--jobs",
            help="Make -j option value to build \n(default=4)",
            default="4",
        )
        parser.add_argument(
            "-a",
            "--auto-skip",
            help="Automatically skip the cases with not enougth ressources.",
            action="store_true",
        )
        parser.add_argument(
            "-N",
            "--no-previous",
            help="Do not load the previous missing results to plot full graph.",
            action="store_true",
        )
        parser.add_argument(
            "--build-ref",
            help="Build the reference directory to be stored once somewhere for validation.",
            type=str,
            default=False,
        )
        parser.add_argument(
            "-u",
            "--use-ref",
            help="Use the reference directory as source of compare instead of seq run.",
            type=str,
            default=False,
        )
        parser.add_argument(
            "--host",
            help="Force host config to use.",
            type=str,
            default=selected_host,
        )
        parser.add_argument(
            "--data-root-path",
            help="Input data path. Needed if cases to run needs external  \n"
            "data. In this case, specified path should contain cases specific directories \n"
            "(croco_files_path in case config file) (default='')",
            type=str,
            default="",
        )
        parser.add_argument(
            "--results",
            help="""Name of the results directory to use. (default="results")""",
            default="results",
        )
        parser.add_argument(
            "-t",
            "--title",
            help="A possible extra name to prepent to the case result directory name\n"
            "(default=None)",
            default=None,
        )
        parser.add_argument(
            "-w",
            "--workdir",
            help="""Working directory name. (default = "rundir")""",
            default="rundir",
        )
        parser.add_argument(
            "--continue",
            help="Continue instead of stopping on first error.",
            dest="continue_on_error",
            action="store_true",
        )
        parser.add_argument(
            "--report", help="Build a status report.", action="store_true"
        )
        parser.add_argument(
            "--html", help="Build a html treeview report.", action="store_true"
        )
        parser.add_argument(
            "--globalhtml",
            help="Build a global html report for all results folders. \n"
            "In this case bench steps are not run !",
            action="store_true",
        )
        parser.add_argument(
            "--debug", help="Run with debug options", action="store_true"
        )
        parser.add_argument(
            "--restart", help="Check restartability", action="store_true"
        )
        parser.add_argument(
            "--no_date_in_result_dir",
            help="Do not add date in result forlder name",
            dest="no_date_in_result_dir",
            action="store_true",
        )
        parser.add_argument(
            "-V",
            "--verbose",
            help="Enable the full logging instead of printing just summaries.",
            action="store_true",
        )

        # parse
        self.args = parser.parse_args()

        # reshape
        self.modes = self.args.modes.split(",")
        self.variant_names = self.args.variants.split(",")
        self.case_names = self.args.cases.split(",")
        self.verbose = self.args.verbose
        self.capture = not self.verbose
        self.workdir = os.path.abspath(self.args.workdir)
        self.rebuild = self.args.rebuild
        self.runs = int(self.args.runs)
        self.load_perf = self.args.load_perf
        self.make_jobs = int(self.args.jobs)
        self.title = self.args.title
        self.auto_skip = self.args.auto_skip
        self.no_previous = self.args.no_previous
        self.build_ref = self.args.build_ref
        self.use_ref = self.args.use_ref
        self.use_host_config = self.args.host
        self.variant_ref_name = self.args.compare_to
        self.cvtk = self.args.cvtk
        self.continue_on_error = self.args.continue_on_error
        self.report = self.args.report
        self.html = self.args.html
        self.globalhtml = self.args.globalhtml
        self.data_root_path = self.args.data_root_path
        self.debug = self.args.debug
        self.restart = self.args.restart

        # swithc
        if self.report:
            self.continue_on_error = True

        # compute clean result subdir name
        use_host_config = self.use_host_config
        if self.title is not None:
            folder_name_base = f"{self.title}-{use_host_config}"
        else:
            folder_name_base = f"{use_host_config}"

        self.run_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        if self.args.no_date_in_result_dir:
            folder_name = f"{folder_name_base}"
        else:
            folder_name = f"{folder_name_base}-{self.run_date}"
        self.results = os.path.join(
            self.args.results,
            folder_name,
        )
        # rm the old one if it exists
        shutil.rmtree(self.results, ignore_errors=True)

        # pattern to search results files (to also take the previous runs if not re-run all)
        if self.no_previous:
            self.results_pattern = self.results
        else:
            self.results_pattern = os.path.join(
                self.args.results,
                f"{folder_name_base}-*",
            )

        # extract some many time used paths
        self.croco_source_dir = os.path.abspath(f"{__file__}/../../../")

        # global report
        self.report = Report(self)

    def pre_checks(self):
        # basic check
        if self.use_ref and not os.path.exists(self.use_ref):
            raise Exception(
                f"You gave --use-ref={self.use_ref}, but directory does not exist !"
            )
        if self.build_ref and self.variant_ref_name not in self.variant_names:
            raise Exception(
                f"You gave --build-ref={self.use_ref}, but not using the required '{self.variant_ref_name}' variant !"
            )

    def filter_variant_ressources(self):
        # extract needed vars
        variants = self.config["variants"]
        host = self.host
        host_ressources = host["ressources"]

        # filter
        for name in self.variant_names.copy():
            variant = variants[name]
            for ressource_name, ressource_min in variant["requires"].items():
                has = host_ressources[ressource_name]
                if has < int(ressource_min):
                    Messaging.step(
                        f"Filter out {name} due to lack of {ressource_name} ({has} < {ressource_min})..."
                    )
                    if self.auto_skip:
                        self.variant_names.remove(name)
                    else:
                        raise Exception(
                            f"Cannot run {name} due to lack of {ressource_name} ({has} < {ressource_min})... You can use -a/--auto-skip to skip silently."
                        )

    def load_config(self) -> None:
        """Load the config file and extract some infos"""
        # load
        self.config = ConfigFile(self.args.config).load()

        # message
        Messaging.section("Applying selection")

        # apply meta : @.....
        self.apply_meta()

        # select host
        self.select_host()

        # filter ressources
        self.filter_variant_ressources()

        # check
        self.pre_checks()

    def apply_meta(self):
        # vars
        avail_variants = self.config["variants"].keys()
        avail_cases = self.config["cases"].keys()

        # apply default
        if self.variant_names == ["@default"]:
            self.variant_names = self.config["default_selection"]["variants"]
        if self.case_names == ["@default"]:
            self.case_names = self.config["default_selection"]["cases"]

        # apply user selection
        if self.variant_names == ["@all"]:
            self.variant_names = avail_variants
        if self.case_names == ["@all"]:
            self.case_names = avail_cases

        # extract variant groups
        final_variants = []
        for variant in self.variant_names:
            if variant.startswith("@"):
                final_variants += self.config["meta_variants"][variant]
            else:
                final_variants.append(variant)
        self.variant_names = final_variants

        # extract cases groups
        final_cases = []
        for case in self.case_names:
            if case.startswith("@"):
                final_cases += self.config["meta_cases"][case]
            else:
                final_cases.append(case)
        self.case_names = final_cases

        # apply whildcard
        self.case_names = self.apply_wildcard_enabling(self.case_names, avail_cases)
        self.variant_names = self.apply_wildcard_enabling(
            self.variant_names, avail_variants
        )

        # filter
        self.case_names = self.apply_minus_on_list(self.case_names)
        self.variant_names = self.apply_minus_on_list(self.variant_names)

    @staticmethod
    def apply_wildcard_enabling(data: list, avails: list) -> list:
        result = []
        pattern: str
        for pattern in data:
            if pattern.startswith("-"):
                result.append(pattern)
            else:
                sublist = fnmatch.filter(avails, pattern)
                sublist.sort()
                result += sublist
        return result

    @staticmethod
    def apply_minus_on_list(data: list) -> list:
        result = copy.copy(data)
        value: str
        for value in data:
            if value.startswith("-"):
                to_remove = fnmatch.filter(result, value[1:])
                for item in to_remove:
                    result.remove(item)
                result.remove(value)
        return result

    def select_host(self):
        # get it
        hostname = platform.node()

        # select tunning
        selected_host = self.use_host_config
        hosts = self.config["hosts"]

        # fallback
        if selected_host not in hosts:
            selected_host = "@default"

        # display
        Messaging.step(f"Hostname : {hostname}")
        Messaging.step(f"Select host : {selected_host}")

        # apply
        self.host = hosts[selected_host]

    def html_global_report(self):
        # global html
        if self.globalhtml:
            generate_global_html(self.args.results)
