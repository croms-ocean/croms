##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import shutil
import sys
import json
from typing import Union
from .config import Config
from .croco import Croco
from .messaging import Messaging
from .system import gen_system_info
from .perf import Performance
from .helpers import run_shell_command
from .htmlreport.htmlreport import generate_html


##########################################################
class Benchmarking:
    def __init__(self, config: Config):
        # set
        self.config = config

        # gen str
        result_dir = self.config.results
        case_names = ", ".join(self.config.case_names)
        variant_names = ", ".join(self.config.variant_names)

        # init
        Messaging.section("Benchmark initialization")
        Messaging.step(f"Results  = {result_dir}")
        Messaging.step(f"Cases    = {case_names}")
        Messaging.step(f"Variants = {variant_names}")

        # create work dir
        os.makedirs(config.workdir, exist_ok=True)

        # create result dir
        os.makedirs(config.results, exist_ok=True)

        # create last directory
        self.create_last_symlink()

        # make sequential first as we use as reference to validate the data produced
        self.make_ref_variant_first()

        # Create croco instances
        self.instances = self.create_croco_instances()

        # dump
        self.dump_bench_infos()

    def create_last_symlink(self):
        # extract
        results = self.config.results

        # create last directory
        dirname = os.path.dirname(results)
        basename = os.path.basename(results)
        link_name = f"{dirname}/last"
        if os.path.exists(link_name):
            os.remove(link_name)
        os.symlink(basename, link_name)

    def make_ref_variant_first(self):
        """
        Make sequential first as we use as reference to validate the data produced
        """
        # get ref
        ref_name = self.config.variant_ref_name

        # ref variant
        if ref_name in self.config.case_names:
            self.config.variant_names.remove(ref_name)
            self.config.variant_names.insert(0, ref_name)

    def supported_variantname(self, variant_name, case_config):
        res = True
        if variant_name in case_config.get("unsupported", []):
            res = False
        elif variant_name.split("-")[0] in case_config.get("unsupported", []):
            # test if openmp is supported
            res = False
        elif (
            variant_name.split("-")[0] == "mpi"
            or variant_name.split("-")[0] == "openmp"
        ):
            # test if mpi with this number of cpu is supported
            res_splitting = False
            for key, value in case_config.get("splitting", []).items():
                if str(key) == variant_name.split("-")[-1]:
                    res_splitting = True
            res = res_splitting
        return res

    def create_croco_instances(self) -> Union[Croco]:
        # extract some
        config = self.config

        # res
        res = []

        # build combination
        for case_name in config.case_names:
            for variant_name in config.variant_names:
                # get config to know how to filter
                case_config = config.config["cases"][case_name]

                # filter
                if self.supported_variantname(variant_name, case_config):
                    if config.restart:
                        if "restart" not in case_config.get("unsupported", []):
                            if variant_name == self.config.variant_ref_name:
                                restarted = False
                                res.append(
                                    Croco(config, case_name, variant_name, restarted)
                                )
                                restarted = True
                                res.append(
                                    Croco(config, case_name, variant_name, restarted)
                                )
                            else:
                                restarted = True
                                res.append(
                                    Croco(config, case_name, variant_name, restarted)
                                )
                        else:
                            Messaging.step(
                                f"Skip unsupported restart option for {case_name}"
                            )
                    else:
                        restarted = False
                        res.append(Croco(config, case_name, variant_name, restarted))
                else:
                    Messaging.step(f"Skip unsupported {case_name}/{variant_name}")

        # ok
        return res

    def run(self):
        self.build_instances()
        self.process_instances()
        self.perf_results()
        self.generate_reports()

    def build_instances(self):
        """Handles the build step."""
        if "build" in self.config.modes:
            for id, instance in enumerate(self.instances):
                instance.build(
                    extra_info=f" - [ {id + 1} / {len(self.instances)} ]",
                    force_rebuild=self.config.rebuild,
                )

    def process_instances(self):
        """Handles running, checking, plotting, and plotraw processing in one loop."""

        list_modes = ["run", "check", "plotphy", "plotraw", "animraw"]
        if any(mode in self.config.modes for mode in list_modes):
            for id, instance in enumerate(self.instances):
                instance: Croco  # Explicit annotation

                if "run" in self.config.modes:
                    instance.run(extra_info=f" - [ {id + 1} / {len(self.instances)} ]")

                if "check" in self.config.modes:
                    instance.check()

                if "plotphy" in self.config.modes:
                    instance.plotphy()

                if "plotraw" in self.config.modes:
                    instance.plotraw()

                if "animraw" in self.config.modes:
                    instance.plotraw(anim=True)

                if (
                    self.config.build_ref
                    and instance.variant_name == self.config.variant_ref_name
                ):
                    instance.make_ref()

    def perf_results(self):
        """Handles perf plotting and tracking if required."""
        if "plotperf" in self.config.modes or "trackperf" in self.config.modes:
            perf = Performance(self.config)
            if "plotperf" in self.config.modes:
                perf.plot()
            if "trackperf" in self.config.modes:
                perf.track()

    def generate_reports(self):
        """Handles reporting and HTML generation."""
        if self.config.report:
            self.config.report.save()
            self.config.report.save_to_junit()
            self.config.report.print()
        if self.config.html:
            generate_html(
                self.config.results,
                output_file=os.path.join(self.config.results, "treeview.html"),
            )
        if self.config.report.contains_false():
            Messaging.step_error("Error: False detected in report")
            sys.exit(1)

    def dump_bench_infos(self):
        # dump bench infos
        Messaging.section("Dumping system infos")

        # get needs
        results = self.config.results

        # gen
        system = gen_system_info()
        with open(f"{results}/system.json", "w+") as fp:
            json.dump(system, fp=fp, indent="\t")

        # dump config
        with open(f"{results}/config.json", "w+") as fp:
            json.dump(self.config.config, fp=fp, indent="\t")

        # display some
        hostname = system["plateform"]["hostname"]
        processor_name = system["plateform"]["processor_name"]
        Messaging.step(f"Hostname  : {hostname}")
        Messaging.step(f"Processor : {processor_name}")

        # dump the CPU infos
        if shutil.which("hwloc-ls"):
            run_shell_command(f"hwloc-ls --of console {results}/cpu.txt")
            run_shell_command(f"hwloc-ls --of svg {results}/cpu.svg")
        else:
            Messaging.step("Warning: 'hwloc-ls' not found.")
            Messaging.step("Skipping hardware topology capture.")
