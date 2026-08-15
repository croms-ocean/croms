##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import shutil
import glob
import json
import numpy
import pandas
from matplotlib import pyplot
import matplotlib.dates as mdates
from .config import Config
from .messaging import Messaging


##########################################################
class Performance:
    def __init__(self, config: Config):
        self.config = config
        self.loaded_files = []

    def get_file_path(self, case_name: str, variant_name: str) -> str:
        """
        Get results file path for a specific case and variant
        """
        # name
        name = f"{case_name}-{variant_name}"
        results_pattern = self.config.results_pattern

        # calc
        fname = f"result-{name}.json"
        avail = glob.glob(f"{results_pattern}/*/*/{fname}")

        # has none
        if avail is None or len(avail) == 0:
            return None

        # we can sort by date and take the last one
        avail.sort()
        fpath = avail[-1]

        # ok
        return fpath

    def load_case_variant_data(self, case_name: str, variant_name: str) -> dict:
        """
        Load results for one variant of a case
        """
        # get path
        fpath = self.get_file_path(case_name, variant_name)

        # if not found skip
        if fpath is None:
            return None

        # progress
        Messaging.step(f"Loading {fpath}")

        # keep track
        self.loaded_files.append(fpath)

        # load it
        with open(fpath, "r") as fp:
            # load
            return json.load(fp)

    def load_variants(self, data: dict, case_name: str) -> bool:
        """
        Load results for all variants of a case
        """
        # extract vars
        variant_names = self.config.config["variants"].keys()

        # loop on all variants
        for variant_name in variant_names:
            # load json
            case_variant_data = self.load_case_variant_data(case_name, variant_name)

            # merge
            if case_variant_data is not None:
                data["variant"].append(variant_name)
                data["variant_cpu"].append(
                    int(case_variant_data["run_config"]["variant"]["requires"]["cpu"])
                )
                data["mean"].append(case_variant_data["results"][0]["mean"])
                data["median"].append(case_variant_data["results"][0]["median"])
                data["stddev"].append(case_variant_data["results"][0]["stddev"])
                data["min"].append(case_variant_data["results"][0]["min"])
                data["max"].append(case_variant_data["results"][0]["max"])
                data["q1"].append(
                    numpy.quantile(case_variant_data["results"][0]["times"], 0.25)
                )
                data["q3"].append(
                    numpy.quantile(case_variant_data["results"][0]["times"], 0.75)
                )

        # ok
        return True

    def build_data(self):
        """
        Build runtime data for all cases and variants in this run
        """
        # extract vars
        case_names = self.config.config["cases"].keys()

        # to fill
        data = {}
        # loop
        for name in case_names:
            # build
            case_data = {
                "variant": [],
                "variant_cpu": [],
                "mean": [],
                "median": [],
                "min": [],
                "max": [],
                "q1": [],
                "q3": [],
                "stddev": [],
            }

            # load
            self.load_variants(case_data, name)

            # attach
            if len(case_data["variant"]) != 0:
                data[name] = case_data

        # return
        return data

    def keep_info_of_previous_runs(self):
        """
        Save information about loaded cases from previous runs
        """
        # vars
        results = self.config.results

        # dump list of loaded paths
        with open(f"{results}/loaded-result-files.json", "w+") as fp:
            json.dump(self.loaded_files, fp=fp, indent="\t")

        # copy the loaded files localy to get a full view
        target = f"{results}/imported_previous_runs"
        os.makedirs(target, exist_ok=True)
        for file in self.loaded_files:
            # only the non local ones
            if results not in file:
                target_path = os.path.join(target, os.path.basename(file))
                shutil.copyfile(file, target_path)

    def plot(self):
        """
        Load result files plots of runtimes per case.
        """
        Messaging.section("Plotting performances")

        # extract
        cases = self.config.config["cases"]

        # get data
        data = self.build_data()

        # plot
        for case_name, case_data in data.items():
            # progress
            Messaging.step(f"Plotting {case_name}")

            # extract title
            case_config = cases[case_name]
            results = self.config.results
            title = case_config.get("title", case_name)

            # start plot
            fig, ax = pyplot.subplots()

            # build graph
            x_pos = numpy.arange(len(case_data["variant"]))
            # ax.bar(x_pos, entry['median'], yerr=[entry['min'],entry['max']], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax.bar(
                x_pos,
                case_data["median"],
                yerr=case_data["stddev"],
                align="center",
                alpha=0.5,
                ecolor="black",
                capsize=10,
            )
            ax.set_ylabel("Average runtime (seconds)")
            ax.set_xticks(x_pos)
            ax.set_xticklabels(case_data["variant"], rotation=45, ha="right")
            ax.set_title(f"CROCO - {title}")
            ax.yaxis.grid(True)

            # Save the figure and show
            pyplot.tight_layout()
            os.makedirs(f"{results}/{case_name}", exist_ok=True)
            pyplot.savefig(f"{results}/{case_name}/plot-{case_name}.png")
            pyplot.savefig(f"{results}/{case_name}/plot-{case_name}.svg")
            pyplot.close()

            # keep track of the previous run used results
            if not self.config.no_previous:
                self.keep_info_of_previous_runs()

    def track(self):
        """
        Load result files and append them to existing summary if provided.
        Returns an updated summary and plots of evolution of runtime.
        """
        Messaging.section("Keep track of performances")

        # Append to existing summary if it exists
        df_updated = self.load_perf_data()
        # Save summary file
        output_summary_file = os.path.join(self.config.results, "trackperf_summary.csv")
        Messaging.step(f"Write summary {output_summary_file}")
        df_updated.to_csv(output_summary_file, index=False)

        # Plot all cases
        Messaging.step("Plot summary")
        self.plot_perf_summary(df_updated)

    def load_perf_data(self):
        """
        Load new result files and append to existing summary if provided.
        Returns a combined DataFrame.
        """
        previous_summary_path = self.config.load_perf
        existing_df = pandas.DataFrame()

        if os.path.exists(previous_summary_path):
            Messaging.step(f"Loading previous summary {previous_summary_path}")
            existing_df = pandas.read_csv(previous_summary_path, parse_dates=["date"])
        else:
            Messaging.step("No previous summary")

        new_data = []
        # Get data
        data = self.build_data()
        # Plot
        for case_name, case_data in data.items():
            for ivariant in range(
                len(
                    case_data["variant"],
                )
            ):
                ncpu = case_data["variant_cpu"][ivariant]
                new_data.append(
                    {
                        "date": self.config.run_date,
                        "case": case_name,
                        "variant": f"{case_data['variant'][ivariant]} x {ncpu}",
                        # mean runtime * ncpu
                        "mean": case_data["mean"][ivariant] * ncpu,
                    }
                )

        new_df = pandas.DataFrame(new_data)
        combined_df = pandas.concat([existing_df, new_df], ignore_index=True)

        return combined_df

    def plot_perf_summary(self, df):
        """Plot the evolution of mean runtime"""
        results = self.config.results

        df["date"] = pandas.to_datetime(df["date"], format="%Y-%m-%d_%H-%M-%S")

        for case in df["case"].unique():
            df_case = df[df["case"] == case]

            if df_case.empty:
                print(f"No data found for case: {case}")
                return

            fig, ax = pyplot.subplots(figsize=(10, 6))

            for variant, group in df_case.groupby("variant"):
                sorted_group = group.sort_values("date")
                ax.plot(
                    sorted_group["date"],
                    sorted_group["mean"],
                    label=variant,
                    marker="o",
                )

            ax.set_title(f"Mean Runtime Evolution for Case: {case}")
            ax.set_xlabel("Date")
            ax.set_ylabel("Mean Runtime (s)")
            ax.legend()
            ax.grid(True)

            # Ticks and date format
            ax.xaxis.set_major_locator(mdates.AutoDateLocator())
            ax.xaxis.set_major_formatter(
                mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
            )

            fig.autofmt_xdate()  # Rotate labels
            pyplot.tight_layout()

            output_dir = os.path.join(results)
            os.makedirs(output_dir, exist_ok=True)
            base_path = os.path.join(output_dir, f"trackperf-{case}")
            Messaging.step(f"Save plot for case : {case}, {base_path}.png/.svg")
            pyplot.savefig(f"{base_path}.png")
            pyplot.savefig(f"{base_path}.svg")
            pyplot.close()
