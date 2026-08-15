##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""
Hyperfine (https://github.com/sharkdp/hyperfine) is a runner to measure perf
of a command.
I provide here a python short replacement as it is not installed by default
systems.
"""

##########################################################
import numpy
import json
from tqdm import tqdm
from tempfile import NamedTemporaryFile
from .messaging import Messaging
from .helpers import run_shell_command_time, run_shell_command


##########################################################
def emulate_hyperfine(
    command: str,
    logfilename=None,
    progress_bar: tqdm = None,
    runs: int = 2,
    retries: int = 8,
    verbose: bool = False,
) -> dict:
    # prepare some vars
    measures = []

    # display
    Messaging.command(command)

    # update progress bar
    if progress_bar is not None:
        progress_bar.refresh()

    # repeat
    for iteration in range(runs):
        # run
        Messaging.step(f"Running iteration - [ {iteration + 1} / {runs} ]")
        logiteration = "%s_%i.log" % (logfilename.split(".")[0], iteration + 1)
        runtime = run_shell_command_time(f"{command}", logiteration, verbose)
        Messaging.step(f"Runtime = {runtime}")
        measures.append(runtime)

        # progress forward
        if progress_bar is not None:
            progress_bar.update(1)

    # build summary
    Messaging.step(f"All = {measures}")
    summary = {
        "results": [
            {
                "command": command,
                "mean": numpy.mean(measures),
                "median": numpy.median(measures),
                "stddev": numpy.std(measures),
                "min": numpy.min(measures),
                "max": numpy.max(measures),
                "times": measures,
            }
        ]
    }

    # return
    return summary


##########################################################
def native_hyperfine(command: str, runs: int = 2, verbose: bool = False) -> dict:
    with NamedTemporaryFile() as fp:
        hyper_command = (
            'hyperfine --export-json=%s --warmup 0 --output /tmp/croco.log --runs %s "%s"'
            % (fp.name, runs, command)
        )
        try:
            run_shell_command(hyper_command, capture=not verbose)
            result = json.load(fp)
            return result
        except Exception:
            with open("/tmp/croco.log", "r") as fp_log:
                content = fp_log.read()
                print(content)
            raise Exception("CROCO failed to run !")


##########################################################
def run_hyperfine(
    command: str,
    logfilename=None,
    progress_bar: tqdm = None,
    runs: int = 2,
    retries: int = 8,
    verbose: bool = False,
):
    hyperfine = False  # shutil.which("hyperfine")
    if hyperfine:
        return native_hyperfine(command, runs, verbose=verbose)
    else:
        return emulate_hyperfine(
            command, logfilename, progress_bar, runs, verbose=verbose
        )
