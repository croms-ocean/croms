##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import platform
import subprocess
import json
from datetime import datetime
import sys
import shutil  # To check if 'lscpu' exists


##########################################################
def run_and_get_first_line(command: str) -> str:
    """Run a shell command and return the first line of the output."""
    executable = command.split()[0]

    if shutil.which(executable) is None:
        return "NOT FOUND"

    try:
        output = subprocess.check_output(command, shell=True, text=True)
        return output.splitlines()[0] if output else "NO OUTPUT"
    except subprocess.CalledProcessError:
        return "UNAVAILABLE"
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return "ERROR"


##########################################################
def get_lscpu_as_json() -> dict:
    """Get lscpu information in JSON format, with a warning if lscpu is not available."""
    if shutil.which("lscpu") is None:
        # If 'lscpu' is not available, issue a warning and return an empty dict
        print(
            "Warning: 'lscpu' command not found. Proceeding without it.",
            file=sys.stderr,
        )
        return {}

    try:
        # Try to execute lscpu command and parse output as JSON
        data = subprocess.check_output(["lscpu", "-J"])
        return json.loads(data)
    except subprocess.CalledProcessError as e:
        # If lscpu fails, log a warning and return an empty dictionary
        print(f"Warning: Failed to run 'lscpu'. Error: {e}", file=sys.stderr)
        return {}


##########################################################
def gen_system_info() -> dict:
    """Generate system information."""
    infos = {}

    # Date and Time
    infos["date"] = datetime.now().strftime("%d/%m/%y")
    infos["time"] = datetime.now().strftime("%H:%M:%S")

    # Platform Information
    system = infos["plateform"] = {}
    system["hostname"] = platform.node()
    system["arch"] = platform.architecture()[
        0
    ]  # Just the architecture (32bit or 64bit)
    system["system"] = platform.system()
    system["uname"] = platform.uname()

    # Get lscpu data if available
    lscpu = get_lscpu_as_json().get("lscpu", [])
    infos["lscpu"] = lscpu

    # Get Processor Name from /proc/cpuinfo
    processor_name = subprocess.getoutput(
        "cat /proc/cpuinfo | grep 'model name' | cut -d ':' -f 2- | head -n 1 | xargs echo"
    )
    system["processor_name"] = processor_name.strip()

    # Extract processor data from lscpu if available
    for entry in lscpu:
        if "children" in entry:
            system["processor"] = entry["children"][0]["data"]

    # Software Information (MPI, GCC, etc.)
    soft = infos["software"] = {}
    soft["mpi"] = run_and_get_first_line("mpirun --version")
    soft["gcc"] = run_and_get_first_line("gcc --version")
    soft["gfortran"] = run_and_get_first_line("gfortran --version")
    soft["nvfortran"] = run_and_get_first_line("nvfortran --version")

    return infos
