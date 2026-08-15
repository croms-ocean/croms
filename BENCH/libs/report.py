##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json
import xml.etree.ElementTree as ET

# termcolor
from termcolor import cprint


##########################################################
class Report:
    def __init__(self, config) -> None:
        self.config = config
        self.report = {}

    def report_status(
        self,
        case: str,
        variant: str,
        restarted: bool,
        step: str,
        status: bool,
        message: str = "",
        log: str = "",
    ) -> None:
        if restarted:
            full_name = f"{case}-{variant}-rst"
        else:
            full_name = f"{case}-{variant}"
        if full_name not in self.report:
            self.report[full_name] = {
                "case": case,
                "variant": variant,
                "restarted": restarted,
                "status": {
                    "build": {"status": None},
                    "run": {"status": None},
                    "check": {"status": None},
                    "plotphy": {"status": None},
                },
            }
        self.report[full_name]["status"][step] = {
            "status": status,
            "message": message,
            "log": log,
        }

    def save(self) -> None:
        status_report_file = os.path.join(self.config.results, "status_report.json")
        with open(status_report_file, "w+") as fp:
            json.dump(self.report, fp, indent="\t")

    # To detect if report contains False values
    def contains_false(self):
        has_false = False
        for key_case, value_case in self.report.items():
            status_dict = value_case.get("status")
            for key, value in status_dict.items():
                if isinstance(value, dict) and value.get("status") is False:
                    has_false = True
        return has_false

    def get_color(self, status):
        if status is True:
            return "green"
        elif status is False:
            return "red"
        elif status is None:
            return "cyan"
        else:
            return "white"

    def print(self) -> None:
        print(
            "------------------------------- REPORT SUMMARY ------------------------------------"
        )

        # header
        full_name = "CASE-VARIANT"
        status_build = "Build"
        status_run = "Run"
        status_check = "Check"
        status_plotphy = "Plotphy"
        print(
            "%s %s %s %s %s"
            % (
                f"{full_name:<40}",
                f"{str(status_build):<8}",
                f"{str(status_run):<8}",
                f"{str(status_check)}",
                f"{str(status_plotphy)}",
            )
        )
        print(
            "..................................................................................."
        )

        for full_name, infos in self.report.items():
            status_build = infos["status"]["build"]["status"]
            status_run = infos["status"]["run"]["status"]
            status_check = infos["status"]["check"]["status"]
            status_plotphy = infos["status"]["plotphy"]["status"]

            cprint(f"{full_name:<40} ", "white", end="")
            cprint(f"{str(status_build):<8} ", self.get_color(status_build), end="")
            cprint(f"{str(status_run):<8} ", self.get_color(status_run), end="")
            cprint(f"{str(status_check):<8} ", self.get_color(status_check), end="")
            cprint(f"{str(status_plotphy):<8}", self.get_color(status_plotphy))

        print(
            "-----------------------------------------------------------------------------------"
        )

    def save_to_junit(self):
        junit_file = os.path.join(self.config.results, "junit_report.xml")

        # Creation of JUnit XML structure
        testsuites = ET.Element("testsuites")
        testsuite = ET.SubElement(testsuites, "testsuite", name="CrocoTestSuite")

        total_tests = 0
        total_failures = 0
        total_errors = 0

        for test_name, test_data in self.report.items():
            for step, result in test_data["status"].items():
                if result is None or not isinstance(result, dict):
                    continue

                if result["status"] is None:  # if skip step
                    continue

                total_tests += 1
                testcase = ET.SubElement(
                    testsuite, "testcase", name=f"{test_name}-{step}"
                )

                if result["status"] is False:  # if failed step
                    total_failures += 1
                    failure = ET.SubElement(
                        testcase,
                        "failure",
                        message=result.get("message", "Test failed"),
                    )
                    failure.text = result.get("log", "")

                    if "exception" in result.get("message", "").lower():
                        total_errors += 1  # count critical errors

        # Add stats
        testsuite.set("tests", str(total_tests))
        testsuite.set("failures", str(total_failures))
        testsuite.set("errors", str(total_errors))

        # Save XML file
        tree = ET.ElementTree(testsuites)
        tree.write(junit_file, encoding="utf-8", xml_declaration=True)
