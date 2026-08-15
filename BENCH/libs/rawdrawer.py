##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
from typing import List

# external libs
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

# internal
from .messaging import Messaging
from .config import Config


##########################################################
class RawDrawer:
    def __init__(self, netcdf_file: str, config: Config, case_name: str) -> None:
        self.netcdf_file = netcdf_file
        self.config = config
        self.case_name = case_name
        self.data = Dataset(netcdf_file, "r")

    def gen_render_mode(self, variable_data, variable: str) -> str:
        # vars
        shape = variable_data.shape

        # mode
        if len(shape) == 4:
            return "2D"
        elif len(shape) <= 3:
            return "SKIP"
        else:
            raise Exception(f"Unsupported shape for {variable} : {shape}")

    def render_data(self, variable_data, step: int, variable: str, ax) -> tuple:
        # dimensions
        cs = None
        mode = self.gen_render_mode(variable_data, variable)
        if mode == "2D":
            cs = ax.imshow(
                variable_data[step, -1, :, :], interpolation="bilinear"
            )  # surface velocity  (along x)  last time step
        elif mode == "SKIP":
            pass
        else:
            raise Exception(f"Unsupported mode for {variable} : {mode}")

        # render
        return cs

    def plot_variable(self, variable: str, out_image: str, step: int = 0) -> None:
        # get variable
        variable_data = self.data.variables[variable]
        shape = variable_data.shape

        # skip
        if len(shape) < 3:
            return

        # case middle time step
        if step == -2:
            step = shape[0] / 2

        # log
        Messaging.step(f"Render var={variable}, step={step}, shape={shape}")

        # build plot
        fig, ax = plt.subplots()
        self.set_ax_labels(ax, variable)

        # plit
        cs = self.render_data(variable_data, step, variable, ax)

        # render
        if cs:
            os.makedirs(os.path.dirname(out_image), exist_ok=True)
            fig.colorbar(cs)
            plt.savefig(out_image, dpi=600)

        # close
        plt.close(fig)

    def plot_variables(self, variables: list, basename: str, step: int = 0) -> None:
        # case
        case_name = self.case_name

        # fix step name
        if step == 0:
            step_name = "t-start"
        elif step == -1:
            step_name = "t-end"
        elif step == -2:
            step_name = "t-middle"
        else:
            step_name = f"t-{step}"

        # loop on vars
        for vname in variables:
            oname = f"{basename}/{vname}/{case_name}-{vname}-{step_name}"
            self.plot_variable(vname, oname, step)

    def plot_all_variables(self, basename: str, steps: List[int] = [0, -2, -1]) -> None:
        for step in steps:
            self.plot_variables(self.data.variables.keys(), basename, step)

    def set_ax_labels(self, ax, variable: str) -> None:
        # case
        case_name = self.case_name
        cases = self.config.config["cases"]

        # extract title
        case_config = cases[case_name]
        title = case_config.get("title", case_name)

        # set labels
        ax.set_ylabel("y")
        ax.set_ylabel("x")
        ax.set_title(f"CROCO - {title} - {variable}")

    def plot_variable_animation(self, variable: str, basename: str) -> None:
        # get variable
        variable_data = self.data.variables[variable]
        shape = variable_data.shape

        # check OK
        mode = self.gen_render_mode(variable_data, variable)
        if mode == "SKIP":
            return

        # plot vars
        artists = []
        fig, ax = plt.subplots()
        self.set_ax_labels(ax, variable)

        # progress over time steps
        for i in range(shape[0]):
            container = self.render_data(
                self.data.variables[variable], i, variable, ax=ax
            )
            artists.append([container])

        fig.colorbar(container, ax=ax)

        ani = ArtistAnimation(fig=fig, artists=artists, interval=40, repeat=True)
        ani.save(f"{basename}.gif", writer="pillow")

        # # ffmpeg has to be install, not listed in dependencies
        # ani.save(f"{basename}.mp4", writer="ffmpeg")

        # close
        plt.close(fig)

    def plot_variables_animation(self, variables: list, basename: str) -> None:
        # loop on vars
        for vname in variables:
            oname = f"{basename}/{vname}"
            self.plot_variable_animation(vname, oname)

    def plot_all_variables_animation(self, basename: str) -> None:
        self.plot_variables_animation(self.data.variables.keys(), basename)
