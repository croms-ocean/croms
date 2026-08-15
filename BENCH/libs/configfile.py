##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import json5
import yaml
import glob
from yaml import Loader
from jsonschema import validate
import copy
from .messaging import Messaging


##########################################################
class ConfigFile:
    def __init__(self, file_path: str):
        self.file_path = file_path

    def load(self) -> dict:
        self._load_config()
        return self.config

    @staticmethod
    def _apply_vars(value: str, var_name: str, var_value: str) -> str:
        return value.replace("{" + var_name + "}", str(var_value))

    @staticmethod
    def _tranverse_and_apply_vars(element, var_name: str, var_value):
        if isinstance(element, list):
            ConfigFile._process_list(element, var_name, var_value)
        elif isinstance(element, dict):
            ConfigFile._process_dict(element, var_name, var_value)
        else:
            raise Exception(f"Unsupported type in tree: {type(element)}")

    @staticmethod
    def _process_list(element, var_name, var_value):
        for key, entry in enumerate(element):
            if isinstance(entry, int) or entry is None:
                pass
            elif isinstance(entry, str):
                element[key] = ConfigFile._apply_vars(entry, var_name, var_value)
            elif isinstance(entry, (dict | list)):
                ConfigFile._tranverse_and_apply_vars(entry, var_name, var_value)
            else:
                raise Exception(f"Unsupported type in tree: {type(entry)}")

    @staticmethod
    def _process_dict(element, var_name, var_value):
        for key, entry in element.items():
            if isinstance(entry, int) or entry is None:
                pass
            elif isinstance(entry, str):
                element[key] = ConfigFile._apply_vars(entry, var_name, var_value)
            elif isinstance(entry, (dict | list)):
                ConfigFile._tranverse_and_apply_vars(entry, var_name, var_value)
            else:
                raise Exception(f"Unsupported type in tree: {type(entry)}")

    def _unpack_variant_vars(self) -> None:
        """Unpack the vars {xxx} in variant names"""
        variants = self.config["variants"]
        new_variants = {}
        to_delete = []
        for name, variant in variants.items():
            if "unpack" in variant:
                # check
                if len(variant["unpack"]) > 1:
                    raise Exception(
                        "Do not yet support multiple variable applying, need to think how to run over & mix !"
                    )

                # loop on vars
                for var_name, var_values in variant["unpack"].items():
                    for var_value in var_values:
                        # copy & make not anymore a template
                        variant_copy = copy.deepcopy(variant)
                        del variant_copy["unpack"]

                        # apply
                        ConfigFile._tranverse_and_apply_vars(
                            variant_copy, var_name, var_value
                        )

                        # calc new name
                        name_copy = ConfigFile._apply_vars(name, var_name, var_value)

                        # inject
                        new_variants[name_copy] = variant_copy

                # remove template
                to_delete.append(name)

        # inject & clean
        for name in to_delete:
            del variants[name]
        for key, value in new_variants.items():
            variants[key] = value

    def _load_config(self) -> None:
        """Load the config file and extract some infos"""
        # message
        Messaging.section("Loading configuration file")
        Messaging.step(f"Reading {self.file_path}")

        # load
        self.config = self.load_json_yaml_file(self.file_path)

        # handle imports
        self._handle_imports()

        # unpack variant templates
        self._unpack_variant_vars()

    def load_json_yaml_file(self, fname: str) -> dict:
        with open(fname, "r") as fp:
            if fname.endswith(".yaml"):
                return yaml.load(fp, Loader)
            else:
                return json5.load(fp)

    def _handle_imports(self):
        # loop on patterns
        for target, patterns in self.config["imports"].items():
            for pattern in patterns:
                # loop on found sub files
                for file in glob.glob(pattern):
                    Messaging.step(f"Importing {file}")
                    # load
                    data = self.load_json_yaml_file(file)

                    # merge in given target key
                    for key, value in data.items():
                        if self.config[target] is None:
                            self.config[target] = {}
                        if key in self.config[target]:
                            raise Exception(
                                f"Try to override an already existing key in {file} : {target}.{key}"
                            )
                        self.config[target][key] = value

        schema_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "schema.json"
        )
        with open(schema_path, "r") as fp:
            schema = json5.load(fp)
        validate(instance=self.config, schema=schema)
