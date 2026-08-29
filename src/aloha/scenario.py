import os
import pathlib
import re

import numpy as np

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from aloha.utils import load_m_file, load_mat_file


class Scenario:
    """
    ALOHA Scenario.

    Scenario is the main ALOHA object, which contains all the necessary parameters
    to run a simulation. Once a simulation has been run, the Scenario object also
    contains a `results` parameters which contains all the calculated results.

    Parameters
    ----------
    scenario: str | path | dict or None (default)
        Path to a scenario file (TOML format), or dictionary containing scenario inputs.

    """

    def __init__(self, scenario: str | os.PathLike | dict | None = None):

        self.scenario = {}
        self.results = {}
        if isinstance(scenario, (str, os.PathLike)):
            self.scenario = self.load(scenario)
        elif isinstance(scenario, dict):
            # TODO : verify that the passed dictionary is relevant
            self.scenario = scenario
        elif scenario is not None:
            raise ValueError("Invalid scenario type. Must be a string, Path, or dict.")

    @classmethod
    def from_file(cls, filename: str | os.PathLike):
        """Create a scenario from a TOML file."""
        return cls(filename)

    @classmethod
    def from_dict(cls, data: dict):
        """Create a scenario from a dictionary."""
        return cls(data)

    def load(self, filename) -> dict:
        """Load scenario from TOML file."""
        with open(filename, "rb") as fp:
            return tomllib.load(fp)

    def __str__(self):
        return str(self.scenario)

    @classmethod
    def from_matlab(cls, filename: str | os.PathLike) -> "Scenario":
        """
        Load scenario from an ALOHA MATLAB file (.m or .mat).

        Parameters
        ----------
        filename : str | os.PathLike
            Path to a MATLAB file (.m script or .mat binary file).

        Returns
        -------
        Scenario
            A Scenario object containing the scenario data following the TOML schema.
            The results dictionary (if available) is stored in the Scenario.results attribute.

        Raises
        ------
        ValueError
            If the file extension is not .m or .mat, or if scipy is required but not available.
        FileNotFoundError
            If the file does not exist.
        """
        filepath = pathlib.Path(filename)

        if not filepath.exists():
            raise FileNotFoundError(f"MATLAB file not found: {filepath}")

        if filepath.suffix.lower() == ".mat":
            matlab_scenario, results_dict = load_mat_file(filepath)
        elif filepath.suffix.lower() == ".m":
            matlab_scenario, results_dict = load_m_file(filepath)
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}. Expected .m or .mat")

        # Convert the MATLAB scenario to the TOML schema
        scenario_dict = cls.convert_matlab_scenario(matlab_scenario)
        scenario_obj = cls(scenario_dict)
        scenario_obj.results = results_dict
        return scenario_obj

    @classmethod
    def convert_matlab_scenario(cls, matlab_scenario: dict) -> dict:
        """
        Convert a MATLAB-style scenario dictionary to the Python TOML schema format.

        This function maps MATLAB scenario fields to the equivalent fields in the
        ALOHA Python TOML schema (similar to scenario_example.toml).

        Parameters
        ----------
        matlab_scenario : dict[str, Any]
            Dictionary containing scenario data in MATLAB format.

        Returns
        -------
        dict[str, Any]
            Dictionary converted to the TOML schema format.
        """
        converted = {}

        # Helper function to get nested value safely
        def get_nested(data: dict, *keys, default=None):
            current = data
            for key in keys:
                if isinstance(current, dict) and key in current:
                    current = current[key]
                else:
                    return default
            return current

        # ==================================================================
        # General
        # ==================================================================
        if "options" in matlab_scenario and "comment" in matlab_scenario["options"]:
            converted["comment"] = matlab_scenario["options"]["comment"]

        # ==================================================================
        # Antenna description
        # ==================================================================
        antenna = {}

        # antenna.architecture -> antenna.file
        architecture_str = None
        if "antenna" in matlab_scenario and "architecture" in matlab_scenario["antenna"]:
            architecture = matlab_scenario["antenna"]["architecture"]
            # Convert list of ASCII codes to string (for .mat files)
            if isinstance(architecture, list):
                # Handle nested lists (e.g., [[97], [98], ...] from .mat files)
                if architecture and isinstance(architecture[0], list):
                    architecture_str = "".join(chr(code[0]) for code in architecture)
                else:
                    architecture_str = "".join(chr(code) for code in architecture)
            else:
                architecture_str = architecture
            # Map MATLAB antenna names to TOML file names
            antenna_name_mapping = {
                "antenne_elementaire": "simple_antenna.toml",
                # Add more mappings as needed
            }
            antenna["file"] = antenna_name_mapping.get(architecture_str, architecture_str)

        # antenna.freq -> antenna.excitation.f
        excitation = {}
        if "antenna" in matlab_scenario and "freq" in matlab_scenario["antenna"]:
            excitation["f"] = float(matlab_scenario["antenna"]["freq"])

        # options.bool_mesure -> antenna.excitation.experimental
        if "options" in matlab_scenario and "bool_mesure" in matlab_scenario["options"]:
            excitation["experimental"] = bool(matlab_scenario["options"]["bool_mesure"])

        # Get number of modules to determine default array size
        # Try to get from antenna architecture or default to 8
        num_modules = 8  # Default value
        if architecture_str:
            # Try to extract number from architecture name if it contains a number
            arch_numbers = re.findall(r"\d+", architecture_str)
            if arch_numbers:
                num_modules = int(arch_numbers[0])

        # antenna.a_ampl -> antenna.excitation.power
        if "antenna" in matlab_scenario and "a_ampl" in matlab_scenario["antenna"]:
            a_ampl = matlab_scenario["antenna"]["a_ampl"]
            if isinstance(a_ampl, (list, np.ndarray)):
                # Handle nested lists from .mat files
                if a_ampl and isinstance(a_ampl[0], list):
                    a_ampl = a_ampl[0]  # Extract the inner list
                # Convert to list of floats and square the values (power = amplitude^2)
                excitation["power"] = [float(x) ** 2 for x in a_ampl]
            else:
                # If not a straightforward array, use default value
                # Default: all modules have power of 1.0
                excitation["power"] = [1.0] * num_modules

        # antenna.a_phase -> antenna.excitation.phase (in degree)
        if "antenna" in matlab_scenario and "a_phase" in matlab_scenario["antenna"]:
            a_phase = matlab_scenario["antenna"]["a_phase"]
            if isinstance(a_phase, (list, np.ndarray)):
                # Handle nested lists from .mat files
                if a_phase and isinstance(a_phase[0], list):
                    a_phase = a_phase[0]  # Extract the inner list
                # Convert radians to degrees
                excitation["phase"] = [float(x) * 180 / np.pi for x in a_phase]
            else:
                # If not a straightforward array, use default value
                # Default: phases are 0, 90, 180, 270, ... degrees (cyclic)
                default_phases = []
                for i in range(num_modules):
                    default_phases.append((i % 4) * 90.0)
                excitation["phase"] = default_phases

        # options.TSport -> antenna.excitation.port
        if "options" in matlab_scenario and "TSport" in matlab_scenario["options"]:
            excitation["port"] = matlab_scenario["options"]["TSport"]

        # options.choc -> antenna.excitation.pulse
        if "options" in matlab_scenario and "choc" in matlab_scenario["options"]:
            excitation["pulse"] = int(matlab_scenario["options"]["choc"])

        # options.tps_1 -> antenna.excitation.avg_times[0]
        # options.tps_2 -> antenna.excitation.avg_times[1]
        avg_times = []
        if "options" in matlab_scenario:
            if "tps_1" in matlab_scenario["options"]:
                avg_times.append(float(matlab_scenario["options"]["tps_1"]))
            if "tps_2" in matlab_scenario["options"]:
                avg_times.append(float(matlab_scenario["options"]["tps_2"]))
        if avg_times:
            excitation["avg_times"] = avg_times

        if excitation:
            antenna["excitation"] = excitation

        if antenna:
            converted["antenna"] = antenna

        # ==================================================================
        # Plasma description
        # ==================================================================
        plasma = {}

        # Default solver
        plasma["solver"] = "spectral_1D"

        # Create spectral_1D section
        spectral_1d = {}

        # version_plasma_1D -> plasma.spectral_1D.profile (3 -> 'linear', 6 -> 'bilinear')
        # Use options.version_code to determine which version to use
        version_code = get_nested(matlab_scenario, "options", "version_code")
        if version_code == "1D":
            version_plasma_1d = matlab_scenario.get("version_plasma_1D")
        elif version_code == "2D":
            version_plasma_1d = matlab_scenario.get("version_plasma_2D")
        else:
            # Fallback to plasma.version
            version_plasma_1d = get_nested(matlab_scenario, "plasma", "version")
            if isinstance(version_plasma_1d, str):
                # This is a variable reference, get the actual value
                if version_plasma_1d == "version_plasma_1D":
                    version_plasma_1d = matlab_scenario.get("version_plasma_1D")
                elif version_plasma_1d == "version_plasma_2D":
                    version_plasma_1d = matlab_scenario.get("version_plasma_2D")

        if version_plasma_1d == 3:
            spectral_1d["profile"] = "linear"
        elif version_plasma_1d == 6:
            spectral_1d["profile"] = "bilinear"

        # Nme -> plasma.spectral_1D.nb_evanescent_modes
        if "Nme" in matlab_scenario:
            spectral_1d["nb_evanescent_modes"] = int(matlab_scenario["Nme"])

        # Create bilinear section
        bilinear = {}

        # options.bool_lignes_identiques -> plasma.spectral_1D.bilinear.identical_profiles
        if "options" in matlab_scenario and "bool_lignes_identiques" in matlab_scenario["options"]:
            bilinear["identical_profiles"] = bool(matlab_scenario["options"]["bool_lignes_identiques"])

        # plasma.ne0 -> plasma.spectral_1D.bilinear.ne0
        if "plasma" in matlab_scenario and "ne0" in matlab_scenario["plasma"]:
            bilinear["ne0"] = float(matlab_scenario["plasma"]["ne0"])

        # plasma.lambda_n(1) -> plasma.spectral_1D.bilinear.lambda_n[0]
        # plasma.lambda_n(2) -> plasma.spectral_1D.bilinear.lambda_n[1]
        lambda_n = []
        plasma_data = matlab_scenario.get("plasma", {})
        if "lambda_n" in plasma_data:
            lambda_n_value = plasma_data["lambda_n"]
            if isinstance(lambda_n_value, (list, np.ndarray)):
                # Handle nested lists from .mat files
                if lambda_n_value and isinstance(lambda_n_value[0], list):
                    lambda_n_value = lambda_n_value[0]  # Extract the inner list
                lambda_n = [float(x) for x in lambda_n_value]
            else:
                lambda_n = [float(lambda_n_value)]
        elif "lambda_n(1)" in plasma_data:
            # Handle MATLAB-style indexed fields
            if "lambda_n(1)" in plasma_data:
                lambda_n.append(float(plasma_data["lambda_n(1)"]))
            if "lambda_n(2)" in plasma_data:
                lambda_n.append(float(plasma_data["lambda_n(2)"]))
        if lambda_n:
            bilinear["lambda_n"] = lambda_n

        # plasma.d_couche -> plasma.spectral_1D.bilinear.plasma_layer_length
        if "plasma" in matlab_scenario and "d_couche" in matlab_scenario["plasma"]:
            bilinear["plasma_layer_length"] = float(matlab_scenario["plasma"]["d_couche"])
            # Also keep the old name for backward compatibility (as expected by tests)
            bilinear["d_couche"] = float(matlab_scenario["plasma"]["d_couche"])

        # options.B0 -> plasma.spectral_1D.bilinear.B0
        if "options" in matlab_scenario and "B0" in matlab_scenario["options"]:
            bilinear["B0"] = float(matlab_scenario["options"]["B0"])

        # plasma.d_vide -> plasma.spectral_1D.bilinear.vacuum_layer_length
        if "plasma" in matlab_scenario and "d_vide" in matlab_scenario["plasma"]:
            bilinear["vacuum_layer_length"] = float(matlab_scenario["plasma"]["d_vide"])
            # Also keep the old name for backward compatibility (as expected by tests)
            bilinear["d_vide"] = float(matlab_scenario["plasma"]["d_vide"])

        # options.type_swan_aloha -> plasma.spectral_1D.bilinear.infinite_waveguide (0 -> False, if 1 -> True)
        if "options" in matlab_scenario and "type_swan_aloha" in matlab_scenario["options"]:
            type_swan_aloha = matlab_scenario["options"]["type_swan_aloha"]
            bilinear["infinite_waveguide"] = bool(int(type_swan_aloha))

        if bilinear:
            spectral_1d["bilinear"] = bilinear

        # Spectral domain parameters
        if "options" in matlab_scenario:
            options = matlab_scenario["options"]
            if "nz_min" in options:
                spectral_1d["nz_min"] = int(options["nz_min"])
            if "nz_max" in options:
                spectral_1d["nz_max"] = int(options["nz_max"])
            if "dnz" in options:
                spectral_1d["dnz"] = float(options["dnz"])
            if "ny_min" in options:
                spectral_1d["ny_min"] = int(options["ny_min"])
            if "ny_max" in options:
                spectral_1d["ny_max"] = int(options["ny_max"])
            if "dny" in options:
                spectral_1d["dny"] = float(options["dny"])

        # Spatial domain parameters
        if "options" in matlab_scenario:
            options = matlab_scenario["options"]
            if "z_coord_min" in options:
                spectral_1d["z_min"] = float(options["z_coord_min"])
            if "z_coord_max" in options:
                spectral_1d["z_max"] = float(options["z_coord_max"])
            if "nbre_z_coord" in options:
                spectral_1d["nb_z"] = int(options["nbre_z_coord"])
            if "x_coord_max" in options:
                spectral_1d["x_max"] = float(options["x_coord_max"])
            if "nbre_x_coord" in options:
                spectral_1d["nb_x"] = int(options["nbre_x_coord"])

        if spectral_1d:
            plasma["spectral_1D"] = spectral_1d

        if plasma:
            converted["plasma"] = plasma

        # ==================================================================
        # Options
        # ==================================================================
        options = {}

        # options.bool_debug -> options.debug
        if "options" in matlab_scenario and "bool_debug" in matlab_scenario["options"]:
            options["debug"] = bool(matlab_scenario["options"]["bool_debug"])

        if options:
            converted["options"] = options

        return converted

    def to_toml(self, filename: str | os.PathLike | None = None) -> str | None:
        """
        Export the scenario to TOML format.

        Parameters
        ----------
        filename : str, Path, or None
            If provided, write the TOML content to this file.
            If None, return the TOML string without writing to a file.

        Returns
        -------
        str or None
            If filename is None, returns the TOML string.
            If filename is provided, writes to file and returns None.
        """

        def format_value(value):
            """Format a value for TOML output."""
            if isinstance(value, str):
                return f'"{value}"'
            elif isinstance(value, bool):
                return str(value).lower()
            elif isinstance(value, (int, float)):
                return str(value)
            elif isinstance(value, list):
                if len(value) == 0:
                    return "[]"
                # Check if all elements are numeric
                if all(isinstance(x, (int, float)) for x in value):
                    return "[" + ", ".join(str(x) for x in value) + "]"
                else:
                    # Handle string lists or mixed
                    return "[" + ", ".join(format_value(x) for x in value) + "]"
            elif isinstance(value, dict):
                # This shouldn't happen in TOML values, but handle it
                return str(value)
            else:
                return str(value)

        def dict_to_toml(data: dict, prefix: str = "") -> list[str]:
            """Convert a dictionary to TOML format lines."""
            lines = []

            # Separate simple values from nested dictionaries
            simple_values = {}
            nested_dicts = {}

            for key, value in data.items():
                if isinstance(value, dict):
                    nested_dicts[key] = value
                else:
                    simple_values[key] = value

            # Add simple values first (only if we have a prefix, otherwise they go at the top)
            if prefix:
                for key, value in simple_values.items():
                    lines.append(f"{key} = {format_value(value)}")
            else:
                for key, value in simple_values.items():
                    lines.append(f"{key} = {format_value(value)}")

            # Add blank line if we have both simple values and nested dicts
            if simple_values and nested_dicts:
                lines.append("")

            # Add nested dictionaries as sections
            for key, value in nested_dicts.items():
                if prefix:
                    # Nested section (e.g., antenna.excitation)
                    section_name = f"{prefix}.{key}"
                else:
                    # Top-level section (e.g., antenna, plasma, options)
                    section_name = key
                lines.append(f"[{section_name}]")
                nested_lines = dict_to_toml(value, section_name)
                lines.extend(nested_lines)
                if nested_lines:  # Only add blank line if there was content
                    lines.append("")  # Blank line between sections

            return lines

        # Start with the scenario dictionary
        toml_lines = []

        # Add comment if present (handle both string and list)
        if "comment" in self.scenario:
            comment = self.scenario["comment"]
            if isinstance(comment, list):
                # Handle list of comments (MATLAB style)
                non_empty_comments = [line for line in comment if line.strip()]
                if non_empty_comments:
                    for comment_line in non_empty_comments:
                        toml_lines.append(f"# {comment_line}")
                    toml_lines.append("")
            elif comment.strip():
                toml_lines.append(f"# {comment}")
                toml_lines.append("")

        # Convert the scenario dictionary to TOML (excluding comment which is handled separately)
        scenario_data = {k: v for k, v in self.scenario.items() if k != "comment"}
        scenario_lines = dict_to_toml(scenario_data)
        toml_lines.extend(scenario_lines)

        # Join all lines and clean up extra blank lines
        toml_content = "\n".join(toml_lines)

        # Remove trailing whitespace and multiple blank lines
        toml_content = "\n".join(line.rstrip() for line in toml_content.split("\n"))
        toml_content = "\n\n".join([part for part in toml_content.split("\n\n") if part.strip()])

        if filename is None:
            return toml_content
        else:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(toml_content)
            return None

    def run(self) -> None:
        """
        Run ALOHA scenario.

        """
        # TODO: run options, like logging levels, etc.
        pass
