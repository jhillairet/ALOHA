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
        Path to a scenario file (TOML, .m, or .mat format), or dictionary containing scenario inputs.

    """

    def __init__(self, scenario: str | os.PathLike | dict | None = None):

        self.scenario = {}
        self.results = {}
        if isinstance(scenario, (str, os.PathLike)):
            filepath = pathlib.Path(scenario)
            if filepath.suffix.lower() in (".m", ".mat"):
                # Handle MATLAB files directly in constructor
                if filepath.suffix.lower() == ".mat":
                    matlab_scenario, results_dict = load_mat_file(filepath)
                else:  # .m file
                    matlab_scenario, results_dict = load_m_file(filepath)

                # Convert the MATLAB scenario to the TOML schema
                self.scenario = self.convert_matlab_scenario(matlab_scenario)
                self.results = results_dict
            else:
                # Assume TOML file
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
                "antenne_elementaire": "8_active_waveguides.toml",
                "tutorial_aloha_antenna_simple_grill_8waveguides": "8_active_waveguides.toml",
                "antenna_8_active_waveguides": "8_active_waveguides.toml",
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

        # Add comment at the beginning of the file
        # Use empty string if no comment is found in the scenario
        comment_value = ""
        if "comment" in self.scenario:
            comment = self.scenario["comment"]
            if isinstance(comment, list):
                # Handle list of comments (MATLAB style) - convert all elements to strings
                non_empty_comments = [str(line) for line in comment if str(line).strip()]
                if non_empty_comments:
                    comment_value = " ".join(non_empty_comments)
            elif isinstance(comment, str) and comment.strip():
                comment_value = comment

        # Always add a comment at the beginning
        if comment_value.strip():
            toml_lines.append(f"# {comment_value}")
        else:
            toml_lines.append("# ALOHA scenario")
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

        This method executes the appropriate plasma coupling calculation based on the
        scenario configuration. Currently, it only supports spectral_1D solver with
        bilinear profile (version 6).
        """
        # Check if we have plasma configuration
        if "plasma" not in self.scenario:
            raise ValueError("Scenario must contain a 'plasma' section to run")

        plasma = self.scenario["plasma"]
        solver = plasma.get("solver", "")

        # Only proceed if solver is spectral_1D
        if solver == "spectral_1D":
            spectral_1D = plasma.get("spectral_1D", {})
            if not spectral_1D:
                raise ValueError("spectral_1D section is required in plasma for spectral_1D solver")

            profile = spectral_1D.get("profile", "")

            # Only run S_plasma_1D if profile is bilinear (version 6)
            if profile == "bilinear":
                from aloha.plasma import S_plasma_1D

                # Execute the plasma coupling calculation
                S_plasma, rac_Zhe = S_plasma_1D(self)

                # Store results in the scenario
                self.results["S_plasma"] = S_plasma
                self.results["rac_Zhe"] = rac_Zhe

                return
            else:
                raise ValueError(
                    f"Unsupported plasma profile '{profile}' for spectral_1D solver. Only 'bilinear' is supported."
                )
        else:
            raise ValueError(f"Unsupported solver '{solver}'. Only 'spectral_1D' is currently supported.")


def _convert_scenario_to_matlab_inputs(scenario: "Scenario") -> dict:
    """
    Convert a Scenario object from TOML schema to MATLAB-style parameter dictionary.

    This function maps the TOML schema parameters to the MATLAB-style parameters
    expected by the Fortran binary and S_plasma_1D_matlab_inputs function.
    """
    # Extract parameters from the scenario
    scenario_dict = scenario.scenario

    # Check that we have the required plasma section
    if "plasma" not in scenario_dict:
        raise ValueError("Scenario must contain a 'plasma' section")

    plasma = scenario_dict["plasma"]

    # Check solver type - we only support spectral_1D for now
    solver = plasma.get("solver", "")
    if solver != "spectral_1D":
        raise ValueError(f"_convert_scenario_to_matlab_inputs only supports 'spectral_1D' solver, got '{solver}'")

    spectral_1D = plasma.get("spectral_1D", {})
    if not spectral_1D:
        raise ValueError("spectral_1D section is required in plasma")

    # Extract antenna parameters
    if "antenna" not in scenario_dict:
        raise ValueError("Scenario must contain an 'antenna' section")

    antenna = scenario_dict["antenna"]
    excitation = antenna.get("excitation", {})

    # Load antenna file to get waveguide parameters
    antenna_file = antenna.get("file")
    antenna_data = {}
    if antenna_file:
        # Try to load the antenna file to get waveguide parameters
        try:
            from pathlib import Path

            from aloha.antenna import Antenna

            # Try to find the antenna file in the antennas directory
            # Try both .toml and .m extensions
            antenna_paths = [
                Path(antenna_file),  # Try as-is
                Path(__file__).parent.parent / "antennas" / antenna_file,  # Try in antennas directory
                Path(__file__).parent.parent.parent / "antennas" / antenna_file,  # Try in parent antennas directory
            ]

            # Also try with .m extension (for MATLAB antenna files)
            if not any(path.exists() for path in antenna_paths):
                antenna_paths.extend(
                    [
                        Path(f"{antenna_file}.m"),
                        Path(__file__).parent.parent / "antennas" / f"{antenna_file}.m",
                        Path(__file__).parent.parent.parent / "antennas" / f"{antenna_file}.m",
                    ]
                )

            for path in antenna_paths:
                if path.exists():
                    antenna_obj = Antenna.from_file(path)
                    antenna_data = antenna_obj.antenna
                    break
        except (FileNotFoundError, ImportError):
            # If we can't load the antenna file, we'll use default values
            pass

    # Get frequency from antenna excitation or antenna default
    freq = excitation.get("f", antenna.get("frequency", antenna_data.get("frequency", None)))
    if freq is None:
        raise ValueError("Frequency (f) is required in antenna.excitation or antenna")

    # Extract plasma profile parameters from spectral_1D section
    profile = spectral_1D.get("profile", "")
    if profile != "bilinear":
        raise ValueError(f"Unsupported plasma profile '{profile}'. Only 'bilinear' is supported.")

    bilinear = spectral_1D.get("bilinear", {})
    if not bilinear:
        raise ValueError("bilinear section is required in plasma.spectral_1D")

    # Extract bilinear profile parameters
    ne0 = bilinear.get("ne0", 0.0)  # edge density [1/m^3]
    lambda_n = bilinear.get("lambda_n", [0.002, 0.02])  # gradients scrape-off lengths [m]
    plasma_layer_length = bilinear.get("plasma_layer_length", 0.002)  # width of the first plasma layer [m]
    vacuum_layer_length = bilinear.get("vacuum_layer_length", 0.0)  # vacuum gap width [m]

    # For version 6, we need to map the bilinear profile to the linear profile parameters
    # Convert to arrays for multiple poloidal rows
    nb_g_pol = 1  # Default, will be calculated from antenna layout

    # Initialize antenna layout parameters with defaults
    nb_mod_phi = 1  # Number of modules in toroidal direction
    nb_mod_theta = 1  # Number of modules in poloidal direction

    # Extract antenna layout parameters
    layout = antenna_data.get("layout", {})
    if layout:
        nb_mod_phi = layout.get("nb_mod_phi", nb_mod_phi)
        nb_mod_theta = layout.get("nb_mod_theta", nb_mod_theta)
        nb_g_pol = nb_mod_theta  # Number of poloidal rows

    # Initialize module parameters with defaults
    nb_wg_theta = 1  # Number of waveguides per module in poloidal direction
    nb_wg_phi = 1  # Number of waveguides per module in toroidal direction
    mask = [1]  # Mask of active/passive waveguides
    nb_pwg_btw_mod_phi = 0  # Number of passive waveguides between modules
    nb_pwg_edge = 1  # Number of passive waveguides on each edge

    module = antenna_data.get("module", {})
    if module:
        nb_wg_theta = module.get("nb_wg_theta", nb_wg_theta)
        nb_wg_phi = module.get("nb_wg_phi", nb_wg_phi)
        mask = module.get("mask", mask)
        nb_pwg_btw_mod_phi = module.get("nb_pwg_btw_mod_phi", nb_pwg_btw_mod_phi)
        nb_pwg_edge = module.get("nb_pwg_edge", nb_pwg_edge)

    # Calculate total number of waveguides per poloidal row
    # Using the MATLAB formula: nb_g_total_ligne = nb_wg_phi * nb_mod_phi + 2 * nb_pwg_edge
    # + nb_pwg_btw_mod_phi * (nb_mod_phi - 1)
    # This matches the waveguide logic in aloha_utils_getAntennaCoordinates.m
    nb_g_total_ligne = nb_wg_phi * nb_mod_phi + 2 * nb_pwg_edge + nb_pwg_btw_mod_phi * (nb_mod_phi - 1)

    # Waveguide dimensions from antenna module parameters
    wg_size_theta = module.get("wg_size_theta", 70e-3)  # Height of waveguides in poloidal direction [m]
    awg_size_phi = module.get("awg_size_phi", 10e-3)  # Width of active waveguides [m]

    # For version 6, we need to provide arrays for each poloidal row
    # Convert scalar values to arrays with length nb_g_pol
    ne0_array = [ne0] * nb_g_pol
    dne0_array = [ne0 / lambda_n[0] if lambda_n and len(lambda_n) > 0 else 0.0] * nb_g_pol
    d_couche_array = [plasma_layer_length] * nb_g_pol
    # Fix dne1 calculation to match MATLAB: (1 + d_couche/lambda_n[0]) * ne0 / lambda_n[1]
    if lambda_n and len(lambda_n) > 1:
        dne1_array = [(1 + plasma_layer_length / lambda_n[0]) * ne0 / lambda_n[1]] * nb_g_pol
    else:
        dne1_array = [0.0] * nb_g_pol
    d_vide_array = [vacuum_layer_length] * nb_g_pol

    # Waveguide parameters using MATLAB waveguide logic from aloha_utils_getAntennaCoordinates.m
    # a = waveguide height in poloidal direction (constant for all waveguides in a line)
    a = wg_size_theta

    # b = array of waveguide widths in toroidal direction
    # For version 6, use active waveguide width for all waveguides
    b = [awg_size_phi] * nb_g_total_ligne

    # z = array of waveguide positions in toroidal direction
    # Calculate positions based on waveguide widths and spacing (e_phi)
    # From MATLAB: z(1,ind) = z(1,ind-1) + b(ind-1) + e_phi
    # Get the spacing between waveguides from antenna module parameters
    e_phi = module.get("e_phi", 1e-3)  # Default spacing if not specified
    z = [0.0] * nb_g_total_ligne
    for ind in range(1, nb_g_total_ligne):
        z[ind] = z[ind - 1] + b[ind - 1] + e_phi

    # Other parameters (using typical defaults for version 6)
    # Match MATLAB defaults from aloha_init.m
    T_grill = 7.0  # Grill periodicity parameter
    D_guide_max = 100.0  # Maximum guide decoupling distance [m]
    erreur_rel = 1e-6  # Relative error tolerance
    pertes = 1e-6  # Loss parameter (match MATLAB default from aloha_init.m)

    # Mode numbers (from spectral_1D or defaults)
    nb_evanescent_modes = spectral_1D.get("nb_evanescent_modes", 2)
    Nmh = 1  # Number of magnetic modes (typical default)
    Nme = nb_evanescent_modes  # Number of electric modes = nb_evanescent_modes

    # Build the MATLAB-style parameter dictionary
    matlab_params = {
        "antenna": {"freq": freq},
        "ne0": ne0_array,
        "dne0": dne0_array,
        "d_couche": d_couche_array,
        "dne1": dne1_array,
        "nb_g_pol": nb_g_pol,
        "nb_g_total_ligne": nb_g_total_ligne,
        "a": a,
        "b": b,
        "z": z,
        "T_grill": T_grill,
        "D_guide_max": D_guide_max,
        "erreur_rel": erreur_rel,
        "pertes": pertes,
        "d_vide": d_vide_array,
        "Nmh": Nmh,
        "Nme": Nme,
    }

    return matlab_params
