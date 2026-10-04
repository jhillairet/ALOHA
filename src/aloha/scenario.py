import os
import pathlib
import re

import numpy as np

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from pathlib import Path
from pprint import pp

from aloha.antenna import Antenna
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

        # Safely set debug flag, defaulting to False if not present
        self.debug = self.scenario.get("options", {}).get("debug", False)

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

    def __eq__(self, other):
        """
        Compare two Scenario objects for equality.

        Parameters
        ----------
        other : Scenario
            The other Scenario object to compare with.

        Returns
        -------
        bool
            True if the scenarios are equal (within numerical tolerance), False otherwise.
        """
        if not isinstance(other, Scenario):
            return False
        try:
            return Scenario.compare_scenarios(self.scenario, other.scenario)
        except AssertionError:
            return False

    @staticmethod
    def compare_scenarios(scenario1, scenario2, path="", _differences=None):
        """
        Recursively compare two scenario dictionaries with NumPy array support.

        Parameters
        ----------
        scenario1 : dict
            First scenario dictionary to compare.
        scenario2 : dict
            Second scenario dictionary to compare.
        path : str, optional
            Current path in the dictionary hierarchy (for error messages).
        _differences : list, optional
            Internal parameter to collect all differences (do not use directly).

        Returns
        -------
        bool
            True if the scenarios are equal (within numerical tolerance), False otherwise.

        Raises
        ------
        AssertionError
            If the scenarios are not equal, with detailed error message listing all differences.
        """
        # Initialize differences list if not provided
        if _differences is None:
            _differences = []

        if isinstance(scenario1, dict) and isinstance(scenario2, dict):
            # Check that both have the same keys
            keys1 = set(scenario1.keys())
            keys2 = set(scenario2.keys())
            missing_in_2 = keys1 - keys2
            missing_in_1 = keys2 - keys1

            if missing_in_2 or missing_in_1:
                msg_parts = []
                if missing_in_2:
                    msg_parts.append(f"keys {sorted(missing_in_2)} missing in scenario2")
                if missing_in_1:
                    msg_parts.append(f"keys {sorted(missing_in_1)} missing in scenario1")
                msg = f"Keys mismatch at '{path}': {', '.join(msg_parts)}"
                _differences.append(msg)
                # Continue to collect all differences, don't raise yet
                return False

            # Recursively compare each value
            for key in scenario1:
                new_path = f"{path}.{key}" if path else key
                # Skip debug option as it may differ between file formats
                if new_path == "options.debug":
                    continue
                Scenario.compare_scenarios(scenario1[key], scenario2[key], new_path, _differences)

        elif isinstance(scenario1, np.ndarray) and isinstance(scenario2, np.ndarray):
            # Compare NumPy arrays with tolerance
            try:
                np.testing.assert_allclose(scenario1, scenario2, rtol=1e-10, atol=1e-10)
            except AssertionError as e:
                msg = f"Array mismatch at '{path}': {str(e).strip()}"
                _differences.append(msg)
                return False

        elif isinstance(scenario1, (list, tuple)) and isinstance(scenario2, (list, tuple)):
            # Check length first
            if len(scenario1) != len(scenario2):
                msg = f"List length mismatch at '{path}': {len(scenario1)} vs {len(scenario2)}"
                _differences.append(msg)
                return False

            # Recursively compare each element
            for i, (item1, item2) in enumerate(zip(scenario1, scenario2, strict=False)):
                new_path = f"{path}[{i}]"
                Scenario.compare_scenarios(item1, item2, new_path, _differences)

        else:
            # Compare scalar values directly
            if scenario1 != scenario2:
                msg = f"Value mismatch at '{path}': {scenario1} != {scenario2}"
                _differences.append(msg)
                return False

        # If we're at the top level and have differences, raise them all
        if path == "" and _differences:
            header = f"Found {len(_differences)} difference(s) between scenarios:"
            all_differences = header + "\n" + "\n".join(f"  - {diff}" for diff in _differences)
            raise AssertionError(all_differences)

        return True

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
            elif isinstance(architecture, np.ndarray):
                # Handle numpy arrays (from .mat files with preserved arrays)
                if architecture.ndim == 2 and architecture.shape[1] == 1:
                    # Column vector of character codes
                    architecture_str = "".join(chr(code) for code in architecture.flatten())
                else:
                    architecture_str = "".join(chr(code) for code in architecture.flatten())
            else:
                architecture_str = architecture
            # Map MATLAB antenna names to TOML file names
            antenna_name_mapping = {
                "antenne_elementaire": "8_active_waveguides.toml",
                "tutorial_aloha_antenna_simple_grill_8waveguides": "8_active_waveguides.toml",
                "antenna_8_active_waveguides": "8_active_waveguides.toml",
                "antenna_C3_ITM": "WEST_LH1_top.toml",
                "antenna_C4_ITM": "WEST_LH2_top.toml",
                "antenna_PA_ITM": "PA_1row.toml",
                # Add more mappings as needed
            }
            antenna["file"] = antenna_name_mapping.get(architecture_str, architecture_str)

        # antenna.freq -> antenna.excitation.f
        excitation = {}
        if "antenna" in matlab_scenario and "freq" in matlab_scenario["antenna"]:
            excitation["f"] = float(matlab_scenario["antenna"]["freq"])

        # options.bool_mesure -> antenna.excitation.experimental
        if "options" in matlab_scenario and "bool_mesure" in matlab_scenario["options"]:
            bool_mesure = matlab_scenario["options"]["bool_mesure"]
            # Handle MATLAB boolean strings properly ("true"/"false")
            if isinstance(bool_mesure, str):
                excitation["experimental"] = bool_mesure.lower() == "true"
            else:
                excitation["experimental"] = bool(bool_mesure)

        # Get number of modules to determine default array size
        # Try to get from antenna architecture
        num_modules = None
        if architecture_str:
            # Try to extract number from architecture name if it contains a number
            arch_numbers = re.findall(r"\d+", architecture_str)
            if arch_numbers:
                num_modules = int(arch_numbers[0])
            else:
                # If no number can be extracted, check if a_ampl or a_phase exists to infer the count
                if "antenna" in matlab_scenario:
                    if "a_ampl" in matlab_scenario["antenna"]:
                        a_ampl = matlab_scenario["antenna"]["a_ampl"]
                        if isinstance(a_ampl, (list, np.ndarray)):
                            num_modules = len(a_ampl)
                        elif isinstance(a_ampl, str):
                            # Try to extract the number from MATLAB expressions like ones(8,1) or (0:7)
                            if ":" in a_ampl:
                                # Handle colon operator (e.g., 0:7 -> 8 elements)
                                colon_parts = a_ampl.split(":")
                                if len(colon_parts) >= 2:
                                    start = int(colon_parts[0].strip())
                                    end = int(colon_parts[1].strip().rstrip(")' "))
                                    num_modules = end - start + 1
                            else:
                                # For expressions like sqrt(1)*ones(8,1), extract the number from ones(8,1)
                                # Use regex to find the last occurrence of ones(N,1) or zeros(N,1)
                                ones_match = re.search(r"ones\((\d+),\s*\d+\)", a_ampl)
                                zeros_match = re.search(r"zeros\((\d+),\s*\d+\)", a_ampl)
                                if ones_match:
                                    num_modules = int(ones_match.group(1))
                                elif zeros_match:
                                    num_modules = int(zeros_match.group(1))
                                else:
                                    # Fallback: use the last number in the expression
                                    a_ampl_numbers = re.findall(r"\d+", a_ampl)
                                    if a_ampl_numbers:
                                        num_modules = int(a_ampl_numbers[-1])
                    elif "a_phase" in matlab_scenario["antenna"]:
                        a_phase = matlab_scenario["antenna"]["a_phase"]
                        if isinstance(a_phase, (list, np.ndarray)):
                            num_modules = len(a_phase)
                        elif isinstance(a_phase, str):
                            # Try to extract the number from MATLAB expressions
                            a_phase_numbers = re.findall(r"\d+", a_phase)
                            if a_phase_numbers:
                                if ":" in a_phase:
                                    # Handle colon operator (e.g., 0:7 -> 8 elements)
                                    colon_parts = a_phase.split(":")
                                    if len(colon_parts) >= 2:
                                        start = int(colon_parts[0].strip())
                                        end = int(colon_parts[1].strip().rstrip(")' "))
                                        num_modules = end - start + 1
                                else:
                                    # Use the first number in the expression
                                    num_modules = int(a_phase_numbers[0])

                if num_modules is None:
                    raise ValueError(
                        f"Cannot determine the number of modules for architecture '{architecture_str}'. "
                        "Please ensure the number of modules exists or provide 'a_ampl' or 'a_phase' as arrays."
                    )

        # antenna.a_ampl -> antenna.excitation.power
        if "antenna" in matlab_scenario and "a_ampl" in matlab_scenario["antenna"]:
            a_ampl = matlab_scenario["antenna"]["a_ampl"]
            if isinstance(a_ampl, (list, np.ndarray)):
                # Handle nested lists/arrays from .mat files
                # Check if this is a 2D array with a single row (common MATLAB format)
                if isinstance(a_ampl, np.ndarray) and a_ampl.ndim == 2 and a_ampl.shape[0] == 1:
                    a_ampl = a_ampl[0]  # Extract the first row
                elif isinstance(a_ampl, list) and a_ampl and isinstance(a_ampl[0], list):
                    a_ampl = a_ampl[0]  # Extract the inner list
                # Convert to list of floats (assume same unit as MATLAB)
                excitation["power"] = [float(x) for x in a_ampl]
            else:
                # If not a straightforward array, raise an error if num_modules is not available
                if num_modules is None:
                    raise ValueError(
                        "Cannot determine the number of modules and 'a_ampl' is not a valid array. "
                        "Please provide 'a_ampl' as a list or array."
                    )
                # Default: all modules have power of 1.0
                excitation["power"] = [1.0] * num_modules

        # antenna.a_phase -> antenna.excitation.phase (in degree)
        if "antenna" in matlab_scenario and "a_phase" in matlab_scenario["antenna"]:
            a_phase = matlab_scenario["antenna"]["a_phase"]
            if isinstance(a_phase, (list, np.ndarray)):
                # Handle nested lists/arrays from .mat files
                # Check if this is a 2D array with a single row (common MATLAB format)
                if isinstance(a_phase, np.ndarray) and a_phase.ndim == 2 and a_phase.shape[0] == 1:
                    a_phase = a_phase[0]  # Extract the first row
                elif isinstance(a_phase, list) and a_phase and isinstance(a_phase[0], list):
                    a_phase = a_phase[0]  # Extract the inner list
                # Convert radians to degrees and normalize modulo 360
                excitation["phase"] = [float(np.rad2deg(x)) % 360 for x in a_phase]
            else:
                # If not a straightforward array, raise an error if num_modules is not available
                if num_modules is None:
                    raise ValueError(
                        "Cannot determine the number of modules and 'a_phase' is not a valid array. "
                        "Please provide 'a_phase' as a list or array."
                    )
                # Default: phases are 0, 90, 180, 270, ... degrees (cyclic)
                default_phases = []
                for i in range(num_modules):
                    default_phases.append((i % 4) * 90.0)
                excitation["phase"] = default_phases

        # options.TSport -> antenna.excitation.port
        if "options" in matlab_scenario and "TSport" in matlab_scenario["options"]:
            TSport = matlab_scenario["options"]["TSport"]
            # Handle TSport as ASCII codes (from .mat files) or as string
            if isinstance(TSport, (list, np.ndarray)):
                # Convert ASCII codes to string
                if isinstance(TSport, list) and TSport and isinstance(TSport[0], list):
                    # Handle nested lists from .mat files
                    TSport = [item[0] if isinstance(item, list) else item for item in TSport]
                elif isinstance(TSport, np.ndarray) and TSport.ndim == 2 and TSport.shape[1] == 1:
                    # Handle 2D numpy array of character codes
                    TSport = TSport.flatten()
                excitation["port"] = "".join(chr(int(code)) for code in TSport)
            else:
                excitation["port"] = str(TSport)

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
        version_plasma_1d = None
        if version_code == "1D":
            version_plasma_1d = matlab_scenario.get("version_plasma_1D") or get_nested(
                matlab_scenario, "plasma", "version"
            )
        elif version_code == "2D":
            version_plasma_1d = matlab_scenario.get("version_plasma_2D") or get_nested(
                matlab_scenario, "plasma", "version"
            )

        if version_plasma_1d is None:
            # Fallback to plasma.version
            version_plasma_1d = get_nested(matlab_scenario, "plasma", "version")

        if isinstance(version_plasma_1d, str):
            # This is a variable reference, get the actual value
            if version_plasma_1d == "version_plasma_1D":
                version_plasma_1d = matlab_scenario.get("version_plasma_1D")
            elif version_plasma_1d == "version_plasma_2D":
                version_plasma_1d = matlab_scenario.get("version_plasma_2D")

        # Convert to int for comparison (MATLAB may store version as float)
        if isinstance(version_plasma_1d, (int, float)):
            version_plasma_1d_int = int(version_plasma_1d)
            if version_plasma_1d_int == 3:
                spectral_1d["profile"] = "linear"
            elif version_plasma_1d_int == 6:
                spectral_1d["profile"] = "bilinear"

        # Nme -> plasma.spectral_1D.nb_evanescent_modes
        if "Nme" in matlab_scenario:
            spectral_1d["nb_evanescent_modes"] = int(matlab_scenario["Nme"])
        elif "options" in matlab_scenario and "modes" in matlab_scenario["options"]:
            # Try to extract nb_evanescent_modes from modes field in options
            modes = matlab_scenario["options"]["modes"]
            if isinstance(modes, (list, np.ndarray)) and len(modes) >= 2:
                # The second element in modes might be the number of evanescent modes
                if isinstance(modes[1], (list, np.ndarray)) and len(modes[1]) > 0:
                    spectral_1d["nb_evanescent_modes"] = int(modes[1][0])
                elif len(modes) >= 2:
                    spectral_1d["nb_evanescent_modes"] = int(modes[1])

        # Create bilinear section
        bilinear = {}

        # options.bool_lignes_identiques -> plasma.spectral_1D.bilinear.identical_profiles
        if "options" in matlab_scenario and "bool_lignes_identiques" in matlab_scenario["options"]:
            bool_lignes_identiques = matlab_scenario["options"]["bool_lignes_identiques"]
            # Handle MATLAB boolean strings properly ("true"/"false")
            if isinstance(bool_lignes_identiques, str):
                bilinear["identical_profiles"] = bool_lignes_identiques.lower() == "true"
            else:
                bilinear["identical_profiles"] = bool(bool_lignes_identiques)

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
                # Handle nested lists from .mat files (e.g., [[0.002], [0.02]])
                if isinstance(lambda_n_value, list) and lambda_n_value and isinstance(lambda_n_value[0], list):
                    # Extract all values from nested lists
                    for item in lambda_n_value:
                        if isinstance(item, list) and len(item) > 0:
                            lambda_n.append(float(item[0]))
                        elif isinstance(item, (int, float, np.number)):
                            lambda_n.append(float(item))
                elif (
                    isinstance(lambda_n_value, np.ndarray) and lambda_n_value.ndim == 2 and lambda_n_value.shape[1] == 1
                ):
                    # Handle 2D numpy array with single column
                    for item in lambda_n_value:
                        if isinstance(item, (int, float, np.number)):
                            lambda_n.append(float(item))
                        elif hasattr(item, "__iter__") and not isinstance(item, str):
                            # It's an iterable (like a 1-element array)
                            lambda_n.append(float(item[0]))
                else:
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

        # options.B0 -> plasma.spectral_1D.bilinear.B0
        if "options" in matlab_scenario and "B0" in matlab_scenario["options"]:
            bilinear["B0"] = float(matlab_scenario["options"]["B0"])

        # plasma.d_vide -> plasma.spectral_1D.bilinear.vacuum_layer_length
        if "plasma" in matlab_scenario and "d_vide" in matlab_scenario["plasma"]:
            bilinear["vacuum_layer_length"] = float(matlab_scenario["plasma"]["d_vide"])

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
            bool_debug = matlab_scenario["options"]["bool_debug"]
            # Handle MATLAB boolean strings properly ("true"/"false")
            if isinstance(bool_debug, str):
                options["debug"] = bool_debug.lower() == "true"
            else:
                options["debug"] = bool(bool_debug)

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
        if self.debug:
            print("[DEBUG] Plasma configuration:")
            pp(plasma)

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

                if self.debug:
                    print("[DEBUG] Sum S_plasma:", np.sum(S_plasma))
                    print("[DEBUG] Sum rac_Zhe", np.sum(rac_Zhe))

                # Store results in the scenario
                self.results["S_plasma"] = S_plasma
                self.results["rac_Zhe"] = rac_Zhe

                # Compute antenna response (reflection coefficients)
                self._compute_antenna_response()

                if self.debug:
                    print("[DEBUG] RC", self.results["RC"])

                return
            else:
                raise ValueError(
                    f"Unsupported plasma profile '{profile}' for spectral_1D solver. Only 'bilinear' is supported."
                )
        else:
            raise ValueError(f"Unsupported solver '{solver}'. Only 'spectral_1D' is currently supported.")

    def _load_sparameters_from_files(
        self,
        filenames: list,
        phases_deembedded: list,
        nb_access_ports: int,
        nb_plasma_ports: int,
        nb_g_total_ligne: int,
        nb_modes_total: int,
        S_plasma: np.ndarray,
        antenna_data: dict,
    ) -> tuple:
        """
        Load S-parameter files and assemble global S_ant matrices.

        This implements the logic from MATLAB's S_antenne.m to:
        1. Load each module's S-parameter file
        2. Extract the S matrix and apply phase deembedding
        3. Assemble the global S_ant matrices
        """
        from pathlib import Path

        from aloha.sparameters import load_sparameter_file

        # Initialize S_ant matrices with transposed convention (to match .mat files)
        S_ant_11 = np.zeros((nb_access_ports, nb_access_ports), dtype=complex)
        S_ant_12 = np.zeros((nb_plasma_ports, nb_access_ports), dtype=complex)
        S_ant_21 = np.zeros((nb_access_ports, nb_plasma_ports), dtype=complex)
        S_ant_22 = np.zeros((nb_plasma_ports, nb_plasma_ports), dtype=complex)

        # Initialize S_ant_22 with passive waveguide values
        # Passive waveguides have short circuits that reflect waves
        # Based on MATLAB S_antenne.m: S_ant_22(1,pass_tot) = -exp(+i*4*pi*lcc)
        # where lcc is the short circuit depth and pass_tot are the passive waveguide indices

        # Get antenna module parameters
        module_data = antenna_data.get("module", {})
        nb_wg_phi = module_data.get("nb_wg_phi", 1)
        nb_wg_theta = module_data.get("nb_wg_theta", 1)
        mask = module_data.get("mask", [1])
        nb_pwg_btw_mod_phi = module_data.get("nb_pwg_btw_mod_phi", 0)
        nb_pwg_edge = module_data.get("nb_pwg_edge", 0)
        pwg_depth = module_data.get("pwg_depth", [0.25])

        # Calculate lcc (short circuit depth in guided wavelengths)
        # For now, use the first pwg_depth value
        lcc_default = pwg_depth[0] if isinstance(pwg_depth, list) and len(pwg_depth) > 0 else 0.25

        # Get the number of modules
        layout_data = antenna_data.get("layout", {})
        nb_modules_tor = layout_data.get("nb_mod_phi", 1)
        nb_modules_pol = layout_data.get("nb_mod_theta", 1)

        # Calculate the number of active waveguides per module
        # mask is a list of 0s and 1s, where 1 = active, 0 = passive
        nb_active_wg_phi = sum(mask) if isinstance(mask, list) else 1

        # Identify passive waveguide indices based on antenna geometry
        # This includes both passive waveguides within modules (from mask) and between modules
        passive_wg_indices = []

        # Calculate waveguide positions
        # Total waveguides per module in toroidal direction: nb_wg_phi
        # Passive waveguides between modules: nb_pwg_btw_mod_phi
        # Passive waveguides on edges: nb_pwg_edge

        # Edge passive waveguides at the beginning
        for i in range(nb_pwg_edge):
            passive_wg_indices.append(i)

        # Active and passive waveguides for each module, and passive waveguides between modules
        current_pos = nb_pwg_edge
        for mod in range(nb_modules_tor):
            # Within each module, identify passive waveguides based on mask
            for wg_offset in range(nb_wg_phi):
                if wg_offset < len(mask) and mask[wg_offset] == 0:
                    # This waveguide is passive
                    passive_wg_indices.append(current_pos + wg_offset)

            # Move past all waveguides in this module (active + passive)
            current_pos += nb_wg_phi

            # Passive waveguides between modules (if any and not at the last module)
            if nb_pwg_btw_mod_phi > 0 and mod < nb_modules_tor - 1:
                for i in range(nb_pwg_btw_mod_phi):
                    passive_wg_indices.append(current_pos + i)
                current_pos += nb_pwg_btw_mod_phi

        # Edge passive waveguides at the end
        for i in range(nb_pwg_edge):
            passive_wg_indices.append(current_pos + i)

        # Set diagonal values for passive waveguides
        # Passive waveguides exist in all poloidal rows
        # Note: MATLAB only sets S_ant_22 for mode 0 of each passive waveguide
        for pol_row in range(nb_wg_theta):
            for wg_idx in passive_wg_indices:
                if wg_idx < nb_g_total_ligne:  # Make sure it's within bounds
                    # Only set mode 0 (matching MATLAB behavior)
                    mode = 0
                    # Calculate plasma port index accounting for poloidal row
                    wg_index = wg_idx + pol_row * nb_g_total_ligne
                    plasma_port = (wg_index + 1) * nb_modes_total + mode - (nb_modes_total - 1) - 1
                    # S_ant_22 is diagonal in MATLAB, so only set diagonal elements
                    # In our transposed convention: S_ant_22[plasma_port, plasma_port]
                    S_ant_22[plasma_port, plasma_port] = -np.exp(1j * 4 * np.pi * lcc_default)

        # Precompute the active waveguide indices for each module
        # This is similar to MATLAB's modules_act
        # For each module, we need to know which plasma ports correspond to its active waveguides
        modules_act = []  # List of lists, where modules_act[mod] = list of plasma port indices for active waveguides
        for mod in range(nb_modules_tor * nb_modules_pol):
            # Calculate the poloidal row and toroidal position for this module
            pol_row = mod // nb_modules_tor
            tor_pos = mod % nb_modules_tor

            # Calculate the starting waveguide index for this module in the toroidal line
            # This includes edge passive waveguides and waveguides from previous modules
            wg_start = nb_pwg_edge + tor_pos * (nb_wg_phi + nb_pwg_btw_mod_phi)

            # For each waveguide in the module, check if it's active (based on mask)
            # The S-parameter files may include waveguides from multiple poloidal rows
            # Note: Following MATLAB's convention, we only use mode 0 for S_ant matrices
            active_plasma_ports = []
            for pol_offset in range(nb_wg_theta):
                for wg_offset in range(nb_wg_phi):
                    if wg_offset < len(mask) and mask[wg_offset] == 1:
                        # This waveguide is active
                        # Calculate the waveguide index accounting for poloidal rows
                        wg_index = wg_start + wg_offset + pol_offset * nb_g_total_ligne
                        # Calculate the plasma port index for mode 0 of this waveguide
                        # This matches MATLAB's convention: modules_act = modules_act * (Nme+Nmh) - (Nme+Nmh-1)
                        # In 0-based indexing: plasma_port = wg_index * nb_modes_total + 0
                        plasma_port = wg_index * nb_modes_total
                        active_plasma_ports.append(plasma_port)

            modules_act.append(active_plasma_ports)

        # For each module, load its S-parameter file
        for ind in range(len(filenames)):
            filename = filenames[ind]
            phase_deembedded = phases_deembedded[ind] if ind < len(phases_deembedded) else 0.0

            # Try to find the S-parameter file
            sparam_file = None
            search_paths = [
                Path(filename),
                Path(filename + ".m"),
                Path(filename + ".mat"),
                # For WEST_LH1
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C3" / filename,
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C3" / (filename + ".m"),
                Path(__file__).parent.parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C3" / filename,
                Path(__file__).parent.parent.parent
                / "aloha_matlab"
                / "S_HFSS"
                / "matrices_HFSS_C3"
                / (filename + ".m"),
                # For WEST_LH2
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C4" / filename,
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C4" / (filename + ".m"),
                Path(__file__).parent.parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_C4" / filename,
                Path(__file__).parent.parent.parent
                / "aloha_matlab"
                / "S_HFSS"
                / "matrices_HFSS_C4"
                / (filename + ".m"),
                # For 8_active_waveguides
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_elem" / filename,
                Path(__file__).parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_elem" / (filename + ".m"),
                Path(__file__).parent.parent.parent / "aloha_matlab" / "S_HFSS" / "matrices_HFSS_elem" / filename,
                Path(__file__).parent.parent.parent
                / "aloha_matlab"
                / "S_HFSS"
                / "matrices_HFSS_elem"
                / (filename + ".m"),
            ]

            for path in search_paths:
                if path.exists():
                    sparam_file = path
                    break

            if sparam_file is None:
                raise FileNotFoundError(f"S-parameter file not found for module {ind}. Filename: {filename}")

            # Load the S-parameter file
            try:
                sparam_data = load_sparameter_file(sparam_file)
                S_module = sparam_data.get("S", None)

                if S_module is None:
                    raise ValueError(f"S-parameter matrix not found in file: {sparam_file}")

                # Extract S_module elements
                if len(S_module.shape) == 1:
                    n = int(np.sqrt(len(S_module)))
                    S_module = S_module.reshape((n, n))

                # S_module is (n_ports, n_ports) where n_ports = 1 + nb_wg_per_module * nb_modes
                # For WEST_LH1: n_ports = 1 + 6 * 3 = 19 (6 active waveguides * 3 poloidal rows * 1 mode)
                n_ports = S_module.shape[0]
                # Calculate the number of modes in the S-parameter file
                # The S-parameter files may include waveguides from multiple poloidal rows
                # Total waveguides in S-parameter file = nb_active_wg_phi * nb_wg_theta
                nb_wg_in_sparam = nb_active_wg_phi * max(nb_wg_theta, 1)
                nb_modes_sparam = (n_ports - 1) // nb_wg_in_sparam
                nb_wg_per_module = (n_ports - 1) // nb_modes_sparam

                S_module_11 = S_module[0, 0]

                # Extract only mode 0 elements from S_module_12, S_module_21, and S_module_22
                # S_module has shape (n_ports, n_ports) where n_ports = 1 + nb_wg_per_module * nb_modes_sparam
                # Waveguide ports are at indices 1, 2, ..., nb_wg_per_module * nb_modes_sparam
                # Mode 0 ports are at indices 1, 1+nb_modes_sparam, 1+2*nb_modes_sparam, ...
                # In Python (0-based): indices 1::nb_modes_sparam
                S_module_12 = S_module[0, 1::nb_modes_sparam]
                S_module_21 = S_module[1::nb_modes_sparam, 0]

                # Extract mode 0 submatrix from S_module_22
                # This selects rows and columns at mode 0 waveguide ports
                S_module_22 = S_module[1::nb_modes_sparam, 1::nb_modes_sparam]

                # Note: We only use mode 0 for S_ant matrices, so no expansion is needed
                # even if nb_modes_sparam != nb_modes_total.
                # The S_module matrices already contain the mode 0 data for all waveguides.

                # Update nb_wg_per_module to match the extracted size
                nb_wg_per_module = S_module_22.shape[0]

                # Note: Phase deembedding is only applied when bool_mesure = true in MATLAB
                # For this scenario (WEST_LH1), bool_mesure = false, so we don't apply it
                # phase_deembedded_rad = np.deg2rad(phase_deembedded)
                # if phase_deembedded != 0:
                #     S_module_11 = S_module_11 * np.exp(1j * 2 * phase_deembedded_rad)
                #     S_module_12 = S_module_12 * np.exp(1j * phase_deembedded_rad)
                #     S_module_21 = S_module_21 * np.exp(1j * phase_deembedded_rad)

                # Place S_module values into global S_ant matrices
                # S_module_12 has shape (nb_wg_per_module * nb_modes,)
                # S_module_21 has shape (nb_wg_per_module * nb_modes, 1)

                # Check if we have modules_act (active waveguide indices per module)
                # If nb_wg_per_module matches the total active waveguides per module (nb_active_wg_phi * nb_wg_theta),
                # then S-parameter files include only active waveguides
                # In this case, we use modules_act to place the S-parameter data
                if nb_wg_per_module == nb_active_wg_phi * nb_wg_theta and len(modules_act) > ind:
                    # S-parameter files include only active waveguides
                    # Use modules_act to get the correct plasma port indices
                    active_plasma_ports = modules_act[ind]

                    # Place S_module_11 on the diagonal
                    S_ant_11[ind, ind] = S_module_11

                    # Place S_module_12: connects access port to active waveguide ports
                    # Note: Python uses transposed convention to match .mat files
                    # MATLAB code: S_ant_12(ind, modules_act(ind,:)) = S_module_12
                    # Python convention: S_ant_12 has shape (nb_plasma_ports, nb_access_ports)
                    # So we need to transpose: S_ant_12[plasma_port, ind] = s_val
                    for i, s_val in enumerate(S_module_12):
                        if i < len(active_plasma_ports):
                            plasma_port = active_plasma_ports[i]
                            S_ant_12[plasma_port, ind] = s_val

                    # Place S_module_21: connects active waveguide ports to access port
                    # Note: Python uses transposed convention to match .mat files
                    # MATLAB code: S_ant_21(modules_act(ind,:), ind) = S_module_21
                    # Python convention: S_ant_21 has shape (nb_access_ports, nb_plasma_ports)
                    # So we need to transpose: S_ant_21[ind, plasma_port] = s_val
                    for i, s_val in enumerate(S_module_21):
                        if i < len(active_plasma_ports):
                            plasma_port = active_plasma_ports[i]
                            S_ant_21[ind, plasma_port] = s_val

                    # Place S_module_22: scattering between active waveguide ports
                    for i in range(len(S_module_21)):
                        for j in range(len(S_module_21)):
                            if i < len(active_plasma_ports) and j < len(active_plasma_ports):
                                plasma_port_i = active_plasma_ports[i]
                                plasma_port_j = active_plasma_ports[j]
                                S_ant_22[plasma_port_i, plasma_port_j] = S_module_22[i, j]
                else:
                    # S-parameter files include all waveguides (active + passive)
                    # Use the original logic with contiguous waveguide indices

                    # Calculate waveguide start index for this module
                    waveguide_start = nb_pwg_edge + ind * (nb_wg_phi + nb_pwg_btw_mod_phi)

                    # Place S_module_11 on the diagonal
                    S_ant_11[ind, ind] = S_module_11

                    # Place S_module_12: connects access port to waveguide ports
                    # Note: Python uses transposed convention to match .mat files
                    # MATLAB code: S_ant_12(ind, modules_act(ind,:)) = S_module_12
                    # Python convention: S_ant_12 has shape (nb_plasma_ports, nb_access_ports)
                    # So we need to transpose: S_ant_12[plasma_port, ind] = s_val
                    for i, s_val in enumerate(S_module_12):
                        mode_index = i // nb_wg_per_module
                        waveguide_offset = i % nb_wg_per_module
                        waveguide_index = waveguide_start + waveguide_offset
                        plasma_port = (waveguide_index + 1) * nb_modes_total + mode_index - (nb_modes_total - 1) - 1
                        S_ant_12[plasma_port, ind] = s_val

                    # Place S_module_21: connects waveguide ports to access port
                    # Note: Python uses transposed convention to match .mat files
                    # MATLAB code: S_ant_21(modules_act(ind,:), ind) = S_module_21
                    # Python convention: S_ant_21 has shape (nb_access_ports, nb_plasma_ports)
                    # So we need to transpose: S_ant_21[ind, plasma_port] = s_val
                    for i, s_val in enumerate(S_module_21):
                        mode_index = i // nb_wg_per_module
                        waveguide_offset = i % nb_wg_per_module
                        waveguide_index = waveguide_start + waveguide_offset
                        plasma_port = (waveguide_index + 1) * nb_modes_total + mode_index - (nb_modes_total - 1) - 1
                        S_ant_21[ind, plasma_port] = s_val

                    # Place S_module_22: scattering between waveguide ports
                    for i in range(len(S_module_21)):
                        for j in range(len(S_module_21)):
                            mode_index_i = i // nb_wg_per_module
                            waveguide_offset_i = i % nb_wg_per_module
                            waveguide_index_i = waveguide_start + waveguide_offset_i
                            plasma_port_i = (
                                (waveguide_index_i + 1) * nb_modes_total + mode_index_i - (nb_modes_total - 1) - 1
                            )

                            mode_index_j = j // nb_wg_per_module
                            waveguide_offset_j = j % nb_wg_per_module
                            waveguide_index_j = waveguide_start + waveguide_offset_j
                            plasma_port_j = (
                                (waveguide_index_j + 1) * nb_modes_total + mode_index_j - (nb_modes_total - 1) - 1
                            )
                            S_ant_22[plasma_port_i, plasma_port_j] = S_module_22[i, j]
            except Exception as e:
                raise RuntimeError(f"Error loading S-parameter file {sparam_file}: {str(e)}")

        return S_ant_11, S_ant_12, S_ant_21, S_ant_22

    def _compute_antenna_response(self) -> None:
        """
        Compute the antenna response (reflection coefficients) from S_plasma and rac_Zhe.

        This implements the logic from MATLAB's reponse_antenne.m and aloha_compute_RC.m
        to connect the plasma S-parameters to the antenna S-parameters and compute
        the reflection coefficients.

        Each scenario will compute its own S_ant matrices:
        - For .toml and .m files: load S-parameter files from antenna definition
        - For .mat files: use S_ant matrices already in the file

        Note: S-parameter files must be available for all modules. If any S-parameter
        file is missing or cannot be loaded, an error will be raised (no simplified
        model fallback is used).

        """
        # Check if S_ant matrices are already computed (e.g., loaded from .mat file)
        if all(key in self.results for key in ["S_ant_11", "S_ant_12", "S_ant_21", "S_ant_22"]):
            # S_ant matrices are already available, use them directly
            S_ant_11 = self.results["S_ant_11"]
            S_ant_12 = self.results["S_ant_12"]
            S_ant_21 = self.results["S_ant_21"]
            S_ant_22 = self.results["S_ant_22"]
        else:
            # Need to compute S_ant matrices
            # Initialize S_ant matrices
            S_ant_11 = None
            S_ant_12 = None
            S_ant_21 = None
            S_ant_22 = None
        # Get S_plasma and rac_Zhe from results
        S_plasma = self.results.get("S_plasma")
        if S_plasma is None:
            raise ValueError("S_plasma must be computed before calling _compute_antenna_response")

        # Get antenna excitation parameters from scenario
        scenario_antenna = self.scenario.get("antenna", {})
        excitation = scenario_antenna.get("excitation", {})

        # Get magnitudes and phases
        a_ampl = np.array(excitation.get("power", []), dtype=float)
        a_phase = np.array(excitation.get("phase", []), dtype=float)

        # Convert phases from degrees to radians if needed
        if len(a_phase) > 0 and np.max(np.abs(a_phase)) > 10:
            a_phase = np.deg2rad(a_phase)

        # Incident wave vector on antenna
        a_acces = a_ampl * np.exp(1j * a_phase)

        # Load antenna file to get layout and module parameters
        antenna_file = scenario_antenna.get("file")
        antenna_data = {}
        if antenna_file:
            try:
                from pathlib import Path

                from aloha.antenna import Antenna

                # Try to find the antenna file
                antenna_paths = [
                    Path(antenna_file),
                    Path(__file__).parent.parent / "antennas" / antenna_file,
                    Path(__file__).parent.parent.parent / "antennas" / antenna_file,
                ]

                # Also try with .m extension
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
            except (FileNotFoundError, ImportError) as e:
                # If we can't load the antenna file, use defaults
                raise (FileNotFoundError, "can't load antenna file")

        # Get layout parameters
        layout = antenna_data.get("layout", {})
        nb_modules_tor = layout.get("nb_mod_phi", 1)
        nb_modules_pol = layout.get("nb_mod_theta", 1)
        total_modules = nb_modules_tor * nb_modules_pol

        # Get module parameters
        module = antenna_data.get("module", {})
        nb_wg_phi = module.get("nb_wg_phi", 1)
        nb_wg_theta = module.get("nb_wg_theta", 1)
        mask = module.get("mask", [1])
        nb_pwg_edge = module.get("nb_pwg_edge", 0)
        nb_pwg_btw_mod_phi = module.get("nb_pwg_btw_mod_phi", 0)

        # Calculate total number of waveguides per poloidal row
        nb_g_total_ligne = nb_wg_phi * nb_modules_tor + 2 * nb_pwg_edge + nb_pwg_btw_mod_phi * (nb_modules_tor - 1)

        # Get the number of modes from S_plasma shape
        # Total waveguides = nb_g_total_ligne * nb_wg_theta (accounting for poloidal rows)
        total_waveguides = nb_g_total_ligne * nb_wg_theta
        nb_modes_total = S_plasma.shape[0] // total_waveguides if total_waveguides > 0 else 1

        # Number of access ports = number of modules
        nb_access_ports = total_modules

        # Number of plasma ports = total_waveguides * nb_modes_total
        nb_plasma_ports = total_waveguides * nb_modes_total

        # Initialize S_ant matrices only if not already loaded
        if S_ant_11 is None:
            # Try to load S-parameter files from antenna definition
            sparameters = antenna_data.get("sparameters", {})
            filenames = sparameters.get("filenames", [])
            phases_deembedded = sparameters.get("phases_deembedded", [])

            if filenames and len(filenames) == total_modules:
                # Load S-parameters from files
                S_ant_11, S_ant_12, S_ant_21, S_ant_22 = self._load_sparameters_from_files(
                    filenames,
                    phases_deembedded,
                    nb_access_ports,
                    nb_plasma_ports,
                    nb_g_total_ligne,
                    nb_modes_total,
                    S_plasma,
                    antenna_data,
                )
            else:
                # S-parameter files are required
                raise ValueError(
                    f"S-parameter files must be defined for all {total_modules} modules. "
                    f"Found {len(filenames) if filenames else 0} filenames in antenna definition."
                )

        # Store antenna S-parameters in results
        self.results["S_ant_11"] = S_ant_11
        self.results["S_ant_12"] = S_ant_12
        self.results["S_ant_21"] = S_ant_21
        self.results["S_ant_22"] = S_ant_22

        # Store access wave vectors
        self.results["a_acces"] = a_acces

        # Compute incident and reflected wave vectors on/from plasma
        # From reponse_antenne.m:
        # a_plasma = inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21*a_acces
        # b_plasma = S_plasma*a_plasma

        # Note: The S_ant matrices in .mat files use a transposed convention
        # relative to MATLAB. So we need to use .T to match the MATLAB formulas.
        # In MATLAB: S_ant_21 is (nb_plasma_ports, nb_access_ports)
        # In our code: S_ant_21 is (nb_access_ports, nb_plasma_ports)
        # So S_ant_21.T gives us (nb_plasma_ports, nb_access_ports) which matches MATLAB

        identity = np.eye(S_plasma.shape[0], dtype=complex)

        # Reshape a_acces to be a column vector
        a_acces_col = a_acces.reshape(-1, 1)

        # Check if S_ant_22 is zero (no passive waveguides)
        if np.allclose(S_ant_22, 0):
            # Simplified case: a_plasma = S_ant_21.T @ a_acces
            # This matches the MATLAB formula when S_ant_21 is properly transposed
            a_plasma = S_ant_21.T @ a_acces_col
        else:
            # General case
            # a_plasma = inv(eye(length(S_plasma)) - S_ant_22*S_plasma) * S_ant_21 * a_acces
            matrix_to_invert = identity - S_ant_22 @ S_plasma
            try:
                inv_matrix = np.linalg.inv(matrix_to_invert)
                # a_plasma = inv(...) @ S_ant_21.T @ a_acces
                a_plasma = inv_matrix @ S_ant_21.T @ a_acces_col
            except np.linalg.LinAlgError:
                # If matrix is singular, use pseudo-inverse
                inv_matrix = np.linalg.pinv(matrix_to_invert)
                a_plasma = inv_matrix @ S_ant_21.T @ a_acces_col

        # Compute b_plasma
        # b_plasma = S_plasma * a_plasma
        b_plasma = S_plasma @ a_plasma

        # Store plasma wave vectors
        self.results["a_plasma"] = a_plasma.flatten()
        self.results["b_plasma"] = b_plasma.flatten()

        # Reflection coefficient at the mouth of the antenna
        # RC_mouth = 100*abs(b_plasma./a_plasma).^2
        with np.errstate(divide="ignore", invalid="ignore"):
            rc_mouth = 100 * np.abs(b_plasma / a_plasma) ** 2
            rc_mouth[~np.isfinite(rc_mouth)] = 0.0
        self.results["RC_mouth"] = rc_mouth.flatten()

        # Compute the plasma-coupled antenna scattering matrix
        # From MATLAB: S_acces = S_ant_11 + S_ant_12*S_plasma*inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21
        # But our S_ant matrices use transposed convention, so we need to use .T
        if np.allclose(S_ant_22, 0):
            # Simplified case: S_acces = S_ant_11 + S_ant_12.T @ S_plasma @ S_ant_21.T
            S_acces = S_ant_11 + S_ant_12.T @ S_plasma @ S_ant_21.T
        else:
            # General case
            inv_matrix = np.linalg.inv(identity - S_ant_22 @ S_plasma)
            S_acces = S_ant_11 + S_ant_12.T @ S_plasma @ inv_matrix @ S_ant_21.T

        # Store S_acces
        self.results["S_acces"] = S_acces

        # Compute reflected wave vector from antenna
        # b_acces = S_acces * a_acces
        b_acces = S_acces @ a_acces_col

        # Store b_acces
        self.results["b_acces"] = b_acces.flatten()

        # Power reflection coefficient at input of a module
        # CoeffRefPuiss = 100*abs(b_acces./a_acces).^2
        with np.errstate(divide="ignore", invalid="ignore"):
            coeff_ref_puiss = 100 * np.abs(b_acces / a_acces_col) ** 2
            coeff_ref_puiss[~np.isfinite(coeff_ref_puiss)] = 0.0

        # Store reflection coefficients
        self.results["CoeffRefPuiss"] = coeff_ref_puiss.flatten()
        self.results["RC"] = coeff_ref_puiss.flatten()


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
        except FileNotFoundError:
            raise (FileNotFoundError, "Antenna S-parameters files not found")

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

    # Fix: nb_g_pol should be the number of poloidal waveguide rows (nb_wg_theta),
    # not the number of poloidal modules (nb_mod_theta).
    # This matches the MATLAB behavior in aloha_utils_ITM2oldAntenna.m line 6:
    # "nb_g_pol = aloha_scenario_get(scenario, 'nwm_theta'); % 11/10/2013 - was nma_theta, but does not work..."
    nb_g_pol = nb_wg_theta

    # Calculate total number of waveguides per poloidal row
    # Using the MATLAB formula: nb_g_total_ligne = nb_wg_phi * nb_mod_phi + 2 * nb_pwg_edge
    # + nb_pwg_btw_mod_phi * (nb_mod_phi - 1)
    # This matches the waveguide logic in aloha_utils_getAntennaCoordinates.m
    nb_g_total_ligne = nb_wg_phi * nb_mod_phi + 2 * nb_pwg_edge + nb_pwg_btw_mod_phi * (nb_mod_phi - 1)

    # Waveguide dimensions from antenna module parameters
    wg_size_theta = module.get("wg_size_theta", 70e-3)  # Height of waveguides in poloidal direction [m]
    awg_size_phi = module.get("awg_size_phi", 10e-3)  # Width of active waveguides [m]
    pwg_size_phi = module.get("pwg_size_phi", 6.5e-3)  # Width of internal passive waveguides [m]
    pwg_size_edge_phi = module.get("pwg_size_edge_phi", 6.5e-3)  # Width of edge passive waveguides [m]
    e_phi = module.get("e_phi", 2e-3)  # Spacing between active waveguides [m]
    e_phi_pwg = module.get("e_phi_pwg", 3e-3)  # Spacing between passive waveguides [m]

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

    # Waveguide parameters using MATLAB waveguide logic from aloha_utils_getAntennaCoordinatesFromCPO.m
    # a = waveguide height in poloidal direction (constant for all waveguides in a line)
    a = wg_size_theta

    # Compute b array (waveguide widths in toroidal direction) to match MATLAB logic
    # From aloha_utils_getAntennaCoordinatesFromCPO.m lines 31-41:
    # b_module = wg.mask.*wg.bwa + not(wg.mask).*wg.biwp;
    # b_edge = repmat(wg.bewp, 1, wg.npwe_phi);
    # b_inter = repmat(wg.biwp, 1, wg.npwbm_phi);
    # b = [b_edge, kron(ones(1,mod.nma_phi-1),[b_module, b_inter]),b_module, b_edge];

    # b_module: waveguide widths within a module (active or internal passive)
    b_module = []
    for m in mask:
        if m == 1:
            b_module.append(awg_size_phi)  # active waveguide
        else:
            b_module.append(pwg_size_phi)  # internal passive waveguide

    # b_edge: passive waveguide widths on each edge
    b_edge = [pwg_size_edge_phi] * nb_pwg_edge

    # b_inter: passive waveguide widths between modules
    b_inter = [pwg_size_phi] * nb_pwg_btw_mod_phi

    # Construct b array following MATLAB logic
    # [b_edge, kron(ones(1, nb_mod_phi-1), [b_module, b_inter]), b_module, b_edge]
    b_parts = []
    b_parts.extend(b_edge)  # leading edge
    for _ in range(nb_mod_phi - 1):
        b_parts.extend(b_module)
        b_parts.extend(b_inter)
    b_parts.extend(b_module)  # last module
    b_parts.extend(b_edge)  # trailing edge

    b = b_parts

    # Compute e array (septum widths between waveguides)
    # From MATLAB: e = wg.e_phi (this is an array)
    # The spacing depends on the waveguide types:
    # - e_phi_pwg (3e-3) when at least one adjacent waveguide is passive
    # - e_phi (2e-3) when both adjacent waveguides are active
    # This matches the pattern seen in MATLAB antenna structures
    e = []
    for i in range(len(b) - 1):
        # Check if current or next waveguide is passive (width != awg_size_phi)
        if b[i] != awg_size_phi or b[i + 1] != awg_size_phi:
            e.append(e_phi_pwg)  # Use passive spacing
        else:
            e.append(e_phi)  # Use active spacing

    # Compute z array (waveguide positions in toroidal direction)
    # From MATLAB: z(1,ind) = z(1,ind-1) + b(ind-1) + e(ind-1)
    z = [0.0] * len(b)
    for ind in range(1, len(b)):
        z[ind] = z[ind - 1] + b[ind - 1] + e[ind - 1]

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
