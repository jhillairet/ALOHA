"""
Utility functions for ALOHA.

This module contains utility functions used throughout the ALOHA package,
"""

import pathlib
import re
from typing import Any, Union

import h5py
import numpy as np
import scipy.io


def parse_matlab_array(array_str: str) -> list[str]:
    """
    Parse a MATLAB array string into a Python list.

    Handles both comma-separated and space-separated arrays.

    Parameters
    ----------
    array_str : str
        String containing MATLAB array elements (e.g., "1 2 3" or "1, 2, 3").

    Returns
    -------
    list[str]
        Python list of parsed values.
    """
    # Remove any whitespace and try to split by commas first
    array_str = array_str.strip()
    if "," in array_str:
        # Comma-separated array
        parts = [x.strip() for x in array_str.split(",")]
    else:
        # Space-separated array
        parts = array_str.split()
    return parts


def load_mat_file(filepath: pathlib.Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Load scenario from a MATLAB .mat binary file.

    Parameters
    ----------
    filepath : pathlib.Path
        Path to the MATLAB .mat binary file.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        A tuple containing:
        - scenario_dict: Dictionary containing the scenario data.
        - results_dict: Dictionary containing the results data.

    Raises
    ------
    ValueError
        If no scenario structure is found in the .mat file.
    """
    # Try scipy.io.loadmat first (for v7 and earlier)
    try:
        mat_data = scipy.io.loadmat(str(filepath))

        # Check if the .mat file contains a scenario structure
        if "scenario" in mat_data:
            scenario_struct = mat_data["scenario"]
            return parse_matlab_structure(scenario_struct)
    except NotImplementedError:
        # v7.3 files require h5py
        pass

    # Try h5py for v7.3 files
    try:
        with h5py.File(filepath, "r") as f:
            if "scenario" in f:
                # Parse the h5py group directly
                return parse_h5py_scenario(f["scenario"])
    except Exception:
        pass

    raise ValueError(f"No scenario structure found in .mat file: {filepath}")


def parse_h5py_scenario(scenario_group) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Parse an h5py scenario group into Python dictionaries.

    Parameters
    ----------
    scenario_group : h5py.Group
        h5py Group containing the scenario data.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        A tuple containing:
        - scenario_dict: Dictionary containing the scenario data.
        - results_dict: Dictionary containing the results data.
    """
    scenario_dict = {}
    results_dict = {}

    for key in scenario_group:
        item = scenario_group[key]

        if key.lower() == "results":
            results_dict = _h5py_group_to_dict(item)
        else:
            scenario_dict[key.lower()] = _h5py_group_to_dict(item)

    return scenario_dict, results_dict


def _h5py_group_to_dict(group):
    """
    Convert an h5py Group to a dictionary-like structure.

    Parameters
    ----------
    group : h5py.Group
        h5py Group to convert.

    Returns
    -------
    dict
        Dictionary representation of the h5py Group.
    """
    result = {}
    for key in group:
        item = group[key]
        if isinstance(item, h5py.Group):
            result[key] = _h5py_group_to_dict(item)
        elif isinstance(item, h5py.Dataset):
            # Convert dataset to Python value
            data = item[()]
            if data.size == 1:
                result[key] = data.item()
            else:
                result[key] = data.tolist() if data.dtype != object else data
        else:
            result[key] = item
    return result


def parse_matlab_structure(scenario_struct) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Parse a MATLAB structure array into Python dictionaries.

    Parameters
    ----------
    scenario_struct : numpy.ndarray
        MATLAB structure array containing scenario data.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        A tuple containing:
        - scenario_dict: Dictionary containing the scenario data.
        - results_dict: Dictionary containing the results data.
    """
    # Extract the first scenario (assuming it's a 1x1 struct array)
    if scenario_struct.ndim == 2 and scenario_struct.shape == (1, 1):
        scenario = scenario_struct[0, 0]
    elif scenario_struct.ndim == 1 and scenario_struct.shape[0] == 1:
        scenario = scenario_struct[0]
    else:
        # Handle array of scenarios - take the first one
        scenario = scenario_struct[0]

    scenario_dict = {}
    results_dict = {}

    # Process each field in the MATLAB structure
    for field_name in scenario.dtype.names:
        field_value = scenario[field_name][0]

        if field_name.lower() == "results":
            results_dict = convert_matlab_value(field_value)
        else:
            scenario_dict[field_name.lower()] = convert_matlab_value(field_value)

    return scenario_dict, results_dict


def convert_matlab_value(value: Any) -> Any:
    """
    Convert a MATLAB value to a Python value.

    Parameters
    ----------
    value : Any
        MATLAB value to convert (e.g., numpy.ndarray, numpy.generic, list, tuple).

    Returns
    -------
    Any
        Python value equivalent of the MATLAB value.
    """
    if isinstance(value, np.ndarray):
        if value.dtype.names:
            # Structure array
            return {name: convert_matlab_value(value[name][0]) for name in value.dtype.names}
        elif value.ndim == 0:
            # Scalar
            return value.item()
        elif value.ndim == 1:
            # Vector
            return value.tolist()
        elif value.ndim == 2:
            # Matrix
            return value.tolist()
        else:
            return value.tolist()
    elif isinstance(value, (list, tuple)):
        return [convert_matlab_value(v) for v in value]
    elif isinstance(value, np.generic):
        return value.item()
    else:
        return value


def load_m_file(filepath: pathlib.Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Load scenario from a MATLAB .m script file.

    Parameters
    ----------
    filepath : pathlib.Path
        Path to the MATLAB .m script file.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        A tuple containing:
        - scenario_dict: Dictionary containing the scenario data.
        - results_dict: Dictionary containing the results data.

    Raises
    ------
    ValueError
        If the scenario cannot be parsed from the .m file.
    """
    with open(filepath, encoding="utf-8") as f:
        content = f.read()

    scenario_dict = parse_matlab_script(content)

    if not scenario_dict:
        raise ValueError(f"Could not parse scenario from .m file: {filepath}")

    return scenario_dict, {}


def parse_matlab_script(content: str) -> dict[str, Any]:
    """
    Parse a MATLAB script file content into scenario and results dictionaries.

    This function parses MATLAB .m files that define ALOHA scenarios, extracting
    variable assignments and structure definitions into Python dictionaries.

    Parameters
    ----------
    content : str
        Content of the MATLAB script file.

    Returns
    -------
    dict[str, Any]
        Dictionary containing all the scenario data from the MATLAB file.
    """
    scenario_dict = {}

    # Remove comments (lines starting with %)
    lines = [line for line in content.split("\n") if not line.strip().startswith("%")]

    # Process each line
    for line in lines:
        line = line.strip()
        if not line or line.startswith("function") or line.startswith("end"):
            continue

        # Skip lines that are just continuation (start with ...)
        if line.startswith("..."):
            continue

        # Handle variable assignments
        if "=" in line:
            # Split on first '=' only
            var_part, value_part = line.split("=", 1)
            var_part = var_part.strip()
            value_part = value_part.strip().rstrip(";")

            # Parse the value
            value = _parse_matlab_value(value_part)

            # Handle nested structure assignments (e.g., modules.nma_theta)
            if "." in var_part:
                parts = var_part.split(".")
                current_dict = scenario_dict

                # Navigate to the parent structure
                for part in parts[:-1]:
                    if part not in current_dict:
                        current_dict[part] = {}
                    current_dict = current_dict[part]

                # Set the final value
                current_dict[parts[-1]] = value
            else:
                # Top-level variable
                scenario_dict[var_part] = value

    return scenario_dict


def _parse_matlab_value(value_str: str) -> Any:
    """
    Parse a MATLAB value string into a Python value.

    Parameters
    ----------
    value_str : str
        MATLAB value string to parse.

    Returns
    -------
    Any
        Python value equivalent of the MATLAB value.
    """
    value_str = value_str.strip()

    # Remove MATLAB comments (everything after %)
    if "%" in value_str:
        value_str = value_str.split("%")[0].strip()

    # Remove trailing semicolons
    if value_str.endswith(";"):
        value_str = value_str[:-1].strip()

    # Handle empty arrays
    if value_str == "[]":
        return []

    # Handle strings (single quotes)
    if value_str.startswith("'") and value_str.endswith("'"):
        return value_str[1:-1]

    # Handle numeric values
    try:
        return int(value_str)
    except ValueError:
        try:
            return float(value_str)
        except ValueError:
            pass

    # Handle MATLAB functions
    if value_str.startswith("ones(") or value_str.startswith("zeros("):
        return value_str

    if value_str.startswith("repmat("):
        return value_str

    if value_str.startswith("linspace("):
        return value_str

    if value_str.startswith("mfilename"):
        return value_str

    if value_str.startswith("pwd"):
        return value_str

    # Handle arrays (e.g., [1, 2, 3] or [1 2 3])
    if value_str.startswith("[") and value_str.endswith("]"):
        inner = value_str[1:-1].strip()
        if "," in inner:
            elements = [x.strip() for x in inner.split(",")]
        else:
            elements = inner.split()

        parsed_elements = []
        for elem in elements:
            if elem:
                parsed_elements.append(_parse_matlab_value(elem))

        return parsed_elements

    # Handle colon operator (e.g., 1:8)
    if ":" in value_str and not value_str.startswith("'") and not value_str.startswith('"'):
        return value_str

    # Handle variables (e.g., modules.nma_phi)
    if "." in value_str and not value_str.startswith("'") and not value_str.startswith('"'):
        return value_str

    # If we can't parse it, return the string as-is
    return value_str


def convert_matlab_scenario(matlab_scenario: dict[str, Any]) -> dict[str, Any]:
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
    if "antenna" in matlab_scenario and "architecture" in matlab_scenario["antenna"]:
        architecture = matlab_scenario["antenna"]["architecture"]
        # Map MATLAB antenna names to TOML file names
        antenna_name_mapping = {
            "antenne_elementaire": "simple_antenna.toml",
            # Add more mappings as needed
        }
        antenna["file"] = antenna_name_mapping.get(architecture, architecture)

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
    if "antenna" in matlab_scenario and "architecture" in matlab_scenario["antenna"]:
        architecture = matlab_scenario["antenna"]["architecture"]
        # Try to extract number from architecture name if it contains a number
        arch_numbers = re.findall(r"\d+", architecture)
        if arch_numbers:
            num_modules = int(arch_numbers[0])

    # antenna.a_ampl -> antenna.excitation.power
    if "antenna" in matlab_scenario and "a_ampl" in matlab_scenario["antenna"]:
        a_ampl = matlab_scenario["antenna"]["a_ampl"]
        if isinstance(a_ampl, (list, np.ndarray)):
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
