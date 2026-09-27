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

    # Handle MATLAB booleans
    if value_str.lower() == "true":
        return True
    if value_str.lower() == "false":
        return False

    # Handle MATLAB functions (raise error for unsupported functions)
    if value_str.startswith("ones(") or value_str.startswith("zeros("):
        raise ValueError(f"MATLAB function '{value_str}' is not supported. Please provide explicit arrays.")

    if value_str.startswith("repmat("):
        raise ValueError(f"MATLAB function '{value_str}' is not supported. Please provide explicit arrays.")

    if value_str.startswith("linspace("):
        raise ValueError(f"MATLAB function '{value_str}' is not supported. Please provide explicit arrays.")

    if value_str.startswith("mfilename"):
        return value_str

    if value_str.startswith("pwd"):
        return value_str

    # Handle arrays (e.g., [1, 2, 3] or [1 2 3])
    # Also handle MATLAB transpose operator (') at the end
    if value_str.startswith("[") and (value_str.endswith("]") or value_str.endswith("]'")):
        # Strip the transpose operator if present
        if value_str.endswith("]'"):
            inner = value_str[1:-2].strip()
        else:
            inner = value_str[1:-1].strip()

        if "," in inner:
            elements = [x.strip() for x in inner.split(",")]
        else:
            elements = inner.split()

        parsed_elements = []
        for elem in elements:
            if elem:
                # Handle expressions like 1/8, etc.
                if any(op in elem for op in ["+", "-", "*", "/", "^"]):
                    # Try to evaluate the expression safely
                    try:
                        # Replace MATLAB-specific constants
                        elem_clean = elem.replace("pi", str(np.pi))
                        # Evaluate the expression
                        parsed_elements.append(float(eval(elem_clean, {"__builtins__": {}}, {})))
                    except (ValueError, SyntaxError, NameError):
                        # If evaluation fails, raise an error
                        raise ValueError(f"Cannot parse MATLAB expression: {elem}")
                else:
                    parsed_elements.append(_parse_matlab_value(elem))

        return parsed_elements

    # Handle colon operator (e.g., 1:8)
    if ":" in value_str and not value_str.startswith("'") and not value_str.startswith('"'):
        return value_str

    # Handle variables (e.g., modules.nma_phi)
    if "." in value_str and not value_str.startswith("'") and not value_str.startswith('"'):
        return value_str

    # Handle expressions involving arrays (e.g., (pi/180)*[0 -90*1 -90*2 ...])
    # Check if the string contains an array and operators
    if "[" in value_str and "]" in value_str and any(op in value_str for op in ["+", "-", "*", "/", "^"]):
        # Extract the array part, including the transpose operator if present
        array_start = value_str.find("[")
        array_end = value_str.rfind("]")
        # Check if there's a transpose operator after the array
        if array_end + 1 < len(value_str) and value_str[array_end + 1] == "'":
            array_str = value_str[array_start : array_end + 2]  # Include the transpose operator
        else:
            array_str = value_str[array_start : array_end + 1]

        # Parse the array
        try:
            array_value = _parse_matlab_value(array_str)
        except ValueError:
            raise ValueError(f"Cannot parse MATLAB array expression: {value_str}")

        if isinstance(array_value, list):
            # Extract the scalar multiplier or expression (the part outside the array)
            prefix = value_str[:array_start].strip()
            # Skip the array and transpose operator (if present) for the suffix
            suffix_start = (
                array_end + 2 if (array_end + 1 < len(value_str) and value_str[array_end + 1] == "'") else array_end + 1
            )
            suffix = value_str[suffix_start:].strip()

            # Combine prefix and suffix for the scalar part
            scalar_expr = (prefix + suffix).strip()

            # Remove any trailing operators (e.g., '*' or '/' at the end)
            scalar_expr = scalar_expr.rstrip("* / + - ^").strip()

            # Replace MATLAB-specific constants in the scalar expression
            scalar_expr_clean = scalar_expr.replace("pi", str(np.pi))

            # Evaluate the scalar expression
            try:
                scalar = float(eval(scalar_expr_clean, {"__builtins__": {}}, {}))
            except (ValueError, SyntaxError, NameError):
                raise ValueError(f"Cannot parse MATLAB scalar expression: {scalar_expr}")

            # Multiply the array by the scalar
            return [scalar * x for x in array_value]
        else:
            raise ValueError(f"Cannot parse MATLAB array expression: {value_str}")

    # If we can't parse it, return the string as-is
    return value_str
