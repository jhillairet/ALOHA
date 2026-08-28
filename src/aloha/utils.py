"""
Utility functions for ALOHA.

This module contains utility functions used throughout the ALOHA package,
"""

import pathlib
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


def load_matlab_mat_file(filepath: pathlib.Path) -> tuple[dict[str, Any], dict[str, Any]]:
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


def load_matlab_scenario_file(filepath: pathlib.Path) -> tuple[dict[str, Any], dict[str, Any]]:
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

    # Try to parse the MATLAB script
    scenario_dict, results_dict = parse_matlab_script(content)

    if not scenario_dict:
        raise ValueError(f"Could not parse scenario from .m file: {filepath}")

    return scenario_dict, results_dict


def parse_matlab_script(content: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Parse a MATLAB script file content into scenario and results dictionaries.

    Parameters
    ----------
    content : str
        Content of the MATLAB script file.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        A tuple containing:
        - scenario_dict: Dictionary containing the scenario data.
        - results_dict: Dictionary containing the results data.
    """
    # TODO: Implement MATLAB script parsing logic
    scenario_dict = {}
    results_dict = {}
    return scenario_dict, results_dict
