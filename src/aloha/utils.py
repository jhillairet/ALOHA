"""
Utility functions for ALOHA.

This module contains utility functions used throughout the ALOHA package,
"""


def parse_matlab_array(array_str: str) -> list:
    """
    Parse a MATLAB array string into a Python list.

    Handles both comma-separated and space-separated arrays.

    Parameters
    ----------
    array_str : str
        String containing MATLAB array elements (e.g., "1 2 3" or "1, 2, 3")

    Returns
    -------
    list
        Python list of parsed values
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
