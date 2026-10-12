"""
S-parameter file loading and parsing utilities.
"""

import re
from pathlib import Path
from typing import Any

import numpy as np


def parse_matlab_complex_value(value_str: str) -> complex:
    """Parse a MATLAB complex value string into a Python complex number."""
    value_str = value_str.strip()
    value_str = value_str.replace("i", "j")
    value_str = re.sub(r"\s*([+\-])\s*", r"\1", value_str)

    try:
        return complex(eval(value_str, {"__builtins__": {}}, {}))
    except (ValueError, SyntaxError, NameError):
        if "j" in value_str:
            j_pos = value_str.index("j")
            complex_expr = value_str[:j_pos]
            split_pos = -1
            for i in range(1, len(complex_expr)):
                if complex_expr[i] in "+-":
                    split_pos = i
            if split_pos > 0:
                real_part = float(complex_expr[:split_pos])
                imag_part = float(complex_expr[split_pos:])
                return complex(real_part, imag_part)
            else:
                imag_part = float(complex_expr)
                return complex(0, imag_part)
        else:
            return complex(float(value_str), 0)


def parse_matlab_array(array_str: str) -> list:
    """Parse a MATLAB array string into a Python list of complex numbers."""
    array_str = array_str.strip()
    if array_str.startswith("["):
        array_str = array_str[1:]
    if array_str.endswith("]"):
        array_str = array_str[:-1]

    # Handle arrays that span multiple lines
    # First, remove all newlines and extra spaces
    array_str = array_str.replace("\n", " ").replace("  ", " ")

    # Split by commas to get individual elements
    # But we need to be careful with complex numbers that have '+' in them
    # Strategy: split by comma, then by whitespace, but handle complex numbers properly

    elements = []
    # First try splitting by commas
    potential_elems = re.split(r"\s*,\s*", array_str)

    for elem in potential_elems:
        elem = elem.strip()
        if not elem:
            continue

        # Check if this might be multiple elements separated by whitespace
        # If it contains '+' or '-' followed by 'i' or 'j', it's likely a complex number
        # Otherwise, try splitting by whitespace
        if "i" in elem or "j" in elem:
            # This is likely a complex number, parse it directly
            elements.append(parse_matlab_complex_value(elem))
        else:
            # Try splitting by whitespace
            sub_elems = elem.split()
            for sub_elem in sub_elems:
                if sub_elem:
                    elements.append(parse_matlab_complex_value(sub_elem))

    return elements


def parse_matlab_sparameter_file(filepath: Path) -> dict[str, Any]:
    """Parse a MATLAB S-parameter file (.m format) and extract S, f, Z matrices."""
    with open(filepath) as f:
        content = f.read()

    result = {"f": None, "S": None, "Z": None}

    # Extract frequency
    f_match = re.search(r"f\s*=\s*([^;]+);", content)
    if f_match:
        f_str = f_match.group(1).strip()
        try:
            result["f"] = float(f_str)
        except ValueError:
            try:
                result["f"] = float(eval(f_str.replace("pi", str(np.pi))))
            except:
                pass

    # Extract S matrix - try to find the entire S = [...] block
    # Handle MATLAB arrays with semicolons (e.g., [0,1; 1,0])
    # Use a non-greedy match that allows semicolons inside brackets
    s_match = re.search(r"S\s*=\s*\[([^\]]*?)\]\s*;", content, re.DOTALL)
    # If no match, try to find S(1,:,:) = [...] (3D array format from HFSS)
    if not s_match:
        s_match = re.search(r"S\(1,:,:\)\s*=\s*\[([^\]]*?)\]\s*;", content, re.DOTALL)
    if s_match:
        s_str = s_match.group(1)
        try:
            # Replace 'i' with 'j' for Python complex numbers
            s_str_clean = s_str.replace("i", "j")
            # Remove MATLAB line continuation markers (...)
            s_str_clean = s_str_clean.replace("...", "")
            # Remove newlines and extra spaces
            s_str_clean = s_str_clean.replace("\n", " ").replace("  ", " ").strip()
            # Remove any trailing semicolons
            s_str_clean = s_str_clean.rstrip(";")
            # Replace MATLAB semicolons with commas for Python
            s_str_clean = s_str_clean.replace(";", ",")
            # Try to evaluate as a numpy array
            s_array = np.array(eval(f"[{s_str_clean}]", {"__builtins__": {}}, {}), dtype=complex)
            n = int(np.sqrt(len(s_array)))
            if n * n == len(s_array):
                result["S"] = s_array.reshape((n, n))
            else:
                result["S"] = s_array
        except Exception as e:
            print(f"Warning: Could not parse S matrix with eval: {e}")
            # Try the old method
            try:
                s_elements = parse_matlab_array(f"[{s_str}]")
                s_array = np.array(s_elements, dtype=complex)
                n = int(np.sqrt(len(s_array)))
                if n * n == len(s_array):
                    result["S"] = s_array.reshape((n, n))
                else:
                    result["S"] = s_array
            except Exception as e2:
                print(f"Warning: Could not parse S matrix: {e2}")
                result["S"] = None

    # Extract Z matrix
    z_match = re.search(r"Z\s*=\s*\[([^\]]*)\]\s*;", content, re.DOTALL)
    if z_match:
        z_str = z_match.group(1)
        try:
            z_str_clean = z_str.replace("i", "j")
            z_str_clean = z_str_clean.replace("\n", " ").replace("  ", " ").strip()
            z_str_clean = z_str_clean.rstrip(";")
            z_array = np.array(eval(f"[{z_str_clean}]", {"__builtins__": {}}, {}), dtype=complex)
            result["Z"] = z_array
        except Exception as e:
            print(f"Warning: Could not parse Z matrix with eval: {e}")
            try:
                z_elements = parse_matlab_array(f"[{z_str}]")
                result["Z"] = np.array(z_elements, dtype=complex)
            except Exception as e2:
                print(f"Warning: Could not parse Z matrix: {e2}")
                result["Z"] = None

    return result


def load_sparameter_file(filepath: Path | str) -> dict[str, Any]:
    """Load an S-parameter file from various formats."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"S-parameter file not found: {filepath}")

    if filepath.suffix.lower() == ".m":
        return parse_matlab_sparameter_file(filepath)
    elif filepath.suffix.lower() == ".mat":
        try:
            import scipy.io

            mat_data = scipy.io.loadmat(str(filepath))
            result = {}
            if "S" in mat_data:
                result["S"] = mat_data["S"].astype(complex)
            if "f" in mat_data:
                result["f"] = mat_data["f"].flatten()[0] if mat_data["f"].size > 0 else None
            if "Z" in mat_data:
                result["Z"] = mat_data["Z"].astype(complex) if "Z" in mat_data else None
            return result
        except NotImplementedError:
            import h5py

            with h5py.File(filepath, "r") as f:
                result = {}
                if "S" in f:
                    result["S"] = f["S"][()].astype(complex)
                if "f" in f:
                    result["f"] = f["f"][()].flatten()[0] if f["f"].size > 0 else None
                if "Z" in f:
                    result["Z"] = f["Z"][()].astype(complex)
                return result
        except Exception as e:
            raise ValueError(f"Could not load .mat file: {e}")
    elif filepath.suffix.lower().endswith("sp"):
        return parse_matlab_sparameter_file(filepath)
    else:
        return parse_matlab_sparameter_file(filepath)
