import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from .constants import pi
from .utils import parse_matlab_array
from .waveguide import Waveguide


class Antenna:
    """
    ALOHA antenna description.

    This class represents an antenna for the ALOHA (Advanced LOwer Hybrid Antenna)
    simulation framework. It can be created from TOML files, MATLAB files,
    or dictionaries containing antenna parameters.
    """

    def __init__(self, source: str | os.PathLike | dict | None = None) -> None:
        """
        Initialize an Antenna instance.

        Parameters
        ----------
        source : str, Path, dict, or None, optional
            Source of antenna data. Can be:
            - A path to a TOML ALOHA antenna file
            - A dictionary containing antenna parameters
            - None to create an empty antenna
        """
        self.antenna = {}
        if source is None:
            return
        elif isinstance(source, (str, os.PathLike)):
            self.antenna = self.load(source)
        elif isinstance(source, dict):
            self.antenna = source
        else:
            raise TypeError("source must be a string, PathLike, dictionary, or None")

    @classmethod
    def load(cls, filename: str | os.PathLike) -> dict:
        """
        Load antenna parameters from a TOML file.

        Parameters
        ----------
        filename : str or Path
            Path to a TOML file containing antenna parameters.

        Returns
        -------
        dict
            Dictionary containing antenna parameters.
        """
        with open(filename, "rb") as fp:
            antenna_data = cls.validate_description(tomllib.load(fp))
        return antenna_data

    @classmethod
    def from_dict(cls, antenna_dict: dict) -> "Antenna":
        """
        Create an Antenna instance from a dictionary.

        Parameters
        ----------
        antenna_dict : dict
            Dictionary containing antenna parameters.

        Returns
        -------
        Antenna
            An instance of the Antenna class.
        """
        _ant = cls()
        _ant.antenna = cls.validate_description(antenna_dict)
        return _ant

    @classmethod
    def from_file(cls, filename: str | os.PathLike) -> "Antenna":
        """
        Create an Antenna instance from a TOML file.

        Parameters
        ----------
        filename : str or Path
            Path to a TOML file containing antenna parameters.

        Returns
        -------
        Antenna
            An instance of the Antenna class.
        """
        antenna_data = cls.load(filename)
        return cls.from_dict(antenna_data)

    @classmethod
    def from_matlab(cls, filename: str | os.PathLike) -> "Antenna":
        """
        Create an Antenna instance from a MATLAB antenna file.

        Parameters
        ----------
        filename : str or Path
            Path to a MATLAB ALOHA antenna file (.m).

        Returns
        -------
        Antenna
            An instance of the Antenna class.
        """
        # Read the MATLAB file
        with open(filename) as f:
            matlab_content = f.read()

        # Initialize antenna dictionary with Python structure
        antenna_dict = {}

        # Extract basic antenna properties
        name_match = re.search(r'antenna_lh\.name\s*=\s*[\'"]([^\'"]+)[\'"]', matlab_content)
        if name_match:
            antenna_dict["name"] = name_match.group(1)

        freq_match = re.search(r"antenna_lh\.frequency\s*=\s*([\d.eE+-]+)", matlab_content)
        if freq_match:
            antenna_dict["frequency"] = float(freq_match.group(1))

        # Extract type. Matlab version only allowed LH antenna
        antenna_dict["type"] = "LH"

        # Extract layout information
        layout = {}

        # Number of modules in toroidal direction (nma_phi)
        nb_mod_phi_match = re.search(r"modules\.nma_phi\s*=\s*(\d+)", matlab_content)
        if nb_mod_phi_match:
            layout["nb_mod_phi"] = int(nb_mod_phi_match.group(1))

        # Number of modules in poloidal direction (nma_theta)
        nb_mod_theta_match = re.search(r"modules\.nma_theta\s*=\s*(\d+)", matlab_content)
        if nb_mod_theta_match:
            layout["nb_mod_theta"] = int(nb_mod_theta_match.group(1))

        # Spacing between modules in toroidal direction (sm_phi)
        spacing_phi_match = re.search(r"modules\.sm_phi\s*=\s*([\d.eE+-]+)", matlab_content)
        if spacing_phi_match:
            layout["spacing_btw_mod_phi"] = float(spacing_phi_match.group(1))
        else:
            layout["spacing_btw_mod_phi"] = 0

        # Spacing between modules in poloidal direction (sm_theta)
        spacing_theta_match = re.search(r"modules\.sm_theta\s*=\s*([\d.eE+-]+)", matlab_content)
        if spacing_theta_match:
            layout["spacing_btw_mod_theta"] = float(spacing_theta_match.group(1))
        else:
            layout["spacing_btw_mod_theta"] = 0

        antenna_dict["layout"] = layout

        # Extract module/waveguide information
        module = {}

        # Number of waveguides per module in toroidal direction (nwm_phi)
        nwm_phi_match = re.search(r"waveguides\.nwm_phi\s*=\s*(\d+)", matlab_content)
        if nwm_phi_match:
            module["nb_wg_phi"] = int(nwm_phi_match.group(1))

        # Number of waveguides per module in poloidal direction (nwm_theta)
        nwm_theta_match = re.search(r"waveguides\.nwm_theta\s*=\s*(\d+)", matlab_content)
        if nwm_theta_match:
            module["nb_wg_theta"] = int(nwm_theta_match.group(1))

        # Mask of active/passive waveguides
        mask_match = re.search(r"waveguides\.mask\s*=\s*\[([^\]]+)\]", matlab_content)
        if mask_match:
            mask_str = mask_match.group(1)
            mask_parts = parse_matlab_array(mask_str)
            mask = [int(x.strip()) for x in mask_parts]
            module["mask"] = mask

        # Number of passive waveguides between modules (npwbm_phi)
        npwbm_phi_match = re.search(r"waveguides\.npwbm_phi\s*=\s*(\d+)", matlab_content)
        if npwbm_phi_match:
            module["nb_pwg_btw_mod_phi"] = int(npwbm_phi_match.group(1))
        else:
            module["nb_pwg_btw_mod_phi"] = 0

        # Number of passive waveguides on each edge (npwe_phi)
        npwe_phi_match = re.search(r"waveguides\.npwe_phi\s*=\s*(\d+)", matlab_content)
        if npwe_phi_match:
            module["nb_pwg_edge"] = int(npwe_phi_match.group(1))
        else:
            module["nb_pwg_edge"] = 0

        # Spacing between waveguides in poloidal direction (sw_theta)
        sw_theta_match = re.search(r"waveguides\.sw_theta\s*=\s*([\d.eE+-]+)", matlab_content)
        if sw_theta_match:
            module["space_btw_wg_theta"] = float(sw_theta_match.group(1))

        # Height of waveguides in poloidal direction (hw_theta)
        hw_theta_match = re.search(r"waveguides\.hw_theta\s*=\s*([\d.eE+-]+)", matlab_content)
        if hw_theta_match:
            module["wg_size_theta"] = float(hw_theta_match.group(1))

        # Width of active waveguides (bwa)
        bwa_match = re.search(r"waveguides\.bwa\s*=\s*([\d.eE+-]+)", matlab_content)
        if bwa_match:
            module["awg_size_phi"] = float(bwa_match.group(1))

        # Width of internal passive waveguides (biwp)
        biwp_match = re.search(r"waveguides\.biwp\s*=\s*([\d.eE+-]+)", matlab_content)
        if biwp_match:
            module["pwg_size_phi"] = float(biwp_match.group(1))

        # Width of edge passive waveguides (bewp)
        bewp_match = re.search(r"waveguides\.bewp\s*=\s*([\d.eE+-]+)", matlab_content)
        if bewp_match:
            module["pwg_size_edge_phi"] = float(bewp_match.group(1))

        # Thickness between waveguides (e_phi)
        e_phi_match = re.search(r"waveguides\.e_phi\s*=\s*\[([^\]]+)\]", matlab_content)
        if e_phi_match:
            e_phi_str = e_phi_match.group(1)
            e_phi_parts = parse_matlab_array(e_phi_str)
            e_phi_values = [float(x.strip()) for x in e_phi_parts]
            if len(e_phi_values) == 1:
                module["e_phi"] = e_phi_values[0]
            else:
                module["e_phi"] = e_phi_values
        else:
            # Try to find ep variable assignment
            ep_match = re.search(r"\be_p\s*=\s*([\d.eE+-]+)", matlab_content)
            if ep_match:
                module["e_phi"] = float(ep_match.group(1))

        # Thickness between passive waveguides (e_phi_pwg)
        e_phi_pwg_match = re.search(r"waveguides\.e_phi_pwg\s*=\s*([\d.eE+-]+)", matlab_content)
        if not e_phi_pwg_match:
            e_phi_pwg_match = re.search(r"\be_phi_pwg\s*=\s*([\d.eE+-]+)", matlab_content)
        if e_phi_pwg_match:
            module["e_phi_pwg"] = float(e_phi_pwg_match.group(1))
        else:
            if "e_phi" in module and isinstance(module["e_phi"], (int, float)):
                module["e_phi_pwg"] = module["e_phi"]
            else:
                module["e_phi_pwg"] = 2e-3

        # Short circuit depth for passive waveguides (scl)
        scl_match = re.search(r"waveguides\.scl\s*=\s*\[([^\]]+)\]", matlab_content)
        if scl_match:
            scl_str = scl_match.group(1)
            scl_parts = parse_matlab_array(scl_str)
            scl_values = [float(x.strip()) for x in scl_parts]
            if len(scl_values) == 1:
                module["pwg_depth"] = scl_values[0]
            else:
                module["pwg_depth"] = scl_values
        else:
            module["pwg_depth"] = 0.25

        antenna_dict["module"] = module

        # Extract S-parameters information
        sparameters = {}

        # S-parameter filenames
        filenames_match = re.search(r"modules\.Sparameters\.SFileNames\s*=\s*\[(.*?)\]", matlab_content, re.DOTALL)
        if not filenames_match:
            filenames_match = re.search(r"SFileNames\s*=\s*\[(.*?)\]", matlab_content, re.DOTALL)
        if not filenames_match:
            filenames_match = re.search(r"filenames\s*=\s*\[(.*?)\]", matlab_content, re.DOTALL)
        if filenames_match:
            filenames_str = filenames_match.group(1)
            filenames_str = filenames_str.replace("'", "").replace('"', "")
            # Handle both comma-separated and semicolon-separated lists (MATLAB strvcat format)
            if ";" in filenames_str:
                filenames_parts = [x.strip() for x in filenames_str.split(";")]
            else:
                filenames_parts = parse_matlab_array(filenames_str)
            filenames = [x.strip() for x in filenames_parts if x.strip()]
            sparameters["filenames"] = filenames

        # Phase deembedding
        phase_deembedded_match = re.search(r"phase_deembedded\s*=\s*\[([^\]]+)\]", matlab_content)
        if phase_deembedded_match:
            phase_str = phase_deembedded_match.group(1)
            phase_parts = parse_matlab_array(phase_str)
            phases = [float(x.strip()) for x in phase_parts]
            sparameters["phases_deembedded"] = phases

        # Extract excitation information
        excitation = {}

        # Magnitudes (amplitude)
        amplitude_match = re.search(r"modules\.amplitude\s*=\s*\[([^\]]+)\]", matlab_content)
        if amplitude_match:
            amp_str = amplitude_match.group(1)
            amp_parts = parse_matlab_array(amp_str)
            try:
                magnitudes = [float(x.strip()) for x in amp_parts]
                excitation["magnitudes"] = magnitudes
            except ValueError:
                # Skip if it's a MATLAB function call like zeros(), ones(), etc.
                pass

        # Phases - try different patterns (skip commented lines)
        # Remove commented lines first to avoid matching them
        matlab_no_comments = re.sub(r"%.+", "", matlab_content)
        phase_match = re.search(r"modules\.phase\s*=\s*([^\n]+)", matlab_no_comments)
        a_phase_match = re.search(r"modules\.a_phase\s*=\s*([^\n]+)", matlab_no_comments)

        if phase_match:
            phase_str = phase_match.group(1).strip().rstrip(";")
            # Skip MATLAB function calls like zeros(), ones(), etc.
            if any(phase_str.startswith(func) for func in ["zeros(", "ones(", "repmat(", "linspace("]):
                pass  # Skip MATLAB function calls
            # Handle expressions like (270*pi/180)*(0:modules.nma_phi-1)
            elif "*" in phase_str and "pi" in phase_str and ":" in phase_str:
                # Extract the multiplier and range
                multiplier_match = re.search(r"([\d.eE+-]+)\s*\*\s*pi\s*/\s*([\d.eE+-]+)", phase_str)
                # Look for the range pattern specifically (with colon)
                range_match = re.search(r"\(([^)]*:[^)]*)\)", phase_str)
                if multiplier_match and range_match:
                    multiplier = float(multiplier_match.group(1)) * np.pi / float(multiplier_match.group(2))
                    range_expr = range_match.group(1)
                    # Parse range like "0:modules.nma_phi-1" or "0:7"
                    range_parts = range_expr.split(":")
                    start = int(range_parts[0].strip())
                    end_str = range_parts[1].strip()
                    # Handle expressions like "modules.nma_phi-1"
                    if "modules.nma_phi" in end_str:
                        nb_modules = layout.get("nb_mod_phi", 1)
                        end = nb_modules - 1
                    else:
                        end = int(end_str)
                    phases = [multiplier * i for i in range(start, end + 1)]
                    excitation["phases"] = phases
            else:
                # Handle simple array
                if "pi" in phase_str:
                    phase_str = phase_str.replace("pi", str(np.pi))
                phase_parts = parse_matlab_array(phase_str)
                try:
                    phases = [float(x.strip()) for x in phase_parts]
                    excitation["phases"] = phases
                except ValueError:
                    # Skip if it's a MATLAB function call
                    pass
        elif a_phase_match:
            phase_str = a_phase_match.group(1).strip()
            # Skip MATLAB function calls
            if any(phase_str.startswith(func) for func in ["zeros(", "ones(", "repmat(", "linspace("]):
                pass
            elif "pi" in phase_str:
                phase_str = phase_str.replace("pi", str(np.pi))
                phase_parts = parse_matlab_array(phase_str)
                try:
                    phases = [float(x.strip()) for x in phase_parts]
                    excitation["phases"] = phases
                except ValueError:
                    pass

        if excitation:
            antenna_dict["excitation"] = excitation
        else:
            # Provide default excitation values
            total_nb_modules = layout.get("nb_mod_phi", 1) * layout.get("nb_mod_theta", 1)
            if total_nb_modules > 0:
                antenna_dict["excitation"] = {
                    "magnitudes": [1.0] * total_nb_modules,
                    "phases": [0.0] * total_nb_modules,
                }

        # Ensure sparameters has required fields
        if "filenames" not in sparameters:
            total_nb_modules = layout.get("nb_mod_phi", 1) * layout.get("nb_mod_theta", 1)
            if total_nb_modules > 0:
                sparameters["filenames"] = ["S_elem"] * total_nb_modules
        if "phases_deembedded" not in sparameters:
            total_nb_modules = layout.get("nb_mod_phi", 1) * layout.get("nb_mod_theta", 1)
            if total_nb_modules > 0:
                sparameters["phases_deembedded"] = [0.0] * total_nb_modules

        antenna_dict["sparameters"] = sparameters

        # Create and return Antenna instance
        return cls.from_dict(antenna_dict)

    def to_toml(self, filename: str | os.PathLike | None = None) -> str | None:
        """
        Export the antenna description to TOML format.

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
            elif isinstance(value, (int, float)):
                return str(value)
            elif isinstance(value, list):
                if len(value) == 0:
                    return "[]"
                # Check if all elements are numeric
                if all(isinstance(x, (int, float)) for x in value):
                    return "[" + ", ".join(str(x) for x in value) + "]"
                else:
                    # Handle string lists
                    return "[" + ", ".join(f'"{x}"' for x in value) + "]"
            elif isinstance(value, bool):
                return str(value).lower()
            else:
                return str(value)

        def format_section(section_name: str, data: dict, comments: dict = None) -> str:
            """Format a TOML section with comments."""
            if comments is None:
                comments = {}

            lines = []
            if section_name:
                lines.append(f"[{section_name}]")
                lines.append("")

            for key, value in data.items():
                # Add comment if available
                if comments and key in comments:
                    lines.append(f"# {comments[key]}")

                if isinstance(value, dict):
                    # Handle nested sections (though TOML doesn't support nested sections directly)
                    # For now, we'll flatten this case
                    for subkey, subvalue in value.items():
                        lines.append(f"{subkey} = {format_value(subvalue)}")
                else:
                    lines.append(f"{key} = {format_value(value)}")
                lines.append("")  # Add blank line between entries

            # Remove the last blank line if the section is not empty
            if lines and lines[-1] == "":
                lines = lines[:-1]

            return "\n".join(lines)

        # Define comments for each field
        comments = {
            "type": 'Antenna type ("LH" or "IC")',
            "name": "Antenna name",
            "frequency": "Antenna default frequency [Hz]",
            "nb_mod_phi": "Number of modules per antenna in the toroidal direction.",
            "nb_mod_theta": "Number of modules per antenna in the poloidal direction.",
            "spacing_btw_mod_phi": "Spacing between toroidally neighboring modules [m]",
            "spacing_btw_mod_theta": "Spacing between poloidally neighboring modules [m]",
            "nb_wg_phi": "Number of waveguides per module in the toroidal direction. (passive and active)",
            "nb_wg_theta": "Number of waveguides per module in the poloidal direction. (passive and active)",
            "mask": "Mask of passive and active waveguides for an internal module. 1 for active -- 0 for passive.",
            "nb_pwg_btw_mod_phi": "Number of passive waveguide between modules in the toroidal direction.",
            "nb_pwg_edge": "Number of passive waveguides on each antenna edge in the toroidal direction.",
            "space_btw_wg_theta": "Spacing between poloidally neighboring waveguides [m]",
            "wg_size_theta": "Height of waveguides in the poloidal direction [m]",
            "awg_size_phi": "Width of active waveguides [m]",
            "pwg_size_phi": "Width of internal passive waveguides [m]",
            "pwg_size_edge_phi": "Width of edge passive waveguides [m]",
            "e_phi": "Thickness between waveguides in the toroidal direction [m]",
            "e_phi_pwg": "Thickness between passive waveguides",
            "pwg_depth": "Short circuit depth for passive waveguides in guided wavelength",
            "filenames": "S-parameters file for each module",
            "phases_deembedded": "Phase deembedding for each module inputs",
            "magnitudes": "Default forward power for all modules [watt]",
            "phases": "Default phase shifts for all modules [deg]",
        }

        # Start building the TOML content
        lines = []

        # Add header comment
        if "name" in self.antenna:
            lines.append(f"# {self.antenna['name']}")
            lines.append("")

        # Add type if present
        if "type" in self.antenna:
            lines.append('# Antenne type ("LH" or "IC")')
            lines.append(f"type = {format_value(self.antenna['type'])}")
            lines.append("")

        # Add name and frequency
        if "name" in self.antenna:
            lines.append("# Antenna name")
            lines.append(f"name = {format_value(self.antenna['name'])}")
            lines.append("")

        if "frequency" in self.antenna:
            lines.append("# Antenna default frequency [Hz]")
            lines.append(f"frequency = {format_value(self.antenna['frequency'])}")
            lines.append("")

        # Add layout section
        if "layout" in self.antenna:
            lines.append("")
            lines.append("# General description of the modules")
            lines.append("")
            layout_lines = format_section("layout", self.antenna["layout"], comments)
            lines.append(layout_lines)

        # Add module section
        if "module" in self.antenna:
            lines.append("")
            lines.append("# module description")
            lines.append("")
            module_lines = format_section("module", self.antenna["module"], comments)
            lines.append(module_lines)

        # Add sparameters section
        if "sparameters" in self.antenna:
            lines.append("")
            lines.append("# S-parameters description")
            lines.append("")
            sparameters_lines = format_section("sparameters", self.antenna["sparameters"], comments)
            lines.append(sparameters_lines)

        # Add excitation section
        if "excitation" in self.antenna:
            lines.append("")
            lines.append("# Antenna default power and phase excitations")
            lines.append("")
            excitation_lines = format_section("excitation", self.antenna["excitation"], comments)
            lines.append(excitation_lines)

        # Join all lines and clean up extra blank lines
        toml_content = "\n".join(lines)
        # Remove consecutive blank lines
        toml_content = re.sub(r"\n{3,}", "\n\n", toml_content)
        # Remove trailing whitespace
        toml_content = "\n".join(line.rstrip() for line in toml_content.split("\n"))

        if filename is not None:
            # Write to file
            with open(filename, "w") as f:
                f.write(toml_content)
            return None
        else:
            return toml_content

    def __str__(self):
        """
        Return a string representation of the Antenna.

        Returns
        -------
        str
            String representation showing antenna name and frequency.
        """
        return f"Antenna(name={self.antenna['name']}, frequency={self.antenna['frequency']})"

    @classmethod
    def validate_description(cls, antenna: dict) -> dict:
        """
        Check that the antenna description is valid.

        Parameters
        ----------
        antenna : dict
            ALOHA antenna description dictionary.

        Returns
        -------
        antenna : dict
            Verified antenna description
        """
        if not antenna:  # empty dict is OK
            return antenna

        # Use layout if available, otherwise fall back to global
        if "layout" in antenna:
            total_nb_modules = antenna["layout"]["nb_mod_phi"] * antenna["layout"]["nb_mod_theta"]
        elif "global" in antenna:
            total_nb_modules = antenna["global"]["nb_mod_phi"] * antenna["global"]["nb_mod_theta"]
        else:
            raise ValueError("Antenna description must have either 'layout' or 'global' section")
        # excitations (skip if not present)
        if "excitation" in antenna:
            if "magnitudes" in antenna["excitation"]:
                if len(antenna["excitation"]["magnitudes"]) != total_nb_modules:
                    raise ValueError(
                        f"Number of excitation magnitudes ({len(antenna['excitation']['magnitudes'])}) "
                        f"does not match total number of modules ({total_nb_modules})"
                    )
            if "phases" in antenna["excitation"]:
                if len(antenna["excitation"]["phases"]) != total_nb_modules:
                    raise ValueError(
                        f"Number of excitation phases ({len(antenna['excitation']['phases'])}) "
                        f"does not match total number of modules ({total_nb_modules})"
                    )

        # module indices (skip if not present, as it's not in TOML format)
        if "global" in antenna and "idx_mod" in antenna["global"]:
            if len(antenna["global"]["idx_mod"]) != total_nb_modules:
                raise ValueError(
                    f"Number of module indices ({len(antenna['global']['idx_mod'])}) "
                    f"does not match total number of modules ({total_nb_modules})"
                )

        # S-parameter (skip if not present)
        if "sparameters" in antenna:
            if "filenames" in antenna["sparameters"]:
                if len(antenna["sparameters"]["filenames"]) != total_nb_modules:
                    raise ValueError(
                        f"Number of S-parameter filenames ({len(antenna['sparameters']['filenames'])}) "
                        f"does not match total number of modules ({total_nb_modules})"
                    )
            if "phases_deembedded" in antenna["sparameters"]:
                if len(antenna["sparameters"]["phases_deembedded"]) != total_nb_modules:
                    raise ValueError(
                        f"Number of deembedded phases ({len(antenna['sparameters']['phases_deembedded'])}) "
                        f"does not match total number of modules ({total_nb_modules})"
                    )
        # so far so good
        return antenna

    def is_valid(self) -> bool:
        """
        Check if the antenna description is valid.

        Returns
        -------
        bool
            True if the antenna description is valid, False otherwise.
        """
        try:
            self.validate_description(self.antenna)
            return True
        except ValueError:
            return False

    def plot(self):
        """
        Plot the antenna architecture.

        Creates a visual representation of the antenna layout showing
        waveguide positions and dimensions.
        """
        if not self.is_valid():
            raise ValueError("Invalid antenna description")

        b, a, z, y, nwr, nwas, act_module_tor = self.antenna_coordinates()
        # Extract parameters from antenna dictionary
        antenna = self.antenna
        plt.figure()
        # TODO : plot all poloidal modules if more than one
        for idx_pol in range(len(y)):
            for idx_tor in range(len(z)):
                rect_pos = [z[idx_tor], y[idx_pol], b[idx_tor], a[idx_pol]]
                # Create passive/active mask
                ar_modules = np.ones(antenna["global"]["nb_mod_phi"])
                ar_pa_mask = np.array(antenna["module"]["mask"])

                # Add passive waveguides between modules
                ar_pa_mask = np.concatenate([ar_pa_mask, np.zeros(antenna["module"]["nb_pwg_btw_mod_phi"])])

                pa_mask = np.kron(ar_modules, ar_pa_mask)

                # Remove last element if there are passive waveguides between modules
                if antenna["module"]["nb_pwg_btw_mod_phi"] > 0:
                    pa_mask = pa_mask[:-1]

                # Add passive waveguides at edges
                pa_mask = np.concatenate(
                    [np.zeros(antenna["module"]["nb_pwg_edge"]), pa_mask, np.zeros(antenna["module"]["nb_pwg_edge"])]
                )

                if pa_mask[idx_tor] == 0:  # passive wg
                    plt.gca().add_patch(
                        plt.Rectangle(
                            (rect_pos[0], rect_pos[1]),
                            rect_pos[2],
                            rect_pos[3],
                            facecolor=[0.8, 0.8, 0.8],
                            edgecolor="k",
                        )
                    )
                elif pa_mask[idx_tor] == 1:  # active wg
                    plt.gca().add_patch(
                        plt.Rectangle(
                            (rect_pos[0], rect_pos[1]), rect_pos[2], rect_pos[3], fill=False, facecolor=[0.8, 0.8, 0.8]
                        )
                    )

        plt.axis("equal")
        plt.xlabel("z [m]")
        plt.ylabel("y [m]")
        plt.title(f"ALOHA antenna: {antenna['name']}\n (view from the plasma)")
        plt.show()

    def antenna_coordinates(self):
        """
        Generate waveguide positions and dimensions from the antenna description.

        Returns
        -------
        b : list
            Waveguide widths (in toroidal direction).
        a : list
            Waveguide heights (in poloidal direction).
        z : list
            Waveguide toroidal positions.
        y : list
            Waveguide poloidal positions.
        nb_wg_per_row : int
            Number of waveguides per row.
        nb_wg_total : int
            Total number of waveguides per antenna.
        act_module_tor : list
            Mask for passive/active waveguides.

        """
        antenna = self.antenna
        mod = antenna["global"]
        wg = antenna["module"]

        # (total) number of waveguides per row
        # = (nb wg in a module) + 2*(nb ext wg) + (nb of wg between modules)
        nb_wg_per_row = (
            mod["nb_mod_phi"] * wg["nb_wg_phi"]
            + 2 * wg["nb_pwg_edge"]
            + (mod["nb_mod_phi"] - 1) * wg["nb_pwg_btw_mod_phi"]
        )

        # (total) number of waveguides per column
        nb_wg_per_col = mod["nb_mod_theta"] * wg["nb_wg_theta"]

        # total number of waveguides
        nb_wg_total = nb_wg_per_row * nb_wg_per_col

        # waveguide height - supposed constant for all the waveguides of the antenna
        a = np.full(nb_wg_total, wg["wg_size_theta"])

        # b
        # Make the array b which contains all the waveguide width of a row of waveguides
        b_module = np.where(wg["mask"], wg["awg_size_phi"], wg["pwg_size_phi"])
        b_edge = np.full(wg["nb_pwg_edge"], wg["pwg_size_edge_phi"])
        b_inter = np.full(wg["nb_pwg_btw_mod_phi"], wg["pwg_size_phi"])

        b = np.concatenate(
            [b_edge, np.tile(np.concatenate([b_module, b_inter]), mod["nb_mod_phi"] - 1), b_module, b_edge]
        )

        # e
        # Make the array e which contains all the waveguide septum width of a row of waveguides
        ne_phi = (
            wg["nb_pwg_btw_mod_phi"] * (mod["nb_mod_phi"] - 1)
            + wg["nb_pwg_edge"] * 2
            + mod["nb_mod_phi"] * wg["nb_wg_phi"]
            - 1
        )
        e = np.tile(wg["pwg_size_phi"], ne_phi)

        # z
        # Make the array z which contains all the waveguide positions in the toroidal direction
        z = np.zeros(nb_wg_per_row)
        for ind in range(1, nb_wg_per_row):
            z[ind] = z[ind - 1] + b[ind - 1] + e[ind - 1]

        # y
        # Make the array y which contains all the waveguide positions in the poloidal direction
        h = np.concatenate(
            [
                np.tile(
                    np.concatenate(
                        [np.full(wg["nb_wg_theta"], wg["awg_size_phi"]), np.full(1, mod["spacing_btw_mod_phi"])]
                    ),
                    mod["nb_mod_theta"] - 1,
                ),
                np.full(wg["nb_wg_theta"], wg["awg_size_phi"]),
            ]
        )
        y = np.zeros(nb_wg_per_col)
        for ind in range(1, nb_wg_per_col):
            y[ind] = y[ind - 1] + h[ind - 1] + a[ind]

        # index of active waveguides in a module
        act_module_tor = wg["mask"]

        return b, a, z, y, nb_wg_per_row, nb_wg_total, act_module_tor
