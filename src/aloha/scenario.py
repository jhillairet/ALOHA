import os
import pathlib

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from aloha.utils import convert_matlab_scenario, load_m_file, load_mat_file


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
    def from_matlab(cls, filename: str | os.PathLike) -> tuple[dict, dict]:
        """
        Load scenario from an ALOHA MATLAB file (.m or .mat).

        Parameters
        ----------
        filename : str | os.PathLike
            Path to a MATLAB file (.m script or .mat binary file).

        Returns
        -------
        tuple[dict, dict]
            A tuple containing:
            - scenario_dict: Dictionary following the TOML schema used in scenario_example.toml
            - results_dict: Dictionary containing results if available, otherwise empty dict

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
        scenario_dict = convert_matlab_scenario(matlab_scenario)
        return scenario_dict, results_dict

    def convert_matlab_scenario(self, matlab_scenario: dict) -> dict:
        """
        Convert the deprecated matlab scenario schema into the new one.

        Parameters
        ----------
        matlab_scenario : dict
            deprecated Matlab scenario parameters.

        Returns
        -------
        scenario_dict : dict
            ALOHA scenario parameters.

        """
        return convert_matlab_scenario(matlab_scenario)

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
