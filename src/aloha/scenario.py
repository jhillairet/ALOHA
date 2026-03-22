import os
import pathlib
import tomllib


class Scenario:
    def __init__(self, scenario: str | os.PathLike | dict | None = None):

        self.scenario = {}
        if isinstance(scenario, (str, os.PathLike)):
            self.scenario = self.load(scenario)
        elif isinstance(scenario, dict):
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
