import os
import tomllib

import numpy as np

from .constants import pi
from .waveguide import Waveguide


class Antenna:
    def __init__(self, source: str | os.PathLike | dict | None = None) -> None:
        """
        ALOHA antenna description.

        Parameters
        ----------
        source : str, Path, dict, or None
            Either a path to a TOML ALOHA antenna file,
            a dictionary containing antenna parameters,
            or None to create an empty antenna.
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
            antenna_data = tomllib.load(fp)
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
        _ant.antenna = antenna_dict
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

    def __str__(self):
        return f"Antenna(name={self.antenna['name']}, frequency={self.antenna['frequency']})"
