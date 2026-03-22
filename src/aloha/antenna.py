import os
import tomllib
from collections import defaultdict

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

    def __str__(self):
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

        total_nb_modules = antenna["global"]["nma_phi"] * antenna["global"]["nma_theta"]
        # excitations
        if len(antenna["excitation"]["magnitudes"]) != total_nb_modules:
            raise ValueError(
                f"Number of excitation magnitudes ({len(antenna['excitation']['magnitudes'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )
        if len(antenna["excitation"]["phases"]) != total_nb_modules:
            raise ValueError(
                f"Number of excitation phases ({len(antenna['excitation']['phases'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )

        # module indices
        if len(antenna["global"]["ima"]) != total_nb_modules:
            raise ValueError(
                f"Number of module indices ({len(antenna['global']['ima'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )

        # S-parameter
        if len(antenna["sparameters"]["filenames"]) != total_nb_modules:
            raise ValueError(
                f"Number of S-parameter filenames ({len(antenna['sparameters']['filenames'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )
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
        state : bool
            True is the antenna description looks valid.
        """
        try:
            self.validate_description(self.antenna)
            return True
        except ValueError:
            return False
