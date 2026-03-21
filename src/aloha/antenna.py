import os
import tomllib

import numpy as np

from .constants import pi
from .waveguide import Waveguide


class Antenna:
    """
    ALOHA antenna description.
    """

    def __init__(self, filename: str | os.PathLike = None):
        """
        ALOHA antenna description.

        Parameters
        ----------
        filename : str | os.PathLike
            Path to a TOML ALOHA antenna file.

        """
        if filename:
            self.antenna = self.load(filename)

    @classmethod
    def load(cls, filename: str | os.PathLike) -> dict:
        """
        Read an ALOHA antenna description (TOML file format).

        Parameters
        ----------
        filename : str
            path to a TOML file.

        Returns
        -------
        antenna : dict
            ALOHA antenna description.

        """
        with open(filename, "rb") as fp:
            return tomllib.load(fp)
