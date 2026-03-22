# Antenna test
import unittest
from glob import glob
from pathlib import Path

from aloha.antenna import Antenna

ANTENNAS_DIR = Path(__file__).resolve().parent.parent.parent / "antennas"


class TestAntenna(unittest.TestCase):
    def test_antenna_constructor(self):
        filename = ANTENNAS_DIR / "simple_antenna.toml"
        # from a file
        Antenna(filename)
        Antenna.from_file(filename)
        # from a dict
        ant_dict = Antenna.load(filename)
        Antenna(ant_dict)
        Antenna.from_dict(ant_dict)

    def test_validate_antenna_files(self):
        """
        Validation of the TOML schema of all the pre-defined antennas.
        """
        # Check if some fields are in the antenna description
        # Some logic (number of modules, etc) is tested when creating an Antenna object.
        values = ["name", "frequency", "global", "module", "sparameters"]

        antenna_filenames = glob("*.toml", root_dir=ANTENNAS_DIR)
        for antenna_filename in antenna_filenames:
            ant = Antenna(ANTENNAS_DIR / antenna_filename)
            # test if the expected elements are defined
            for value in values:
                assert value in ant.antenna
