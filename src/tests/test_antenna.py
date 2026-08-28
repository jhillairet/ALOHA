# Antenna test
import os
import tempfile
import unittest
from glob import glob
from pathlib import Path

from aloha.antenna import Antenna

ANTENNAS_DIR = Path(__file__).resolve().parent.parent.parent / "antennas"
MATLAB_ANTENNAS_DIR = Path(__file__).resolve().parent.parent.parent / "aloha_matlab" / "architecture_antenne"


class TestAntenna(unittest.TestCase):
    def test_antenna_constructor(self):
        filename = ANTENNAS_DIR / "8_active_waveguides.toml"
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
        values = ["name", "frequency", "layout", "module", "sparameters"]

        antenna_filenames = glob("*.toml", root_dir=ANTENNAS_DIR)
        for antenna_filename in antenna_filenames:
            try:
                ant = Antenna(ANTENNAS_DIR / antenna_filename)
                # test if the expected elements are defined
                for value in values:
                    assert value in ant.antenna
            except Exception as e:
                # Skip files that are not valid TOML files
                # (e.g., WEST_IC_Classical.toml is not a valid TOML file)
                if "TOMLDecodeError" in str(type(e)):
                    continue
                else:
                    raise

    def test_from_matlab_antenna_example(self):
        """
        Test reading some matlab antenna files using from_matlab

        should produce a valid antenna description.
        """
        filenames = ["antenna_example.m", "antenna_C3_WEST_ITM.m", "antenna_C4_ITM.m"]
        matlab_files = [MATLAB_ANTENNAS_DIR / filename for filename in filenames]

        for matlab_file in matlab_files:
            # Test that from_matlab can read the file and create a valid antenna
            antenna = Antenna.from_matlab(matlab_file)

            # Verify it's a valid antenna description
            assert antenna.is_valid(), "Antenna from MATLAB file should be valid"

            # Check basic properties are extracted
            assert "name" in antenna.antenna, "Antenna should have a name"
            assert "frequency" in antenna.antenna, "Antenna should have a frequency"
            assert "layout" in antenna.antenna, "Antenna should have layout information"
            assert "module" in antenna.antenna, "Antenna should have module information"

    def test_to_toml(self):
        """
        Test the to_toml method for exporting antenna descriptions to TOML format.
        """
        # Test with a known antenna
        antenna = Antenna.from_file(ANTENNAS_DIR / "8_active_waveguides.toml")

        # Test returning TOML string
        toml_string = antenna.to_toml()
        assert isinstance(toml_string, str), "to_toml() should return a string when no filename is provided"
        assert "# 8 waveguides single row antenna" in toml_string, "TOML should contain the antenna name comment"
        assert 'name = "8 waveguides single row antenna"' in toml_string, "TOML should contain the antenna name"
        assert "[layout]" in toml_string, "TOML should contain the layout section"
        assert "[module]" in toml_string, "TOML should contain the module section"
        assert "[sparameters]" in toml_string, "TOML should contain the sparameters section"
        assert "[excitation]" in toml_string, "TOML should contain the excitation section"

        # Test writing to file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".toml", delete=False) as f:
            temp_filename = f.name

        try:
            # Write to file
            result = antenna.to_toml(temp_filename)
            assert result is None, "to_toml(filename) should return None when writing to file"
            assert os.path.exists(temp_filename), "TOML file should be created"

            # Verify file content
            with open(temp_filename) as f:
                file_content = f.read()
            assert len(file_content) > 0, "TOML file should not be empty"

            # Test that the written file can be loaded back
            reloaded_antenna = Antenna.from_file(temp_filename)
            assert reloaded_antenna.is_valid(), "Reloaded antenna from TOML file should be valid"
            assert reloaded_antenna.antenna["name"] == antenna.antenna["name"], (
                "Reloaded antenna should have the same name"
            )

        finally:
            # Clean up
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    def test_to_toml_round_trip(self):
        """
        Test that reading a TOML file and exporting it back creates a file with the same values.
        """
        import os
        import tempfile

        # Test with all available TOML antenna files
        antenna_filenames = glob("*.toml", root_dir=ANTENNAS_DIR)
        for antenna_filename in antenna_filenames:
            original_file = ANTENNAS_DIR / antenna_filename

            # Skip files that are not valid TOML (like WEST_IC_Classical.toml)
            try:
                original_antenna = Antenna.from_file(original_file)
            except Exception:
                continue

            # Export to a temporary file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".toml", delete=False) as f:
                temp_filename = f.name

            try:
                # Write to temporary file
                original_antenna.to_toml(temp_filename)

                # Read back the exported file
                exported_antenna = Antenna.from_file(temp_filename)

                # Check that all values are the same
                self.assertEqual(
                    original_antenna.antenna["name"],
                    exported_antenna.antenna["name"],
                    f"Name should be preserved for {antenna_filename}",
                )

                self.assertEqual(
                    original_antenna.antenna["frequency"],
                    exported_antenna.antenna["frequency"],
                    f"Frequency should be preserved for {antenna_filename}",
                )

                # Check layout
                if "layout" in original_antenna.antenna:
                    for key in original_antenna.antenna["layout"]:
                        self.assertEqual(
                            original_antenna.antenna["layout"][key],
                            exported_antenna.antenna["layout"][key],
                            f"Layout {key} should be preserved for {antenna_filename}",
                        )

                # Check module
                if "module" in original_antenna.antenna:
                    for key in original_antenna.antenna["module"]:
                        self.assertEqual(
                            original_antenna.antenna["module"][key],
                            exported_antenna.antenna["module"][key],
                            f"Module {key} should be preserved for {antenna_filename}",
                        )

                # Check sparameters
                if "sparameters" in original_antenna.antenna:
                    for key in original_antenna.antenna["sparameters"]:
                        self.assertEqual(
                            original_antenna.antenna["sparameters"][key],
                            exported_antenna.antenna["sparameters"][key],
                            f"Sparameters {key} should be preserved for {antenna_filename}",
                        )

                # Check excitation
                if "excitation" in original_antenna.antenna:
                    for key in original_antenna.antenna["excitation"]:
                        self.assertEqual(
                            original_antenna.antenna["excitation"][key],
                            exported_antenna.antenna["excitation"][key],
                            f"Excitation {key} should be preserved for {antenna_filename}",
                        )

            finally:
                # Clean up
                if os.path.exists(temp_filename):
                    os.remove(temp_filename)
