import pathlib
import unittest

from aloha.utils import load_m_file, load_mat_file, parse_matlab_array

MATLAB_TEST_CASES_DIR = pathlib.Path(__file__).parent / "matlab_reference_cases"


class TestUtils(unittest.TestCase):
    def test_parse_matlab_array_space_separated(self):
        """Test parsing space-separated MATLAB array string."""
        result = parse_matlab_array("1 2 3")
        self.assertEqual(result, ["1", "2", "3"])

    def test_parse_matlab_array_comma_separated(self):
        """Test parsing comma-separated MATLAB array string."""
        result = parse_matlab_array("1, 2, 3")
        self.assertEqual(result, ["1", "2", "3"])

    def test_parse_matlab_array_comma_no_space(self):
        """Test parsing comma-separated MATLAB array string without spaces."""
        result = parse_matlab_array("1,2,3")
        self.assertEqual(result, ["1", "2", "3"])

    def test_parse_matlab_array_empty(self):
        """Test parsing empty string."""
        result = parse_matlab_array("")
        self.assertEqual(result, [])

    def test_parse_matlab_array_single_element(self):
        """Test parsing single element."""
        result = parse_matlab_array("42")
        self.assertEqual(result, ["42"])

    def test_parse_matlab_array_mixed_types(self):
        """Test parsing array with mixed types (as strings)."""
        result = parse_matlab_array("1 a 3.14")
        self.assertEqual(result, ["1", "a", "3.14"])

    def test_parse_matlab_array_whitespace_handling(self):
        """Test parsing with leading/trailing whitespace."""
        result = parse_matlab_array("  1 2 3  ")
        self.assertEqual(result, ["1", "2", "3"])

    def test_load_mat_files(self):
        """Test loading a .mat file from 8_active_waveguides case."""
        mat_files = [
            MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat",
            MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.mat",
            MATLAB_TEST_CASES_DIR / "WEST_LH2" / "scenario_WEST_LH2.mat",
        ]

        for mat_file in mat_files:
            scenario_dict, results_dict = load_mat_file(mat_file)

            # Check that scenario_dict contains expected keys
            self.assertIn("antenna", scenario_dict)
            self.assertIn("antenna_lh", scenario_dict)
            self.assertIn("options", scenario_dict)
            self.assertIn("plasma", scenario_dict)

            # Check that results_dict contains expected keys
            self.assertIn("Pin", results_dict)


if __name__ == "__main__":
    unittest.main()
