import pathlib
import unittest

from aloha.utils import convert_matlab_scenario, load_m_file, load_mat_file, parse_matlab_array

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

    def test_parse_matlab_scenario_scripts(self):
        """Test parsing antenna_8_active_waveguides.m file."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "antenna_8_active_waveguides.m"
        scenario_dict, _ = load_m_file(m_file)

        # Check that scenario_dict contains expected keys
        self.assertIn("antenna_lh", scenario_dict)
        self.assertIn("modules", scenario_dict)
        self.assertIn("waveguides", scenario_dict)

        # Check antenna_lh structure
        antenna_lh = scenario_dict["antenna_lh"]
        self.assertIn("name", antenna_lh)
        self.assertIn("frequency", antenna_lh)
        self.assertIn("power", antenna_lh)
        self.assertEqual(antenna_lh["name"], "Elementary antenna of 8 active waveguides")
        self.assertEqual(antenna_lh["frequency"], 3.7e9)

        # Check modules structure
        modules = scenario_dict["modules"]
        self.assertIn("nma_theta", modules)
        self.assertIn("nma_phi", modules)
        self.assertEqual(modules["nma_theta"], 1)
        self.assertEqual(modules["nma_phi"], 8)

        # Check waveguides structure
        waveguides = scenario_dict["waveguides"]
        self.assertIn("nwm_theta", waveguides)
        self.assertIn("nwm_phi", waveguides)

    def test_convert_matlab_scenario_basic(self):
        """Test basic conversion of MATLAB scenario to TOML format."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        matlab_scenario, _ = load_m_file(m_file)
        converted = convert_matlab_scenario(matlab_scenario)

        # Check that converted scenario has expected top-level keys
        self.assertIn("antenna", converted)
        self.assertIn("plasma", converted)
        self.assertIn("options", converted)

    def test_convert_matlab_scenario_antenna(self):
        """Test conversion of antenna section from MATLAB to TOML format."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        matlab_scenario, _ = load_m_file(m_file)
        converted = convert_matlab_scenario(matlab_scenario)

        antenna = converted["antenna"]

        # Check antenna file mapping
        self.assertIn("file", antenna)
        self.assertEqual(antenna["file"], "simple_antenna.toml")

        # Check excitation parameters
        self.assertIn("excitation", antenna)
        excitation = antenna["excitation"]

        # Check frequency
        self.assertIn("f", excitation)
        self.assertEqual(excitation["f"], 3700000000.0)

        # Check experimental flag
        self.assertIn("experimental", excitation)
        self.assertTrue(excitation["experimental"])

        # Check port
        self.assertIn("port", excitation)
        self.assertEqual(excitation["port"], "Q6B")

        # Check pulse
        self.assertIn("pulse", excitation)
        self.assertEqual(excitation["pulse"], 45155)

        # Check avg_times
        self.assertIn("avg_times", excitation)
        self.assertEqual(excitation["avg_times"], [7.0, 8.0])

    def test_convert_matlab_scenario_plasma(self):
        """Test conversion of plasma section from MATLAB to TOML format."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        matlab_scenario, _ = load_m_file(m_file)
        converted = convert_matlab_scenario(matlab_scenario)

        plasma = converted["plasma"]

        # Check solver
        self.assertIn("solver", plasma)
        self.assertEqual(plasma["solver"], "spectral_1D")

        # Check spectral_1D section
        self.assertIn("spectral_1D", plasma)
        spectral_1d = plasma["spectral_1D"]

        # Check profile
        self.assertIn("profile", spectral_1d)
        self.assertEqual(spectral_1d["profile"], "bilinear")

        # Check bilinear section
        self.assertIn("bilinear", spectral_1d)
        bilinear = spectral_1d["bilinear"]

        # Check bilinear parameters
        self.assertIn("ne0", bilinear)
        self.assertEqual(bilinear["ne0"], 5e17)

        self.assertIn("lambda_n", bilinear)
        self.assertEqual(bilinear["lambda_n"], [0.002, 0.02])

        self.assertIn("d_couche", bilinear)
        self.assertEqual(bilinear["d_couche"], 0.002)

        self.assertIn("d_vide", bilinear)
        self.assertEqual(bilinear["d_vide"], 0.0)

        self.assertIn("B0", bilinear)
        self.assertEqual(bilinear["B0"], 2.95)

        # Check spectral domain parameters
        self.assertIn("nz_min", spectral_1d)
        self.assertEqual(spectral_1d["nz_min"], -20)

        self.assertIn("nz_max", spectral_1d)
        self.assertEqual(spectral_1d["nz_max"], 20)

        self.assertIn("dnz", spectral_1d)
        self.assertEqual(spectral_1d["dnz"], 0.01)

        # Check spatial domain parameters
        self.assertIn("z_min", spectral_1d)
        self.assertEqual(spectral_1d["z_min"], -0.015)

        self.assertIn("z_max", spectral_1d)
        self.assertEqual(spectral_1d["z_max"], 0.075)

    def test_convert_matlab_scenario_options(self):
        """Test conversion of options section from MATLAB to TOML format."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        matlab_scenario, _ = load_m_file(m_file)
        converted = convert_matlab_scenario(matlab_scenario)

        options = converted["options"]

        # Check debug flag
        self.assertIn("debug", options)
        self.assertTrue(options["debug"])


if __name__ == "__main__":
    unittest.main()
