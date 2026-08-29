import pathlib
import tempfile
import unittest
from pathlib import Path

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from aloha.scenario import Scenario
from aloha.utils import load_m_file

SCENARIOS_DIR = Path(__file__).parent.parent.parent / "scenarios"
TOML_SCENARIO_FILES = list(SCENARIOS_DIR.glob("*.toml"))
MATLAB_TEST_CASES_DIR = pathlib.Path(__file__).parent / "matlab_reference_cases"


class TestScenario(unittest.TestCase):
    def test_scenario_constructor_no_args(self):
        scenario = Scenario()
        self.assertEqual(scenario.scenario, {})

    def test_scenario_constructor_dict(self):
        data = {"options": {"test": True}}
        scenario = Scenario(data)
        self.assertEqual(scenario.scenario, data)

    def test_scenario_constructor_str(self):
        for scenario_file in TOML_SCENARIO_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(str(scenario_file))
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_path(self):
        for scenario_file in TOML_SCENARIO_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(scenario_file)
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_invalid(self):
        with self.assertRaises(ValueError):
            Scenario(123)

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
        converted = Scenario.convert_matlab_scenario(matlab_scenario)

        # Check that converted scenario has expected top-level keys
        self.assertIn("antenna", converted)
        self.assertIn("plasma", converted)
        self.assertIn("options", converted)

    def test_convert_matlab_scenario_antenna(self):
        """Test conversion of antenna section from MATLAB to TOML format."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        matlab_scenario, _ = load_m_file(m_file)
        converted = Scenario.convert_matlab_scenario(matlab_scenario)

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
        converted = Scenario.convert_matlab_scenario(matlab_scenario)

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
        converted = Scenario.convert_matlab_scenario(matlab_scenario)

        options = converted["options"]

        # Check debug flag
        self.assertIn("debug", options)
        self.assertTrue(options["debug"])

    def test_matlab_to_toml_roundtrip(self):
        """Test that reading MATLAB scenario, saving to TOML, and reading back produces the same scenario dictionary."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"

        # Step 1: Read the MATLAB scenario using from_matlab method
        scenario_obj = Scenario.from_matlab(m_file)

        # Step 2: Save the scenario to a temporary TOML file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".toml", delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)

        try:
            # Save to TOML
            scenario_obj.to_toml(tmp_path)

            # Step 3: Read the TOML file using Scenario constructor
            scenario_from_toml = Scenario(tmp_path)

            # Step 4: Assert that the scenario dictionaries are the same
            # Note: The 'comment' field is intentionally not preserved in TOML format
            # as it's converted to actual TOML comments, so we exclude it from comparison
            scenario_dict_no_comment = {k: v for k, v in scenario_obj.scenario.items() if k != "comment"}
            self.assertEqual(scenario_from_toml.scenario, scenario_dict_no_comment)
        finally:
            # Clean up the temporary file
            tmp_path.unlink()

    def test_matlab_files_consistency(self):
        """Test that loading scenario from .m and .mat files produces similar Scenario objects."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        mat_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat"

        # Load both files using from_matlab
        scenario_from_m = Scenario.from_matlab(m_file)
        scenario_from_mat = Scenario.from_matlab(mat_file)

        # Verify both are Scenario objects
        self.assertIsInstance(scenario_from_m, Scenario)
        self.assertIsInstance(scenario_from_mat, Scenario)

        # Verify both have scenario and results attributes
        self.assertTrue(hasattr(scenario_from_m, "scenario"))
        self.assertTrue(hasattr(scenario_from_m, "results"))
        self.assertTrue(hasattr(scenario_from_mat, "scenario"))
        self.assertTrue(hasattr(scenario_from_mat, "results"))

        # Verify both have the expected top-level keys in their scenario
        expected_keys = {"antenna", "plasma", "options"}
        self.assertTrue(expected_keys.issubset(scenario_from_m.scenario.keys()))
        self.assertTrue(expected_keys.issubset(scenario_from_mat.scenario.keys()))

        # Verify that the MAT file has results (since it's a computed scenario)
        self.assertTrue(len(scenario_from_mat.results) > 0)

        # Compare the structure of the scenario dictionaries (excluding comment)
        m_scenario = {k: v for k, v in scenario_from_m.scenario.items() if k != "comment"}
        mat_scenario = {k: v for k, v in scenario_from_mat.scenario.items() if k != "comment"}

        # Check that both have the same top-level keys
        self.assertEqual(set(m_scenario.keys()), set(mat_scenario.keys()))

        # Check that frequency is the same
        self.assertEqual(m_scenario["antenna"]["excitation"]["f"], mat_scenario["antenna"]["excitation"]["f"])

        # Check that power arrays have the same length
        m_power = m_scenario["antenna"]["excitation"]["power"]
        mat_power = mat_scenario["antenna"]["excitation"]["power"]
        self.assertEqual(len(m_power), len(mat_power))

        # Check that phase arrays have the same length
        m_phase = m_scenario["antenna"]["excitation"]["phase"]
        mat_phase = mat_scenario["antenna"]["excitation"]["phase"]
        self.assertEqual(len(m_phase), len(mat_phase))

        # Check that both have the same plasma solver
        self.assertEqual(m_scenario["plasma"]["solver"], mat_scenario["plasma"]["solver"])

        # Check that both have spectral_1D in plasma
        self.assertIn("spectral_1D", m_scenario["plasma"])
        self.assertIn("spectral_1D", mat_scenario["plasma"])

        # Check that both have the same spectral_1D profile
        self.assertEqual(
            m_scenario["plasma"]["spectral_1D"]["profile"], mat_scenario["plasma"]["spectral_1D"]["profile"]
        )


if __name__ == "__main__":
    unittest.main()
