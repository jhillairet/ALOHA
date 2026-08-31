import pathlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

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

    def test_scenario_constructor_matlab_files(self):
        """Test that the constructor can handle .m and .mat files directly."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"
        mat_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat"

        # Test .m file
        scenario_from_m = Scenario(m_file)
        self.assertIsInstance(scenario_from_m, Scenario)
        self.assertTrue(scenario_from_m.scenario)

        # Test .mat file
        scenario_from_mat = Scenario(mat_file)
        self.assertIsInstance(scenario_from_mat, Scenario)
        self.assertTrue(scenario_from_mat.scenario)

    def test_matlab_to_toml_roundtrip(self):
        """Test that reading MATLAB scenario, saving to TOML, and reading back produces the same scenario dictionary."""
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"

        # Step 1: Read the MATLAB scenario using constructor
        scenario_obj = Scenario(m_file)

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

        # Load both files using constructor
        scenario_from_m = Scenario(m_file)
        scenario_from_mat = Scenario(mat_file)

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

    def test_run_method_consistency_across_file_formats(self):
        """Test that run() method generates the same results for scenarios from different file formats."""
        # Paths to the different file formats
        toml_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.toml"
        mat_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat"
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"

        # Check if the Fortran binary exists before running tests
        from aloha.plasma import get_binary_path

        binary_path = get_binary_path(6, "glnxa64")
        if not binary_path.exists():
            self.skipTest(f"Fortran binary not found: {binary_path}")

        # Create scenarios from different file formats
        scenario_from_toml = Scenario.from_file(toml_file)
        scenario_from_mat = Scenario(mat_file)
        scenario_from_m = Scenario(m_file)

        # check that the scenario dictionnary is the same between formats
        # "comment" is missing in the matlab version -- adding it to pass the test
        scenario_from_m.scenario["comment"] = scenario_from_toml.scenario["comment"]
        scenario_from_mat.scenario["comment"] = scenario_from_toml.scenario["comment"]

        # Compare all three scenario dictionaries using Scenario equality test
        self.assertEqual(scenario_from_toml, scenario_from_m)
        self.assertEqual(scenario_from_toml, scenario_from_mat)

        # Run each scenario to generate results
        try:
            scenario_from_toml.run()
            scenario_from_mat.run()
            scenario_from_m.run()
        except RuntimeError as e:
            if "Binary execution failed" in str(e):
                self.skipTest(f"Fortran binary execution failed: {e}")
            else:
                raise

        # Verify all scenarios have results
        self.assertIn("S_plasma", scenario_from_toml.results)
        self.assertIn("rac_Zhe", scenario_from_toml.results)
        self.assertIn("S_plasma", scenario_from_mat.results)
        self.assertIn("rac_Zhe", scenario_from_mat.results)
        self.assertIn("S_plasma", scenario_from_m.results)
        self.assertIn("rac_Zhe", scenario_from_m.results)

        # Compare S_plasma matrices
        S_plasma_toml = scenario_from_toml.results["S_plasma"]
        S_plasma_mat = scenario_from_mat.results["S_plasma"]
        S_plasma_m = scenario_from_m.results["S_plasma"]

        # Check shapes are the same
        self.assertEqual(
            S_plasma_toml.shape,
            S_plasma_mat.shape,
            f"S_plasma shape mismatch: TOML={S_plasma_toml.shape}, MAT={S_plasma_mat.shape}",
        )
        self.assertEqual(
            S_plasma_toml.shape,
            S_plasma_m.shape,
            f"S_plasma shape mismatch: TOML={S_plasma_toml.shape}, M={S_plasma_m.shape}",
        )

        # Compare rac_Zhe matrices
        rac_Zhe_toml = scenario_from_toml.results["rac_Zhe"]
        rac_Zhe_mat = scenario_from_mat.results["rac_Zhe"]
        rac_Zhe_m = scenario_from_m.results["rac_Zhe"]

        # Check shapes are the same
        self.assertEqual(
            rac_Zhe_toml.shape,
            rac_Zhe_mat.shape,
            f"rac_Zhe shape mismatch: TOML={rac_Zhe_toml.shape}, MAT={rac_Zhe_mat.shape}",
        )
        self.assertEqual(
            rac_Zhe_toml.shape,
            rac_Zhe_m.shape,
            f"rac_Zhe shape mismatch: TOML={rac_Zhe_toml.shape}, M={rac_Zhe_m.shape}",
        )

        # Compare values with tolerance (due to potential numerical differences)
        np.testing.assert_allclose(
            S_plasma_toml,
            S_plasma_mat,
            rtol=1e-10,
            atol=1e-10,
            err_msg="S_plasma values differ between TOML and MAT files",
        )
        np.testing.assert_allclose(
            S_plasma_toml, S_plasma_m, rtol=1e-10, atol=1e-10, err_msg="S_plasma values differ between TOML and M files"
        )

        np.testing.assert_allclose(
            rac_Zhe_toml,
            rac_Zhe_mat,
            rtol=1e-10,
            atol=1e-10,
            err_msg="rac_Zhe values differ between TOML and MAT files",
        )
        np.testing.assert_allclose(
            rac_Zhe_toml, rac_Zhe_m, rtol=1e-10, atol=1e-10, err_msg="rac_Zhe values differ between TOML and M files"
        )


if __name__ == "__main__":
    unittest.main()
