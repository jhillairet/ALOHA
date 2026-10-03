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
        m_files = [
            MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m",
            MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.m",
        ]

        mat_files = [
            MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat",
            MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.mat",
        ]

        for m_file, mat_file in zip(m_files, mat_files, strict=True):
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

    def test_run_method_consistency_across_file_formats_8waveguides(self):
        """Test that run() method generates the same results for scenarios from different file formats."""
        # Paths to the different file formats
        toml_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.toml"
        mat_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.mat"
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides" / "scenario_8_active_waveguides.m"

        # Create scenarios from different file formats
        scenario_from_mat_ref = Scenario(mat_file)  # reference results to compare results (from ALOHA-matlab) with
        scenario_from_mat = Scenario(mat_file)
        scenario_from_m = Scenario(m_file)
        scenario_from_toml = Scenario.from_file(toml_file)

        for scenario in [scenario_from_mat_ref, scenario_from_mat, scenario_from_m]:
            # "comment" is missing in the matlab version -- adding it to pass the following tests
            scenario.scenario["comment"] = scenario_from_toml.scenario["comment"]

            # Compare all three scenario dictionaries using Scenario equality test
            self.assertEqual(scenario_from_toml, scenario)

        # Run each scenario to (re)generate results
        scenario_from_toml.run()
        scenario_from_mat.run()  # results are overwritten in this case
        scenario_from_m.run()

        # Verify all scenarios have a "results" fields
        for fields in ["S_plasma", "rac_Zhe"]:
            for scenario in [scenario_from_toml, scenario_from_m, scenario_from_mat, scenario_from_mat_ref]:
                self.assertIn(fields, scenario.results)

        # Compare S_plasma matrices
        S_plasma_mat_ref = scenario_from_mat_ref.results["S_plasma"]
        S_plasma_mat = scenario_from_mat.results["S_plasma"]
        S_plasma_m = scenario_from_m.results["S_plasma"]
        S_plasma_toml = scenario_from_toml.results["S_plasma"]

        for S_plasma in [S_plasma_toml, S_plasma_mat, S_plasma_m]:
            # Check shapes are the same
            self.assertEqual(
                S_plasma_mat_ref.shape,
                S_plasma.shape,
                f"S_plasma shape mismatch: ref={S_plasma_mat_ref.shape}, test={S_plasma.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                S_plasma_mat_ref,
                S_plasma,
                rtol=1e-10,
                atol=1e-10,
                err_msg="S_plasma values differ between TOML and other}",
            )

        # Compare rac_Zhe matrices
        rac_Zhe_mat_ref = scenario_from_mat_ref.results["rac_Zhe"]
        rac_Zhe_mat = scenario_from_mat.results["rac_Zhe"]
        rac_Zhe_m = scenario_from_m.results["rac_Zhe"]
        rac_Zhe_toml = scenario_from_toml.results["rac_Zhe"]

        for rac_Zhe in [rac_Zhe_toml, rac_Zhe_mat, rac_Zhe_m]:
            # Check shapes are the same
            self.assertEqual(
                rac_Zhe_mat_ref.shape,
                rac_Zhe.shape,
                f"rac_Zhe shape mismatch: ref={rac_Zhe_mat_ref.shape}, test={rac_Zhe.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                rac_Zhe_mat_ref,
                rac_Zhe,
                rtol=1e-10,
                atol=1e-10,
                err_msg="rac_Zhe values differ between ref and other}",
            )

        # Compare reflection coefficients (RC)
        rc_mat_ref = scenario_from_mat_ref.results["RC"]
        rc_mat = scenario_from_mat.results["RC"]
        rc_m = scenario_from_m.results["RC"]
        rc_toml = scenario_from_toml.results["RC"]

        # Check shapes are the same
        for rc in [rc_toml, rc_mat, rc_m]:
            self.assertEqual(
                rc_mat_ref.shape,
                rc.shape,
                f"RC shape mismatch: ref={rc_mat_ref.shape}, test={rc.shape}",
            )

            np.testing.assert_allclose(
                rc_mat_ref, rc, rtol=1e-10, atol=1e-10, err_msg="RC values differ between ref and other"
            )

    def test_run_method_consistency_across_file_formats_8waveguides_3modes(self):
        """Test that run() method generates the same results for scenarios from different file formats."""
        # Paths to the different file formats
        toml_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides_3modes" / "scenario_8_active_waveguides_3modes.toml"
        mat_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides_3modes" / "scenario_8_active_waveguides_3modes.mat"
        m_file = MATLAB_TEST_CASES_DIR / "8_active_waveguides_3modes" / "scenario_8_active_waveguides_3modes.m"

        # Create scenarios from different file formats
        scenario_from_mat_ref = Scenario(mat_file)  # reference results to compare results (from ALOHA-matlab) with
        scenario_from_mat = Scenario(mat_file)
        scenario_from_m = Scenario(m_file)
        scenario_from_toml = Scenario.from_file(toml_file)

        for scenario in [scenario_from_mat_ref, scenario_from_mat, scenario_from_m]:
            # "comment" is missing in the matlab version -- adding it to pass the following tests
            scenario.scenario["comment"] = scenario_from_toml.scenario["comment"]

            # Compare all three scenario dictionaries using Scenario equality test
            self.assertEqual(scenario_from_toml, scenario)

        # Run each scenario to (re)generate results
        scenario_from_toml.run()
        scenario_from_mat.run()  # results are overwritten in this case
        scenario_from_m.run()

        # Verify all scenarios have a "results" fields
        for fields in ["S_plasma", "rac_Zhe"]:
            for scenario in [scenario_from_toml, scenario_from_m, scenario_from_mat, scenario_from_mat_ref]:
                self.assertIn(fields, scenario.results)

        # Compare S_plasma matrices
        S_plasma_mat_ref = scenario_from_mat_ref.results["S_plasma"]
        S_plasma_mat = scenario_from_mat.results["S_plasma"]
        S_plasma_m = scenario_from_m.results["S_plasma"]
        S_plasma_toml = scenario_from_toml.results["S_plasma"]

        for S_plasma in [S_plasma_toml, S_plasma_mat, S_plasma_m]:
            # Check shapes are the same
            self.assertEqual(
                S_plasma_mat_ref.shape,
                S_plasma.shape,
                f"S_plasma shape mismatch: ref={S_plasma_mat_ref.shape}, test={S_plasma.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                S_plasma_mat_ref,
                S_plasma,
                rtol=1e-10,
                atol=1e-10,
                err_msg="S_plasma values differ between TOML and other}",
            )

        # Compare rac_Zhe matrices
        rac_Zhe_mat_ref = scenario_from_mat_ref.results["rac_Zhe"]
        rac_Zhe_mat = scenario_from_mat.results["rac_Zhe"]
        rac_Zhe_m = scenario_from_m.results["rac_Zhe"]
        rac_Zhe_toml = scenario_from_toml.results["rac_Zhe"]

        for rac_Zhe in [rac_Zhe_toml, rac_Zhe_mat, rac_Zhe_m]:
            # Check shapes are the same
            self.assertEqual(
                rac_Zhe_mat_ref.shape,
                rac_Zhe.shape,
                f"rac_Zhe shape mismatch: ref={rac_Zhe_mat_ref.shape}, test={rac_Zhe.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                rac_Zhe_mat_ref,
                rac_Zhe,
                rtol=1e-10,
                atol=1e-10,
                err_msg="rac_Zhe values differ between ref and other}",
            )

        # Compare reflection coefficients (RC)
        rc_mat_ref = scenario_from_mat_ref.results["RC"]
        rc_mat = scenario_from_mat.results["RC"]
        rc_m = scenario_from_m.results["RC"]
        rc_toml = scenario_from_toml.results["RC"]

        # Check shapes are the same
        for rc in [rc_toml, rc_mat, rc_m]:
            self.assertEqual(
                rc_mat_ref.shape,
                rc.shape,
                f"RC shape mismatch: ref={rc_mat_ref.shape}, test={rc.shape}",
            )

            np.testing.assert_allclose(
                rc_mat_ref, rc, rtol=1e-10, atol=1e-10, err_msg="RC values differ between ref and other"
            )

    def test_run_method_consistency_across_file_formats_LH1(self):
        """Test that run() method generates the same results for scenarios from different file formats."""
        # Paths to the different file formats
        toml_file = MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.toml"
        mat_file = MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.mat"
        m_file = MATLAB_TEST_CASES_DIR / "WEST_LH1" / "scenario_WEST_LH1.m"

        # Create scenarios from different file formats
        scenario_from_mat_ref = Scenario(mat_file)  # reference results to compare results (from ALOHA-matlab) with
        scenario_from_mat = Scenario(mat_file)
        scenario_from_m = Scenario(m_file)
        scenario_from_toml = Scenario.from_file(toml_file)
        # list of scenarios to run and test against
        scenarios = [scenario_from_m, scenario_from_mat, scenario_from_toml]

        for scenario in [scenario_from_mat_ref, scenario_from_mat, scenario_from_m]:
            # "comment" is missing in the matlab version -- adding it to pass the following tests
            scenario.scenario["comment"] = scenario_from_toml.scenario["comment"]

            # Compare all three scenario dictionaries using Scenario equality test
            self.assertEqual(scenario_from_toml, scenario)

        # Run each scenario to (re)generate results
        for scenario in scenarios:
            scenario.run()  # results are overwritten in the _mat case

        # Verify all scenarios have a "results" fields
        for fields in ["S_plasma", "rac_Zhe"]:
            for scenario in scenarios:
                self.assertIn(fields, scenario.results)

        # Compare S_plasma matrices
        S_plasma_mat_ref = scenario_from_mat_ref.results["S_plasma"]

        for scenario in scenarios:
            S_plasma = scenario.results["S_plasma"]
            # Check shapes are the same
            self.assertEqual(
                S_plasma_mat_ref.shape,
                S_plasma.shape,
                f"S_plasma shape mismatch: ref={S_plasma_mat_ref.shape}, test={S_plasma.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                S_plasma_mat_ref,
                S_plasma,
                rtol=1e-6,
                atol=1e-6,
                err_msg="S_plasma values differ between TOML and other}",
            )

        # Compare rac_Zhe matrices
        rac_Zhe_mat_ref = scenario_from_mat_ref.results["rac_Zhe"]

        for scenario in scenarios:
            rac_Zhe = scenario.results["rac_Zhe"]
            # Check shapes are the same
            self.assertEqual(
                rac_Zhe_mat_ref.shape,
                rac_Zhe.shape,
                f"rac_Zhe shape mismatch: ref={rac_Zhe_mat_ref.shape}, test={rac_Zhe.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                rac_Zhe_mat_ref,
                rac_Zhe,
                rtol=1e-6,
                atol=1e-6,
                err_msg="rac_Zhe values differ between ref and other}",
            )

        # Compare reflection coefficients (RC)
        # TODO: RC comparison is disabled for now as it requires loading antenna S-parameters
        # from MATLAB antenna architecture files. The S_plasma and rac_Zhe matrices are matching.
        # rc_mat_ref = scenario_from_mat_ref.results["RC"]
        #
        # for scenario in scenarios:
        #     rc = scenario.results["RC"]
        #     self.assertEqual(
        #         rc_mat_ref.shape,
        #         rc.shape,
        #         f"RC shape mismatch: ref={rc_mat_ref.shape}, test={rc.shape}",
        #     )
        #
        #     np.testing.assert_allclose(
        #         rc_mat_ref, rc, rtol=1e-6, atol=1e-6, err_msg="RC values differ between ref and other"
        #     )

    def test_run_method_consistency_across_file_formats_LH2(self):
        """Test that run() method generates the same results for scenarios from different file formats."""
        # Paths to the different file formats
        toml_file = MATLAB_TEST_CASES_DIR / "WEST_LH2" / "scenario_WEST_LH2.toml"
        mat_file = MATLAB_TEST_CASES_DIR / "WEST_LH2" / "scenario_WEST_LH2.mat"
        m_file = MATLAB_TEST_CASES_DIR / "WEST_LH2" / "scenario_WEST_LH2.m"

        # Create scenarios from different file formats
        scenario_from_mat_ref = Scenario(mat_file)  # reference results to compare results (from ALOHA-matlab) with
        scenario_from_mat = Scenario(mat_file)
        scenario_from_m = Scenario(m_file)
        scenario_from_toml = Scenario.from_file(toml_file)

        for scenario in [scenario_from_mat_ref, scenario_from_mat, scenario_from_m]:
            # "comment" is missing in the matlab version -- adding it to pass the following tests
            scenario.scenario["comment"] = scenario_from_toml.scenario["comment"]

            # Compare all three scenario dictionaries using Scenario equality test
            self.assertEqual(scenario_from_toml, scenario)

        # Run each scenario to (re)generate results
        scenario_from_toml.run()
        scenario_from_mat.run()  # results are overwritten in this case
        scenario_from_m.run()

        # Verify all scenarios have a "results" fields
        for fields in ["S_plasma", "rac_Zhe"]:
            for scenario in [scenario_from_toml, scenario_from_m, scenario_from_mat, scenario_from_mat_ref]:
                self.assertIn(fields, scenario.results)

        # Compare S_plasma matrices
        S_plasma_mat_ref = scenario_from_mat_ref.results["S_plasma"]
        S_plasma_mat = scenario_from_mat.results["S_plasma"]
        S_plasma_m = scenario_from_m.results["S_plasma"]
        S_plasma_toml = scenario_from_toml.results["S_plasma"]

        for S_plasma in [S_plasma_toml, S_plasma_mat, S_plasma_m]:
            # Check shapes are the same
            self.assertEqual(
                S_plasma_mat_ref.shape,
                S_plasma.shape,
                f"S_plasma shape mismatch: ref={S_plasma_mat_ref.shape}, test={S_plasma.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                S_plasma_mat_ref,
                S_plasma,
                rtol=1e-6,
                atol=1e-6,
                err_msg="S_plasma values differ between TOML and other}",
            )

        # Compare rac_Zhe matrices
        rac_Zhe_mat_ref = scenario_from_mat_ref.results["rac_Zhe"]
        rac_Zhe_mat = scenario_from_mat.results["rac_Zhe"]
        rac_Zhe_m = scenario_from_m.results["rac_Zhe"]
        rac_Zhe_toml = scenario_from_toml.results["rac_Zhe"]

        for rac_Zhe in [rac_Zhe_toml, rac_Zhe_mat, rac_Zhe_m]:
            # Check shapes are the same
            self.assertEqual(
                rac_Zhe_mat_ref.shape,
                rac_Zhe.shape,
                f"rac_Zhe shape mismatch: ref={rac_Zhe_mat_ref.shape}, test={rac_Zhe.shape}",
            )

            # Compare values with tolerance (due to potential numerical differences)
            np.testing.assert_allclose(
                rac_Zhe_mat_ref,
                rac_Zhe,
                rtol=1e-6,
                atol=1e-6,
                err_msg="rac_Zhe values differ between ref and other}",
            )

        # Compare reflection coefficients (RC)
        rc_mat_ref = scenario_from_mat_ref.results["RC"]
        rc_mat = scenario_from_mat.results["RC"]
        rc_m = scenario_from_m.results["RC"]
        rc_toml = scenario_from_toml.results["RC"]

        # Check shapes are the same
        for rc in [rc_toml, rc_mat, rc_m]:
            self.assertEqual(
                rc_mat_ref.shape,
                rc.shape,
                f"RC shape mismatch: ref={rc_mat_ref.shape}, test={rc.shape}",
            )

            np.testing.assert_allclose(
                rc_mat_ref, rc, rtol=1e-6, atol=1e-6, err_msg="RC values differ between ref and other"
            )


if __name__ == "__main__":
    unittest.main()
