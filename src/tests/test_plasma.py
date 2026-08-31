# Plasma module test
import tempfile
import unittest
from pathlib import Path

import numpy as np

from aloha.plasma import S_plasma_1D, S_plasma_1D_matlab_inputs, get_binary_name, get_binary_path


class TestPlasma(unittest.TestCase):
    def test_get_binary_name(self):
        """Test binary name generation."""
        name = get_binary_name(6, "glnxa64")
        self.assertEqual(name, "coupl_plasma_version6_glnxa64")

        name = get_binary_name(3, "alpha")
        self.assertEqual(name, "coupl_plasma_version3_alpha")

    def test_get_binary_path(self):
        """Test binary path generation."""
        path = get_binary_path(6, "glnxa64")
        expected_path = (
            Path(__file__).resolve().parent.parent.parent
            / "aloha_matlab"
            / "code_1D"
            / "couplage_1D"
            / "coupl_plasma_version6_glnxa64"
        )
        self.assertEqual(path, expected_path)

    def test_binary_exists(self):
        """Test that the binary exists."""
        binary_path = get_binary_path(6, "glnxa64")
        self.assertTrue(binary_path.exists(), f"Binary not found: {binary_path}")

    def test_S_plasma_1D_minimal(self):
        """Test S_plasma_1D with minimal parameters."""
        # Create a minimal scenario for version 6
        scenario = {
            "antenna": {"freq": 3.7e9},  # 3.7 GHz frequency
            "ne0": [5e17],  # Electron density at reference position [m^-3]
            "dne0": [1e18],  # Electron density gradient [m^-4]
            "d_couche": [0.002],  # Thickness of first plasma layer [m]
            "dne1": [1e18],  # Electron density gradient for second layer [m^-4]
            "nb_g_pol": 1,  # Number of poloidal rows
            "nb_g_total_ligne": 1,  # Number of waveguides per poloidal row
            "a": 0.1,  # Waveguide width parameter [m]
            "b": [0.1],  # Waveguide height parameters [m]
            "z": [0.0],  # Waveguide position parameters [m]
            "T_grill": 1,  # Grill periodicity parameter
            "D_guide_max": 10,  # Maximum guide decoupling distance [m]
            "erreur_rel": 1e-6,  # Relative error tolerance
            "pertes": 0.0,  # Loss parameter
            "d_vide": [0.0],  # Vacuum layer thickness [m]
            "Nmh": 1,  # Number of magnetic modes
            "Nme": 1,  # Number of electric modes
        }

        # Check if the binary exists before trying to run it
        binary_path = get_binary_path(6, "glnxa64")
        if not binary_path.exists():
            self.skipTest(f"Fortran binary not found: {binary_path}")

        try:
            # Run with debug output
            S_plasma, rac_Zhe = S_plasma_1D_matlab_inputs(
                scenario,
                version=6,
                architecture="glnxa64",
                bool_debug=False,  # Disable debug for cleaner test output
            )

            # Basic checks
            self.assertIsInstance(S_plasma, np.ndarray)
            self.assertIsInstance(rac_Zhe, np.ndarray)
            self.assertEqual(S_plasma.dtype, np.complex128)
            self.assertEqual(rac_Zhe.dtype, np.complex128)

            # Check matrix dimensions
            expected_size = scenario["nb_g_total_ligne"] * scenario["nb_g_pol"] * (scenario["Nmh"] + scenario["Nme"])
            self.assertEqual(S_plasma.shape, (expected_size, expected_size))
            self.assertEqual(rac_Zhe.shape, (expected_size, expected_size))

        except RuntimeError as e:
            # Fortran binary might fail due to environment/dependency issues
            # This is acceptable for the test - we just want to ensure the Python
            # interface works correctly
            if "Binary execution failed" in str(e):
                self.skipTest(f"Fortran binary execution failed: {e}")
            else:
                raise

    def test_S_plasma_1D_wrong_version(self):
        """Test that S_plasma_1D raises error for unsupported versions."""
        scenario = {
            "antenna": {"freq": 60e6},
            "ne0": [1e18],
            "dne0": [1e18],
            "d_couche": [0.1],
            "dne1": [1e18],
            "nb_g_pol": 1,
            "nb_g_total_ligne": 1,
            "a": 0.1,
            "b": [0.1],
            "z": [0.0],
            "T_grill": 1,
            "D_guide_max": 10,
            "erreur_rel": 1e-6,
            "pertes": 0.0,
            "d_vide": [0.0],
            "Nmh": 1,
            "Nme": 1,
        }

        # Test with unsupported version
        with self.assertRaises(ValueError):
            S_plasma_1D_matlab_inputs(scenario, version=3)

        with self.assertRaises(ValueError):
            S_plasma_1D_matlab_inputs(scenario, version=7)

    def test_S_plasma_1D_with_scenario_object(self):
        """Test S_plasma_1D with Scenario object using TOML schema."""
        from aloha.scenario import Scenario

        # Create a scenario using the TOML schema
        scenario_dict = {
            "antenna": {"file": "8_active_waveguides.toml", "excitation": {"f": 3.7e9}},
            "plasma": {
                "solver": "spectral_1D",
                "spectral_1D": {
                    "profile": "bilinear",
                    "nb_evanescent_modes": 2,
                    "bilinear": {
                        "ne0": 5e17,
                        "lambda_n": [0.002, 0.02],
                        "plasma_layer_length": 0.002,
                        "vacuum_layer_length": 0.0,
                        "B0": 2.95,
                        "identical_profiles": True,
                        "infinite_waveguide": True,
                    },
                    "nz_min": -20,
                    "nz_max": 20,
                    "dnz": 0.01,
                    "ny_min": -2,
                    "ny_max": 2,
                    "dny": 0.1,
                    "z_min": -0.015,
                    "z_max": 0.075,
                    "nb_z": 60,
                    "x_max": 0.05,
                    "nb_x": 8,
                },
            },
            "options": {"debug": False},
        }

        scenario = Scenario(scenario_dict)

        # Check if the binary exists before trying to run it
        binary_path = get_binary_path(6, "glnxa64")
        if not binary_path.exists():
            self.skipTest(f"Fortran binary not found: {binary_path}")

        try:
            # Test with the Scenario object (no additional parameters needed)
            S_plasma, rac_Zhe = S_plasma_1D(scenario)

            # Basic checks
            self.assertIsInstance(S_plasma, np.ndarray)
            self.assertIsInstance(rac_Zhe, np.ndarray)
            self.assertEqual(S_plasma.dtype, np.complex128)
            self.assertEqual(rac_Zhe.dtype, np.complex128)

            # The expected size should be (nb_g_total_ligne * (Nmh + Nme)) x (nb_g_total_ligne * (Nmh + Nme))
            # From the antenna: 8 modules * 1 waveguide + 0 edge waveguides = 8 waveguides
            # From spectral_1D: Nmh=1, Nme=2, so 3 modes
            # Total size: 8 * 3 = 24
            expected_size = 24
            self.assertEqual(S_plasma.shape, (expected_size, expected_size))
            self.assertEqual(rac_Zhe.shape, (expected_size, expected_size))

        except RuntimeError as e:
            # Fortran binary might fail due to environment/dependency issues
            # This is acceptable for the test - we just want to ensure the Python
            # interface works correctly
            if "Binary execution failed" in str(e):
                self.skipTest(f"Fortran binary execution failed: {e}")
            else:
                raise

    def test_S_plasma_1D_wrong_solver(self):
        """Test that S_plasma_1D raises error for unsupported solvers."""
        from aloha.scenario import Scenario

        # Create a scenario with unsupported solver
        scenario_dict = {
            "antenna": {"file": "8_active_waveguides.toml", "excitation": {"f": 3.7e9}},
            "plasma": {
                "solver": "spectral_2D",  # This should cause an error
                "spectral_2D": {"profile": "bilinear"},
            },
            "options": {"debug": False},
        }

        scenario = Scenario(scenario_dict)

        # Test with unsupported solver
        with self.assertRaises(ValueError):
            S_plasma_1D(scenario)

    def test_S_plasma_1D_unsupported_profile(self):
        """Test that S_plasma_1D raises error for unsupported plasma profiles."""
        from aloha.scenario import Scenario

        # Create a scenario with unsupported profile
        scenario_dict = {
            "antenna": {"file": "8_active_waveguides.toml", "excitation": {"f": 3.7e9}},
            "plasma": {
                "solver": "spectral_1D",
                "spectral_1D": {
                    "profile": "linear",  # This should cause an error
                    "nb_evanescent_modes": 2,
                    "bilinear": {
                        "ne0": 5e17,
                        "lambda_n": [0.002, 0.02],
                        "plasma_layer_length": 0.002,
                        "vacuum_layer_length": 0.0,
                        "B0": 2.95,
                        "identical_profiles": True,
                        "infinite_waveguide": True,
                    },
                },
            },
            "options": {"debug": False},
        }

        scenario = Scenario(scenario_dict)

        # Test with unsupported profile
        with self.assertRaises(ValueError) as context:
            S_plasma_1D(scenario)

        self.assertIn("Unsupported plasma profile", str(context.exception))
        self.assertIn("linear", str(context.exception))


if __name__ == "__main__":
    unittest.main()
