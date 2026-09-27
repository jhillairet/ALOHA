# Plasma module test
import tempfile
from pathlib import Path

import numpy as np
import pytest

from aloha.plasma import S_plasma_1D, S_plasma_1D_matlab_inputs, get_binary_name, get_binary_path
from aloha.scenario import Scenario


@pytest.fixture
def scenario_dict():
    """Fixture to provide a fresh copy of the scenario dictionary for each test."""
    return {
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


def test_get_binary_name():
    """Test binary name generation."""
    name = get_binary_name(6, "glnxa64")
    assert name == "coupl_plasma_version6_glnxa64"

    name = get_binary_name(3, "alpha")
    assert name == "coupl_plasma_version3_alpha"


def test_get_binary_path():
    """Test binary path generation."""
    path = get_binary_path(6, "glnxa64")
    expected_path = (
        Path(__file__).resolve().parent.parent.parent
        / "aloha_matlab"
        / "code_1D"
        / "couplage_1D"
        / "coupl_plasma_version6_glnxa64"
    )
    assert path == expected_path


def test_binary_exists():
    """Test that the binary exists."""
    binary_path = get_binary_path(6, "glnxa64")
    assert binary_path.exists(), f"Binary not found: {binary_path}"


def test_S_plasma_1D_wrong_version(scenario_dict):
    """Test that S_plasma_1D raises error for unsupported versions."""
    scenario = Scenario(scenario_dict)

    # Test with unsupported version
    with pytest.raises(ValueError):
        S_plasma_1D_matlab_inputs(scenario.scenario, version=3)

    with pytest.raises(ValueError):
        S_plasma_1D_matlab_inputs(scenario.scenario, version=7)


def test_S_plasma_1D_with_scenario_object(scenario_dict):
    """Test S_plasma_1D with Scenario object using TOML schema."""
    scenario = Scenario(scenario_dict)

    # Test with the Scenario object (no additional parameters needed)
    S_plasma, rac_Zhe = S_plasma_1D(scenario)

    # Basic checks
    assert isinstance(S_plasma, np.ndarray)
    assert isinstance(rac_Zhe, np.ndarray)
    assert S_plasma.dtype == np.complex128
    assert rac_Zhe.dtype == np.complex128

    # The expected size should be (nb_g_total_ligne * (Nmh + Nme)) x (nb_g_total_ligne * (Nmh + Nme))
    # From the antenna: 8 modules * 1 waveguide + 0 edge waveguides = 8 waveguides
    # From spectral_1D: Nmh=1, Nme=2, so 3 modes
    # Total size: 8 * 3 = 24
    expected_size = 24
    assert S_plasma.shape == (expected_size, expected_size)
    assert rac_Zhe.shape == (expected_size, expected_size)


def test_S_plasma_1D_wrong_solver(scenario_dict):
    """Test that S_plasma_1D raises error for unsupported solvers."""
    # Create a scenario with unsupported solver
    scenario_dict["plasma"]["solver"] = "spectral_2D"  # This should cause an error
    scenario = Scenario(scenario_dict)

    # Test with unsupported solver
    with pytest.raises(ValueError):
        S_plasma_1D(scenario)


def test_S_plasma_1D_unsupported_profile(scenario_dict):
    """Test that S_plasma_1D raises error for unsupported plasma profiles."""
    # Create a scenario with unsupported profile
    scenario_dict["plasma"]["spectral_1D"]["profile"] = "linear"  # This should cause an error
    scenario = Scenario(scenario_dict)

    # Test with unsupported profile
    with pytest.raises(ValueError) as exc_info:
        S_plasma_1D(scenario)

    assert "Unsupported plasma profile" in str(exc_info.value)
    assert "linear" in str(exc_info.value)
