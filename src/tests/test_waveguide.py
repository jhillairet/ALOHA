"""
Testing for the Waveguide class
"""

import unittest

import numpy as np
from scipy.signal import gammatone
from skrf import Frequency
from skrf.media import RectangularWaveguide

from aloha.constants import Z0, c, pi
from aloha.waveguide import Waveguide


class WaveguideTests(unittest.TestCase):
    def test_waveguide_constructor(self):
        wg = Waveguide(a=1e-6, b=0.5e-6)
        assert wg.a == 1e-6
        assert wg.b == 0.5e-6

    def test_waveguide_properties(self):
        "Test various waveguide properties for a few modes"
        a, b = 76e-3, 14e-3
        f0 = 3.7e9
        wg = Waveguide(a, b)
        # scikit-rf reference
        freq = Frequency(f0, f0, npoints=1, unit="Hz")
        for m in range(1, 3):
            for n in range(0, 3):
                for mode in ("te", "tm"):
                    # warning is not logged below
                    # to hide the warnings which appears 1/0 (evanescent modes)
                    with np.errstate(divide="ignore"):
                        wg_skrf = RectangularWaveguide(freq, a=a, b=b, m=m, n=n, mode_type=mode)
                        # cut-off properties
                        assert np.isclose(wg.cutoff_wavenumber(m, n), wg_skrf.kc)
                        assert np.isclose(wg.cutoff_frequency(m, n), wg_skrf.f_cutoff)
                        assert np.isclose(wg.cutoff_wavelength(m, n), wg_skrf.lambda_cutoff)
                        # guided properties
                        assert np.isclose(wg.guided_wavenumber(f0, m, n), wg_skrf.gamma)
                        assert np.isclose(wg.guided_wavelength(f0, m, n), wg_skrf.lambda_guide)
                        assert np.isclose(wg.phase_velocity(f0, m, n), wg_skrf.v_p)
                        assert np.isclose(wg.characteristic_impedance(f0, m, n, mode), wg_skrf.z0_characteristic)
                        assert np.isclose(wg.characteristic_admittance(f0, m, n, mode), 1 / wg_skrf.z0_characteristic)

    def test_characteristic_admittance(self):
        """
        Test characteristic admittances for TE and TM modes.

        Formulas from ALOHA doc:
        Y_c_mn = gamma_mn / (j k_0 Z_0) for TE modes
        Y_c_mn = j k_0 / (gamma_mn Z_0) fpr TM modes
        """
        a, b = 76e-3, 14e-3
        f0 = 3.7e9
        wg = Waveguide(a, b)
        k0 = 2 * pi * f0 / c

        for m in range(1, 3):
            for n in range(0, 3):
                gamma_mn = wg.guided_wavenumber(f0, m, n)
                for mode in ("te", "tm"):
                    if mode == "te":
                        Y_c_mn = gamma_mn / (1j * k0 * Z0)
                    else:
                        Y_c_mn = 1j * k0 / (gamma_mn * Z0)

                    assert np.isclose(wg.characteristic_admittance(f0, m, n, mode), Y_c_mn)
