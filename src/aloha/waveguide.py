import numpy as np

from aloha.constants import c, epsilon_0, mu_0, pi


class Waveguide:
    """
    Rectangular Waveguide.
    """

    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

    def cutoff_wavenumber(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff wavenumber for TE_mn or TM_mn mode.
        """
        return np.sqrt((m * pi / self.a) ** 2 + (n * pi / self.b) ** 2)

    def cutoff_frequency(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff frequency for TE_mn or TM_mn mode.
        """
        return c / (2 * pi) * self.cutoff_wavenumber(m=m, n=n)

    def cutoff_wavelength(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff wavelength for TE_mn or TM_mn mode.
        """
        return c / self.cutoff_frequency(m=m, n=n)

    def guided_wavenumber(self, f: float, m: int = 1, n: int = 0) -> complex:
        r"""
        Guided wavenumber `gamma` for TE_mn or TM_mn mode.

        `gamma` follows the convention:

        * positive real(gamma) = attenuation
        * positive imag(gamma) = forward propagation

        Defined as

        .. math::

            gamma_{mn} = \pm j \sqrt {k_0^2 - k_{c,mn}^2}

        This is:

        * IMAGINARY for propagating modes
        * REAL for non-propagating modes

        """
        k0 = 2 * pi * f / c
        kc = self.cutoff_wavenumber(m=m, n=n)
        gamma = (
            1j * np.sqrt(np.abs(k0**2 - kc**2)) * (k0 > kc)
            + np.sqrt(np.abs(kc**2 - k0**2)) * (k0 < kc)
            + 0 * (kc == k0)
        )
        return gamma

    def guided_wavelength(self, f: float, m: int = 1, n: int = 0) -> complex:
        """
        Guided wavelength for TE_mn or TM_mn mode.
        """
        return 2 * pi / self.guided_wavenumber(f, m, n).imag

    def phase_velocity(self, f: float, m: int = 1, n: int = 0) -> complex:
        r"""
        Complex wave phase velocity for TE_mn or TM_mn mode (in m/s).

        .. math::
            j \cdot \omega / \gamma

        Notes
        -----
        The `j` is used so that real phase velocity corresponds to propagation

        where:

        * :math:`\omega` is angular frequency (rad/s),
        * :math:`\gamma` is complex propagation constant (rad/m)

        Returns
        -------
        v_p : :class:`numpy.ndarray`

        """
        return 1j * 2 * pi * f / self.guided_wavenumber(f, m, n)

    def characteristic_impedance(self, f: float, m: int = 1, n: int = 1, mode: str = "te") -> complex:
        r"""
        The characteristic impedance, :math:`z_{0,mn}`.

        The characteristic impedance depends of the mode ('TE' or 'TM').

        Returns
        -------
        z0_characteristic : np.ndarray
            Characteristic Impedance in units of ohms
        """
        mode = mode.lower()
        omega = 2 * pi * f
        gamma = self.guided_wavenumber(f, m, n)
        impedance_dict = {
            "te": 1j * omega * mu_0 / gamma,
            "tm": -1j * gamma / (omega * epsilon_0),
        }
        return impedance_dict[mode]

    def characteristic_admittance(self, f: float, m: int = 1, n: int = 1, mode: str = "te") -> complex:
        r"""The characteristic admittance, :math:`y_{0,mn}`.

        The characteristic admittance depends of the mode ('TE' or 'TM').

        Returns
        -------
        y0_characteristic : np.ndarray
            Characteristic Impedance in units of ohms
        """
        return 1 / self.characteristic_impedance(f, m, n, mode)
