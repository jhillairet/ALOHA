import numpy as np

from aloha.constants import c, epsilon_0, mu_0, pi


class Waveguide:
    def __init__(self, a: float, b: float, L: float = 0, active: bool = True) -> None:
        """
        Rectangular waveguide.

        Parameters
        ----------
        a : float
            Length of the large side in meter.
        b : float
            Length of the small side in meter.
        L : float, optional
            Guide (longitudinal) length in meter. Default is 0.
        active : boolean, optional
            True if active (directly fed waveguide), False if passive (short-circuited).
            Default is True.

        """
        self.a = a
        self.b = b
        self.L = L
        self.active = active

    def cutoff_wavenumber(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff wavenumber for TE_mn or TM_mn mode.

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        k_c : float
            Cutoff wavenumber in 1/m.

        """
        return np.sqrt((m * pi / self.a) ** 2 + (n * pi / self.b) ** 2)

    def cutoff_frequency(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff frequency for TE_mn or TM_mn mode.

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        f_c : float
            Cutoff frequency in Hz.

        """
        return c / (2 * pi) * self.cutoff_wavenumber(m=m, n=n)

    def cutoff_wavelength(self, m: int = 1, n: int = 0) -> float:
        """
        Cutoff wavelength for TE_mn or TM_mn mode.

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        lambda_c : float
            Cutoff wavelength in m.

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

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        gamma_mn : complex
            Complex wavenumber in 1/m.

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

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        lambda_g : complex
            Guided wavelength in meter.

        """
        return 2 * pi / self.guided_wavenumber(f, m, n).imag

    def phase_velocity(self, f: float, m: int = 1, n: int = 0) -> complex:
        r"""
        Complex wave phase velocity for TE_mn or TM_mn mode (in m/s).

        .. math::
            j \cdot \omega / \gamma

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.

        Results
        -------
        v_phi : complex
            Phase velocity in m/s.

        Notes
        -----
        The `j` is used so that real phase velocity corresponds to propagation

        where:

        * :math:`\omega` is angular frequency (rad/s),
        * :math:`\gamma` is complex propagation constant (rad/m)

        """
        return 1j * 2 * pi * f / self.guided_wavenumber(f, m, n)

    def characteristic_impedance(self, f: float, m: int = 1, n: int = 0, mode: str = "te") -> complex:
        r"""
        The characteristic impedance, :math:`z_{0,mn}`.

        The characteristic impedance depends of the mode ('TE' or 'TM').

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.
        mode : str
            Electromagnetic mode. 'te' (default) or 'tm'

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

    def characteristic_admittance(self, f: float, m: int = 1, n: int = 0, mode: str = "te") -> complex:
        r"""The characteristic admittance, :math:`y_{0,mn}`.

        The characteristic admittance depends of the mode ('TE' or 'TM').

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int
            Large side mode index. Default is 1.
        n : int
            Small side mode index. Default is 0.
        mode : str
            Electromagnetic mode. 'te' (default) or 'tm'

        Returns
        -------
        y0_characteristic : np.ndarray
            Characteristic Impedance in units of ohms
        """
        return 1 / self.characteristic_impedance(f, m, n, mode)

    def electric_length(self, f: float, m: int = 1, n: int = 0, L: float | None = None) -> complex:
        r"""
        Electric length of the waveguide for the TE or TM mode mn.

        .. math::
            L_e = \gamma_{mn} L


        Expressing the physical length in terms of the phase shift
        that the wave experiences as it propagates down the line.

        The convention has been chosen that forward propagation is
        represented by the positive imaginary part of the value
        returned by the gamma function.

        Parameters
        ----------
        f : float
            Frequency of the wave in Hz.
        m : int, optional
            Large side mode index. Default is 1.
        n : int, optional
            Small side mode index. Default is 0.
        L : float, optional
            Longitudinal length of the guide. Default is the object parameter.

        Results
        -------
        L_e : complex
            Electrical length in meter.

        """
        _L = L if L else self.L
        return self.guided_wavenumber(f, m, n) * _L
