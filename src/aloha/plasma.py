"""
ALOHA Plasma coupling calculation module.
"""

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.constants import c, epsilon_0, m_e
from scipy.constants import e as q_e

if TYPE_CHECKING:
    from aloha.scenario import Scenario


def get_binary_name(version: int, architecture: str = "glnxa64") -> str:
    """
    Determine the Fortran binary name based on version and architecture.

    Parameters
    ----------
    version : int
        Version number of the binary (e.g., 6)
    architecture : str, optional
        Target architecture (default: 'glnxa64')

    Returns
    -------
    str
        Binary filename string
    """
    return f"coupl_plasma_version{version}_{architecture}"


def get_binary_path(version: int, architecture: str = "glnxa64") -> Path:
    """
    Get the full path to the Fortran binary.

    Parameters
    ----------
    version : int
        Version number of the binary
    architecture : str, optional
        Target architecture

    Returns
    -------
    Path
        Path to the binary executable
    """
    # Assuming the binaries are in the aloha_matlab/code_1D/couplage_1D directory
    # relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    binary_dir = project_root / "aloha_matlab" / "code_1D" / "couplage_1D"
    binary_name = get_binary_name(version, architecture)
    return binary_dir / binary_name


def save_fortran_input_file(filename: str | Path, parameters: list) -> None:
    """
    Save parameters to a Fortran input file.

    Each parameter is written on a separate line with scientific notation.
    If a parameter is an array-like object (numpy array or list), all its elements
    are written on the same line with space separation, matching the MATLAB behavior.

    Parameters
    ----------
    filename : str or Path
        Output filename
    parameters : list
        List of parameter values to write. Can contain scalars, numpy arrays, or lists.
    """
    with open(filename, "w") as f:
        for param in parameters:
            # Check if the parameter is array-like (numpy array or list)
            if isinstance(param, (np.ndarray, list)):
                # Write all elements of the array on the same line with space separation
                param_array = np.asarray(param)
                line = "  " + "  ".join([f"{x:1.7e}" for x in param_array.flat])
                f.write(line + "\n")
            else:
                # Write scalar parameter on its own line
                f.write(f"  {param:1.7e}\n")


def read_fortran_output_file(filename: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read Fortran output file containing S_plasma, rac_Zhe, and K_cpl matrices.

    The file format is:
    - First 2 lines: header
    - Remaining lines: S(i,j), RZheg(i,j), K_cpl(i,j) for each i,j

    Parameters
    ----------
    filename : str or Path
        Input filename

    Returns
    -------
    tuple of np.ndarray
        Tuple of (S_plasma, rac_Zhe, K_cpl) as complex numpy arrays
    """
    # Read all lines
    with open(filename) as f:
        lines = f.readlines()

    # Skip header lines (first 2 lines)
    data_lines = lines[2:]

    # Parse complex numbers from each line
    # Format: (S_real,S_imag) (RZheg_real,RZheg_imag) (K_cpl_real,K_cpl_imag)
    S_real = []
    S_imag = []
    RZheg_real = []
    RZheg_imag = []
    K_cpl_real = []
    K_cpl_imag = []

    for line in data_lines:
        line = line.strip()
        if not line:
            continue

        # Parse the line with format: (r1,i1) (r2,i2) (r3,i3)
        # Remove parentheses and split
        line_clean = line.replace("(", "").replace(")", "").replace(",", " ")
        parts = line_clean.split()

        if len(parts) >= 6:
            S_real.append(float(parts[0]))
            S_imag.append(float(parts[1]))
            RZheg_real.append(float(parts[2]))
            RZheg_imag.append(float(parts[3]))
            K_cpl_real.append(float(parts[4]))
            K_cpl_imag.append(float(parts[5]))

    # Convert to numpy arrays
    S_plasma = np.array(S_real) + 1j * np.array(S_imag)
    rac_Zhe = np.array(RZheg_real) + 1j * np.array(RZheg_imag)
    K_cpl = np.array(K_cpl_real) + 1j * np.array(K_cpl_imag)

    # Reshape to square matrices
    # The matrices are (Nme+Nmh)*nb_g_total_ligne x (Nme+Nmh)*nb_g_total_ligne
    n = int(np.sqrt(len(S_plasma)))
    S_plasma = S_plasma.reshape((n, n))
    rac_Zhe = rac_Zhe.reshape((n, n))
    K_cpl = K_cpl.reshape((n, n))

    return S_plasma, rac_Zhe, K_cpl


def S_plasma_1D_matlab_inputs(
    scenario: dict,
    version: int = 6,
    architecture: str = "glnxa64",
    bool_lignes_identiques: bool = False,
    bool_debug: bool = False,
    max_nz: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the plasma S-parameter matrix using ALOHA Fortran binary.

    This function only implements version 6 of the plasma coupling calculation.
    Version 6: linear profile, two plasma layers

    Parameters
    ----------
    scenario : dict
        Dictionary containing scenario parameters. Required keys:
        - antenna.freq: Frequency in Hz
        - ne0: Electron density at reference position (array)
        - dne0: Electron density gradient (array)
        - d_couche: Thickness of first plasma layer (array)
        - dne1: Electron density gradient for second layer (array)
        - nb_g_pol: Number of poloidal rows
        - nb_g_total_ligne: Number of waveguides per poloidal row
        - a: Waveguide width parameter
        - b: Waveguide height parameters (array)
        - z: Waveguide position parameters (array)
        - T_grill: Grill periodicity parameter
        - D_guide_max: Maximum guide decoupling distance
        - erreur_rel: Relative error tolerance
        - pertes: Loss parameter
        - d_vide: Vacuum layer thickness (array)
        - Nmh: Number of magnetic modes
        - Nme: Number of electric modes
    version : int, optional
        Version of the calculation (only 6 is implemented)
    architecture : str, optional
        Target architecture for binary selection
    bool_lignes_identiques : bool, optional
        If True, all poloidal rows are identical
    bool_debug : bool, optional
        If True, print debug information
    max_nz : int, optional
        Maximum integration limit for nz

    Returns
    -------
    tuple of np.ndarray
        Tuple of (S_plasma, rac_Zhe) numpy arrays

    Raises
    ------
    ValueError
        If version is not 6
    RuntimeError
        If Fortran binary execution fails
    """
    if version != 6:
        raise ValueError(f"Only version 6 is implemented. Got version {version}")

    # Extract parameters from scenario
    freq = scenario["antenna"]["freq"]
    ne0 = np.array(scenario["ne0"])
    dne0 = np.array(scenario["dne0"])
    d_couche = np.array(scenario["d_couche"])
    dne1 = np.array(scenario["dne1"])
    nb_g_pol = scenario["nb_g_pol"]
    nb_g_total_ligne = scenario["nb_g_total_ligne"]
    a = scenario["a"]
    b = np.array(scenario["b"])
    z = np.array(scenario["z"])
    T_grill = scenario["T_grill"]
    D_guide_max = scenario["D_guide_max"]
    erreur_rel = scenario["erreur_rel"]
    pertes = scenario["pertes"]
    d_vide = np.array(scenario["d_vide"])
    Nmh = scenario["Nmh"]
    Nme = scenario["Nme"]

    # Physical constants (currently unused but kept for potential future use)
    # k0 = 2 * np.pi * freq / c
    # nc = (2 * np.pi * freq) ** 2 * m_e * epsilon_0 / q_e**2

    # Initialize output matrices
    matrix_size = nb_g_total_ligne * nb_g_pol * (Nmh + Nme)
    S_plasma = np.zeros((matrix_size, matrix_size), dtype=complex)
    rac_Zhe = np.zeros((matrix_size, matrix_size), dtype=complex)

    # Get binary path
    binary_path = get_binary_path(version, architecture)

    if not binary_path.exists():
        raise FileNotFoundError(f"Fortran binary not found: {binary_path}")

    # Working directory for Fortran execution
    # Use the directory containing the binary
    working_dir = binary_path.parent

    # Input and output filenames
    fortran_input_file = working_dir / "par_grill.dat"
    fortran_output_file = working_dir / "S_plasma2.dat"

    # Process each poloidal row
    for ind in range(nb_g_pol):
        text_ind = f"[row #{ind + 1}/{nb_g_pol}]: "

        # Skip binary execution if lines are identical and this is not the first row
        # But still need to read results and update matrices
        if bool_lignes_identiques and ind > 0:
            if bool_debug:
                print(f"{text_ind}[bool_lignes_identiques=true]")  # noqa: T201
                print(f"{text_ind}--> No need to calculate this row! taking previous result")  # noqa: T201
            # Skip binary execution but continue to read results
            skip_binary = True
        else:
            skip_binary = False

        if not skip_binary:
            if bool_debug:
                print(f"{text_ind}Writing parameters file for the binary")  # noqa: T201

            # Prepare parameters for version 6
            # Order: Nmh, Nme, freq, ne0(ind), dne0(ind), d_couche(ind), dne1(ind),
            #        nb_g_total_ligne, a, b, z, T_grill, D_guide_max, erreur_rel, pertes, max_nz, d_vide(ind)

            # For version 6, we need to pass arrays b and z as arrays, not flattened
            # First, write scalar parameters
            var_list = [
                float(Nmh),
                float(Nme),
                float(freq),
                float(ne0[ind]),
                float(dne0[ind]),
                float(d_couche[ind]),
                float(dne1[ind]),
                float(nb_g_total_ligne),
                float(a),
                b,  # Keep as numpy array for proper formatting
                z,  # Keep as numpy array for proper formatting
                float(T_grill),
                float(D_guide_max),
                float(erreur_rel),
                float(pertes),
                float(max_nz),
                float(d_vide[ind]),
            ]

            # Write parameters to file
            save_fortran_input_file(fortran_input_file, var_list)

            # Execute Fortran binary
            if bool_debug:
                print(f"{text_ind}Run binary {binary_path}")  # noqa: T201

            try:
                result = subprocess.run([str(binary_path)], cwd=working_dir, capture_output=True, text=True, check=True)

                if bool_debug:
                    print(result.stdout)  # noqa: T201
                    if result.stderr:
                        print(result.stderr)  # noqa: T201

            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Binary execution failed: {e}\nStdout: {e.stdout}\nStderr: {e.stderr}") from e

        # Read results (even if we skipped binary execution, we still need to read the previous results)
        if bool_debug:
            if skip_binary:
                print(f"{text_ind}Reading previous result file")  # noqa: T201
            else:
                print(f"{text_ind}Reading ascii result file")  # noqa: T201

        if not fortran_output_file.exists():
            raise FileNotFoundError(f"Output file not created: {fortran_output_file}")

        S_plasma_line, rac_Zhe_line, K_cpl = read_fortran_output_file(fortran_output_file)

        # Update the output matrices
        # For version 6, the matrices are (Nme+Nmh)*nb_g_total_ligne x (Nme+Nmh)*nb_g_total_ligne

        # Create a mask for this poloidal row
        tableau = np.zeros(nb_g_pol)
        tableau[ind] = 1

        # Update matrices using kronecker product
        S_plasma += np.kron(np.diag(tableau), S_plasma_line)
        rac_Zhe += np.kron(np.diag(tableau), rac_Zhe_line)

        if bool_debug:
            print(f"{text_ind}Results reading: OK")  # noqa: T201
            print(f"{text_ind}Sum(S_plasma)={np.sum(S_plasma)}")  # noqa: T201
            print(f"{text_ind}Sum(rac_Zhe)={np.sum(rac_Zhe)}")  # noqa: T201

    return S_plasma, rac_Zhe


def S_plasma_1D(scenario: "Scenario") -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the plasma S-parameter matrix using ALOHA Fortran binary with a Scenario object.

    This function takes a Scenario object following the TOML schema and converts it to the
    MATLAB-style parameter format expected by S_plasma_1D_matlab_inputs.

    Parameters
    ----------
    scenario : Scenario
        A Scenario object containing scenario parameters in TOML schema format.
        The function extracts the following parameters from the scenario:
        - version: from plasma.spectral_1D.profile (only 'bilinear' supported)
        - architecture: default 'glnxa64'
        - bool_lignes_identiques: from plasma.spectral_1D.bilinear.identical_profiles
        - bool_debug: from options.debug
        - max_nz: from plasma.spectral_1D.max_nz or default 100

    Returns
    -------
    tuple of np.ndarray
        Tuple of (S_plasma, rac_Zhe) numpy arrays

    Raises
    ------
    ValueError
        If scenario conversion fails
    RuntimeError
        If Fortran binary execution fails
    """
    # Import the conversion function from scenario module
    from aloha.scenario import _convert_scenario_to_matlab_inputs

    # Extract parameters from scenario
    scenario_dict = scenario.scenario

    # Get version from plasma.spectral_1D.profile
    # For now, we only support version 6 (linear profile)
    # The profile type is mapped to version: 'bilinear' -> version 6
    version = 6  # Default version for bilinear profile

    # Get architecture (default to glnxa64)
    architecture = "glnxa64"

    # Get bool_lignes_identiques from plasma.spectral_1D.bilinear.identical_profiles
    plasma = scenario_dict.get("plasma", {})
    spectral_1D = plasma.get("spectral_1D", {})
    bilinear = spectral_1D.get("bilinear", {})
    bool_lignes_identiques = bilinear.get("identical_profiles", False)

    # Get bool_debug from options.debug
    options = scenario_dict.get("options", {})
    bool_debug = options.get("debug", False)

    # Get max_nz from plasma.spectral_1D or default to 100
    max_nz = spectral_1D.get("max_nz", 100)

    # Convert the Scenario object to MATLAB-style parameter dictionary
    matlab_scenario = _convert_scenario_to_matlab_inputs(scenario)

    # Call the existing S_plasma_1D_matlab_inputs function with the converted parameters
    return S_plasma_1D_matlab_inputs(
        matlab_scenario,
        version=version,
        architecture=architecture,
        bool_lignes_identiques=bool_lignes_identiques,
        bool_debug=bool_debug,
        max_nz=max_nz,
    )


def compute_spectrum1D(scenario: "Scenario") -> None:
    """
    Compute the power spectrum of the antenna (1D version).

    This function implements the logic from MATLAB's aloha_compute_spectrum1D.m to
    calculate the power spectrum dP and its nz-integrated version dP_nz.

    Parameters
    ----------
    scenario : Scenario
        A Scenario object containing the scenario parameters and results.
        Must have the following in scenario.results:
        - S_plasma: Plasma scattering matrix
        - rac_Zhe: Characteristic impedance matrix
        - a_plasma: Incident wave vector on plasma
        - b_plasma: Reflected wave vector from plasma

    Results are stored in scenario.results:
        - ny: Array of ny values
        - nz: Array of nz values
        - nz0: Main nz peak value
        - dny: Step in ny
        - dnz: Step in nz
        - dP: Power spectrum density (2D array: ny x nz)
        - dP_nz: Power spectrum integrated over ny (1D array: nz)

    Raises
    ------
    ValueError
        If required parameters are missing from scenario.results or scenario.scenario
    """
    _compute_spectrum1D_impl(scenario)


def _compute_spectrum_modes(
    scenario,
    ny,
    nz,
    nyt,
    nzt,
    k0,
    nc_val,
    X0,
    D0,
    X1,
    D1,
    nb_g_pol,
    nb_g_total_ligne,
    a,
    b,
    y,
    z,
    Nme,
    Nmh,
    ne0_arr,
    lambda_n_arr,
    d_couche_arr,
    d_vide_arr,
    version,
    poids_E,
    poids_H,
    Y0,
):
    """Compute the spectrum modes and return dP matrix."""
    import numpy as np
    from scipy.special import airy

    modes_E_sp_list = []
    modes_H_sp_list = []

    if version == 6:
        for ind in range(nb_g_pol):
            neta_1 = (
                ((nzt**2 - 1) ** (1 / 3)) * ((D1[ind] / X1[ind]) ** (2 / 3)) * (X1[ind] - 1) * np.exp(-1j * np.pi / 3)
            )

            airy_0_neta_1, airy_1_neta_1 = airy(neta_1)[0:2]
            airy_0_neta_1[np.abs(airy_0_neta_1) < 1e-30] = 1e-30

            # Avoid division by zero for (nzt**2 - 1)**(2/3)
            nzt_squared_minus_1 = nzt**2 - 1
            nzt_squared_minus_1[np.abs(nzt_squared_minus_1) < 1e-30] = 1e-30
            nzt_squared_minus_1_power_2_3 = nzt_squared_minus_1 ** (2 / 3)
            nzt_squared_minus_1_power_2_3[np.abs(nzt_squared_minus_1_power_2_3) < 1e-30] = 1e-30

            Y_L = (
                -(airy_1_neta_1 / airy_0_neta_1)
                / nzt_squared_minus_1_power_2_3
                * ((X1[ind] / D1[ind]) ** (1 / 3))
                * np.exp(1j * np.pi / 6)
            )

            neta_0 = (
                (nzt_squared_minus_1 ** (1 / 3))
                * ((D0[ind] / X0[ind]) ** (2 / 3))
                * (X0[ind] - 1)
                * np.exp(-1j * np.pi / 3)
            )

            airy_0_neta_0, airy_1_neta_0, airy_2_neta_0, airy_3_neta_0 = airy(neta_0)[0:4]
            airy_0_neta_0[np.abs(airy_0_neta_0) < 1e-30] = 1e-30
            airy_2_neta_0[np.abs(airy_2_neta_0) < 1e-30] = 1e-30

            Ya = (
                -(airy_1_neta_0 / airy_0_neta_0)
                / nzt_squared_minus_1_power_2_3
                * ((X0[ind] / D0[ind]) ** (1 / 3))
                * np.exp(1j * np.pi / 6)
            )
            Yb = (
                -(airy_3_neta_0 / airy_2_neta_0)
                / nzt_squared_minus_1_power_2_3
                * ((X0[ind] / D0[ind]) ** (1 / 3))
                * np.exp(1j * np.pi / 6)
            )

            neta_d = neta_0 + (nzt_squared_minus_1 ** (1 / 3)) * ((X0[ind] / D0[ind]) ** (1 / 3)) * k0 * d_couche_arr[
                ind
            ] * np.exp(-1j * np.pi / 3)

            airy_0_neta_d, airy_1_neta_d, airy_2_neta_d, airy_3_neta_d = airy(neta_d)[0:4]

            # neta_d and neta_0 have the same length, so no tiling is needed
            # Use the airy values directly
            airy_0_neta_0_tiled = airy_0_neta_0
            airy_2_neta_0_tiled = airy_2_neta_0
            airy_1_neta_0_tiled = airy_1_neta_0
            airy_3_neta_0_tiled = airy_3_neta_0

            airy_0_neta_0_tiled[np.abs(airy_0_neta_0_tiled) < 1e-30] = 1e-30
            airy_2_neta_0_tiled[np.abs(airy_2_neta_0_tiled) < 1e-30] = 1e-30

            num = Y_L * (
                Yb * airy_0_neta_d / airy_0_neta_0_tiled - Ya * airy_2_neta_d / airy_2_neta_0_tiled
            ) + Ya * Yb * (airy_3_neta_d / airy_3_neta_0_tiled - airy_1_neta_d / airy_1_neta_0_tiled)

            denom = (
                Y_L * (airy_0_neta_d / airy_0_neta_0_tiled - airy_2_neta_d / airy_2_neta_0_tiled)
                - Ya * airy_1_neta_d / airy_1_neta_0_tiled
                + Yb * airy_3_neta_d / airy_3_neta_0_tiled
            )

            # Avoid division by zero
            denom[np.abs(denom) < 1e-30] = 1e-30
            # Handle NaN/Inf in num and denom before division
            # num[np.isnan(num) | np.isinf(num)] = 0  # REMOVED: let NaN/Inf propagate
            # denom[np.isnan(denom) | np.isinf(denom)] = 1e-30  # REMOVED: let NaN/Inf propagate
            denom[np.abs(denom) < 1e-30] = 1e-30
            Ys = num / denom
            # Handle any remaining NaN/Inf
            # Ys[np.isnan(Ys) | np.isinf(Ys)] = 0  # REMOVED: let NaN/Inf propagate

            alph = 1.0
            gamm = k0 * (nzt**2 - alph) ** (1 / 2)
            # Avoid division by zero for gamm
            gamm[np.abs(gamm) < 1e-30] = 1e-30
            u = gamm * d_vide_arr[ind]
            tanh_d_vide = np.tanh(u)

            y2 = 1j * k0 / gamm
            # Avoid division by zero in the final Ys calculation
            denominator_Ys = Ys * tanh_d_vide + y2
            denominator_Ys[np.abs(denominator_Ys) < 1e-30] = 1e-30
            Ys = y2 * (Ys + y2 * tanh_d_vide) / denominator_Ys

            for k in range(nb_g_total_ligne):
                m = 1
                n = 0

                coeff_h = -np.sqrt(2 / (a * b[k]))

                term1 = 1 - ((-1) ** m) * np.exp(-1j * k0 * nyt * a)
                denominator_fy = nyt**2 - (m * np.pi / (k0 * a)) ** 2
                denominator_fy[np.abs(denominator_fy) < 1e-30] = 1e-30
                fy = -(m * np.pi / (a * (k0**2))) * (term1 / denominator_fy) * np.exp(-1j * k0 * nyt * y[ind])

                term2 = 1 - ((-1) ** n) * np.exp(-1j * k0 * nzt * b[k])
                denominator_fz = nzt**2 - (n * np.pi / (k0 * b[k])) ** 2
                denominator_fz[np.abs(denominator_fz) < 1e-30] = 1e-30
                fz = -1j * (nzt / k0) * (term2 / denominator_fz) * np.exp(-1j * k0 * nzt * z[k])

                modes_E_sp_list.append(coeff_h * fy * fz)
                modes_H_sp_list.append(-Y0 * Ys * coeff_h * fy * fz)

                for n_val in range(1, Nme + 1):
                    coeff_e = -(2 * n_val / np.sqrt(b[k] / a + (n_val**2) * a / b[k])) / b[k]

                    term2_tm = 1 - ((-1) ** n_val) * np.exp(-1j * k0 * nzt * b[k])
                    denominator_fz_tm = nzt**2 - (n_val * np.pi / (k0 * b[k])) ** 2
                    denominator_fz_tm[np.abs(denominator_fz_tm) < 1e-30] = 1e-30
                    fz_tm = -1j * (nzt / k0) * (term2_tm / denominator_fz_tm) * np.exp(-1j * k0 * nzt * z[k])

                    modes_E_sp_list.append(coeff_e * fy * fz_tm)
                    modes_H_sp_list.append(-Y0 * Ys * coeff_e * fy * fz_tm)
    else:
        raise ValueError(f"Plasma version {version} not supported for spectrum calculation")

    modes_E_sp = np.array(modes_E_sp_list)
    modes_H_sp = np.array(modes_H_sp_list)

    # MATLAB: E_spc = modes_E_sp.*(poids_E*ones(1,length(nyt)))
    # poids_E is a column vector with length equal to number of modes
    # ones(1,length(nyt)) is a row vector
    # The result is a matrix where each row is a copy of poids_E

    # In Python: reshape poids_E to (N, 1) to match MATLAB's column vector behavior
    # Then broadcast: modes_E_sp (N, M) * poids_E.reshape(N, 1) -> (N, M)
    # This multiplies each row of modes_E_sp by the corresponding element of poids_E
    # NOTE: MATLAB uses poids_E for BOTH E_spc and H_spc (line 228-229)
    poids_E_reshaped = poids_E.reshape(-1, 1)

    E_spc = modes_E_sp * poids_E_reshaped
    H_spc = modes_H_sp * poids_E_reshaped

    sum_E_spc = np.sum(E_spc, axis=0)
    sum_H_spc = np.sum(H_spc, axis=0)
    dP_ligne = -0.5 * ((k0 / (2 * np.pi)) ** 2) * (sum_E_spc * np.conj(sum_H_spc))

    dP = np.reshape(dP_ligne, (len(ny), len(nz)))
    # Transpose to match MATLAB reference shape (nz, ny)
    # MATLAB's dP has shape (length(nz), length(ny)), not (length(ny), length(nz))
    dP = dP.T

    return dP


def _compute_spectrum1D_impl(scenario: "Scenario") -> None:
    """Implementation of compute_spectrum1D."""
    import numpy as np
    from scipy.constants import e as qe
    from scipy.constants import mu_0

    # Extract scenario parameters
    scenario_dict = scenario.scenario
    results = scenario.results

    # Calculate Y0 (characteristic admittance of free space)
    Z0 = np.sqrt(mu_0 / epsilon_0)
    Y0 = 1 / Z0

    # Check if required results exist
    required_results = ["S_plasma", "rac_Zhe", "a_plasma", "b_plasma"]
    for field in required_results:
        if field not in results:
            raise ValueError(f"Required field '{field}' not found in scenario.results")

    # Extract spectrum parameters
    plasma = scenario_dict.get("plasma", {})
    spectral_1D = plasma.get("spectral_1D", {})
    bilinear = spectral_1D.get("bilinear", {})

    nz_min = spectral_1D.get("nz_min", -20)
    nz_max = spectral_1D.get("nz_max", 20)
    dnz = spectral_1D.get("dnz", 0.01)
    ny_min = spectral_1D.get("ny_min", -2)
    ny_max = spectral_1D.get("ny_max", 2)
    dny = spectral_1D.get("dny", 0.1)

    # Get frequency
    antenna = scenario_dict.get("antenna", {})
    excitation = antenna.get("excitation", {})
    freq = excitation.get("f", 3.7e9)

    # Convert scenario to get geometry
    from aloha.scenario import _convert_scenario_to_matlab_inputs

    matlab_scenario = _convert_scenario_to_matlab_inputs(scenario)

    # Extract geometry
    a = matlab_scenario.get("a", 1.0)
    b = matlab_scenario.get("b", [])
    h = matlab_scenario.get("h", [])
    z = matlab_scenario.get("z", [])
    y = matlab_scenario.get("y", [])
    nb_g_total_ligne = matlab_scenario.get("nb_g_total_ligne", 0)
    nb_g_pol = matlab_scenario.get("nb_g_pol", 1)

    # Ensure y has the correct length (nb_g_pol elements)
    # Convert to lists first to handle the extension
    y_list = list(y) if not isinstance(y, list) else y
    z_list = list(z) if not isinstance(z, list) else z
    b_list = list(b) if not isinstance(b, list) else b

    # Ensure arrays have the correct lengths
    if len(y_list) < nb_g_pol:
        y_list = y_list + [0.0] * (nb_g_pol - len(y_list))
    elif len(y_list) == 0:
        y_list = [0.0] * nb_g_pol

    if len(z_list) < nb_g_total_ligne:
        z_list = z_list + [0.0] * (nb_g_total_ligne - len(z_list))
    elif len(z_list) == 0:
        z_list = [0.0] * nb_g_total_ligne

    if len(b_list) < nb_g_total_ligne:
        b_list = (
            b_list + [b_list[-1]] * (nb_g_total_ligne - len(b_list)) if len(b_list) > 0 else [0.008] * nb_g_total_ligne
        )
    elif len(b_list) == 0:
        b_list = [0.008] * nb_g_total_ligne

    # Convert to numpy arrays for easier indexing
    y = np.array(y_list)
    z = np.array(z_list)
    b = np.array(b_list)

    # Get plasma parameters
    ne0 = bilinear.get("ne0", 5e17)
    lambda_n = bilinear.get("lambda_n", [0.002, 0.02])
    d_couche = bilinear.get("plasma_layer_length", 0.002)
    d_vide = bilinear.get("vacuum_layer_length", 0.0)
    bool_lignes_identiques = bilinear.get("identical_profiles", False)
    version = bilinear.get("version", 6)

    # Get number of modes
    Nme = matlab_scenario.get("Nme", 1)
    Nmh = matlab_scenario.get("Nmh", 1)

    # Extract results
    S_plasma = results["S_plasma"]
    rac_Zhe = results["rac_Zhe"]
    a_plasma = results["a_plasma"]
    b_plasma = results["b_plasma"]

    # Ensure a_plasma and b_plasma are column vectors
    if a_plasma.ndim == 1:
        a_plasma = a_plasma.reshape(-1, 1)
    if b_plasma.ndim == 1:
        b_plasma = b_plasma.reshape(-1, 1)

    # Physical constants
    k0 = 2 * np.pi * freq / c
    nc_val = ((2 * np.pi * freq) ** 2) * m_e * epsilon_0 / (qe**2)

    # Compute modes and spectrum
    # Translate MATLAB aloha_compute_spectrum1D.m logic

    # Create ny and nz arrays
    # MATLAB: ny = ny_min:dny:ny_max, nz = nz_min:dnz:nz_max
    # Use np.arange with the correct endpoint handling
    ny = np.arange(ny_min, ny_max + dny / 2, dny)
    nz = np.arange(nz_min, nz_max + dnz / 2, dnz)

    # Avoid pure 0 and +/-1 discontinuities by adding a small real part
    # Use a tolerance for floating-point comparison
    nz = np.where(np.abs(nz) < 1e-10, 1e-4, nz)
    nz = np.where(np.abs(nz - 1) < 1e-10, 1 - 1e-4, nz)
    nz = np.where(np.abs(nz + 1) < 1e-10, -1 + 1e-4, nz)

    # Create nyt and nzt arrays (kron in MATLAB)
    nyt = np.kron(np.ones(len(nz)), ny)
    nzt = np.kron(nz, np.ones(len(ny)))

    # Handle identical profiles
    if bool_lignes_identiques:
        ne0_arr = np.ones(nb_g_pol) * ne0
        lambda_n_arr = np.ones(nb_g_pol) * lambda_n
        d_couche_arr = np.ones(nb_g_pol) * d_couche
        d_vide_arr = np.ones(nb_g_pol) * d_vide
    else:
        ne0_arr = np.array(ne0) if isinstance(ne0, list) else np.ones(nb_g_pol) * ne0
        lambda_n_arr = np.array(lambda_n) if isinstance(lambda_n, list) else np.ones(nb_g_pol) * lambda_n
        d_couche_arr = np.array(d_couche) if isinstance(d_couche, list) else np.ones(nb_g_pol) * d_couche
        d_vide_arr = np.array(d_vide) if isinstance(d_vide, list) else np.ones(nb_g_pol) * d_vide

    # Calculate X0, D0, X1, D1 for each polarization line
    X0 = ne0_arr / nc_val

    # Calculate dne0 and dne1 following MATLAB logic
    if len(lambda_n_arr) > 0:
        dne0_arr = ne0_arr / lambda_n_arr[0]
        D0 = k0 * ne0_arr / dne0_arr
    else:
        dne0_arr = ne0_arr / 0.002
        D0 = k0 * ne0_arr / dne0_arr

    if len(lambda_n_arr) > 0:
        n1 = ne0_arr + lambda_n_arr[0] * d_couche_arr
    else:
        n1 = ne0_arr + 0.002 * d_couche_arr

    X1 = n1 / nc_val

    if len(lambda_n_arr) > 1:
        # MATLAB: dne1 = (1 + d_couche/lambda_n0) * ne0 / lambda_n1
        dne1_arr = (1 + d_couche_arr / lambda_n_arr[0]) * ne0_arr / lambda_n_arr[1]
        D1 = k0 * n1 / dne1_arr
    else:
        D1 = k0 * n1 / 0.02

    # Calculate weights from power waves (MATLAB lines 53-57)
    poids_E = rac_Zhe @ (a_plasma + b_plasma)
    poids_H = np.linalg.inv(rac_Zhe) @ (a_plasma - b_plasma)

    poids_E = poids_E * np.sqrt(2)
    poids_H = poids_H * np.sqrt(2)

    # Call the mode calculation function
    dP = _compute_spectrum_modes(
        scenario,
        ny,
        nz,
        nyt,
        nzt,
        k0,
        nc_val,
        X0,
        D0,
        X1,
        D1,
        nb_g_pol,
        nb_g_total_ligne,
        a,
        b,
        y,
        z,
        Nme,
        Nmh,
        ne0_arr,
        lambda_n_arr,
        d_couche_arr,
        d_vide_arr,
        version,
        poids_E,
        poids_H,
        Y0,
    )

    # Handle NaN/Inf in dP first
    dP = np.where(np.isnan(dP) | np.isinf(dP), 1e-10, dP)

    # Calculate dP_nz (sum over ny dimension)
    dP_nz = (dny * np.sum(dP, axis=1)).reshape(-1, 1)

    # Find main nz peak
    max_nz_idx = np.argmax(np.real(dP_nz))
    nz0 = float(nz[max_nz_idx].real)

    # Store results
    scenario.results["ny"] = ny
    scenario.results["nz"] = np.real(nz)
    scenario.results["nz0"] = nz0
    scenario.results["dny"] = dny
    scenario.results["dnz"] = dnz
    scenario.results["dP"] = dP
    scenario.results["dP_nz"] = dP_nz
