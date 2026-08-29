"""
ALOHA Plasma coupling calculation module.
"""

import subprocess
from pathlib import Path

import numpy as np
from scipy.constants import c, epsilon_0, m_e
from scipy.constants import e as q_e


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

    Parameters
    ----------
    filename : str or Path
        Output filename
    parameters : list
        List of parameter values to write
    """
    with open(filename, "w") as f:
        for param in parameters:
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

        # Skip calculation if lines are identical and this is not the first row
        if bool_lignes_identiques and ind > 0:
            if bool_debug:
                print(f"{text_ind}[bool_lignes_identiques=true]")  # noqa: T201
                print(f"{text_ind}--> No need to calculate this row! taking previous result")  # noqa: T201
            continue

        if bool_debug:
            print(f"{text_ind}Writing parameters file for the binary")  # noqa: T201

        # Prepare parameters for version 6
        # Order: Nmh, Nme, freq, ne0(ind), dne0(ind), d_couche(ind), dne1(ind),
        #        nb_g_total_ligne, a, b, z, T_grill, D_guide_max, erreur_rel, pertes, max_nz, d_vide(ind)

        # For version 6, we need to write arrays element by element
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
        ]

        # Then write b array
        var_list.extend([float(bi) for bi in b])

        # Then write z array
        var_list.extend([float(zi) for zi in z])

        # Then write remaining scalar parameters
        var_list.extend(
            [float(T_grill), float(D_guide_max), float(erreur_rel), float(pertes), float(max_nz), float(d_vide[ind])]
        )

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

        # Read results
        if bool_debug:
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
