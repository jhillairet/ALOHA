import os
import tomllib

import matplotlib.pyplot as plt
import numpy as np

from .constants import pi
from .waveguide import Waveguide


class Antenna:
    def __init__(self, source: str | os.PathLike | dict | None = None) -> None:
        """
        ALOHA antenna description.

        Parameters
        ----------
        source : str, Path, dict, or None
            Either a path to a TOML ALOHA antenna file,
            a dictionary containing antenna parameters,
            or None to create an empty antenna.
        """
        self.antenna = {}
        if source is None:
            return
        elif isinstance(source, (str, os.PathLike)):
            self.antenna = self.load(source)
        elif isinstance(source, dict):
            self.antenna = source
        else:
            raise TypeError("source must be a string, PathLike, dictionary, or None")

    @classmethod
    def load(cls, filename: str | os.PathLike) -> dict:
        """
        Load antenna parameters from a TOML file.

        Parameters
        ----------
        filename : str or Path
            Path to a TOML file containing antenna parameters.

        Returns
        -------
        dict
            Dictionary containing antenna parameters.
        """
        with open(filename, "rb") as fp:
            antenna_data = cls.validate_description(tomllib.load(fp))
        return antenna_data

    @classmethod
    def from_dict(cls, antenna_dict: dict) -> "Antenna":
        """
        Create an Antenna instance from a dictionary.

        Parameters
        ----------
        antenna_dict : dict
            Dictionary containing antenna parameters.

        Returns
        -------
        Antenna
            An instance of the Antenna class.
        """
        _ant = cls()
        _ant.antenna = cls.validate_description(antenna_dict)
        return _ant

    @classmethod
    def from_file(cls, filename: str | os.PathLike) -> "Antenna":
        """
        Create an Antenna instance from a TOML file.

        Parameters
        ----------
        filename : str or Path
            Path to a TOML file containing antenna parameters.

        Returns
        -------
        Antenna
            An instance of the Antenna class.
        """
        antenna_data = cls.load(filename)
        return cls.from_dict(antenna_data)

    def __str__(self):
        return f"Antenna(name={self.antenna['name']}, frequency={self.antenna['frequency']})"

    @classmethod
    def validate_description(cls, antenna: dict) -> dict:
        """
        Check that the antenna description is valid.

        Parameters
        ----------
        antenna : dict
            ALOHA antenna description dictionary.

        Returns
        -------
        antenna : dict
            Verified antenna description
        """
        if not antenna:  # empty dict is OK
            return antenna

        total_nb_modules = antenna["global"]["nb_mod_phi"] * antenna["global"]["nb_mod_theta"]
        # excitations
        if len(antenna["excitation"]["magnitudes"]) != total_nb_modules:
            raise ValueError(
                f"Number of excitation magnitudes ({len(antenna['excitation']['magnitudes'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )
        if len(antenna["excitation"]["phases"]) != total_nb_modules:
            raise ValueError(
                f"Number of excitation phases ({len(antenna['excitation']['phases'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )

        # module indices
        if len(antenna["global"]["idx_mod"]) != total_nb_modules:
            raise ValueError(
                f"Number of module indices ({len(antenna['global']['idx_mod'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )

        # S-parameter
        if len(antenna["sparameters"]["filenames"]) != total_nb_modules:
            raise ValueError(
                f"Number of S-parameter filenames ({len(antenna['sparameters']['filenames'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )
        if len(antenna["sparameters"]["phases_deembedded"]) != total_nb_modules:
            raise ValueError(
                f"Number of deembedded phases ({len(antenna['sparameters']['phases_deembedded'])}) "
                f"does not match total number of modules ({total_nb_modules})"
            )
        # so far so good
        return antenna

    def is_valid(self) -> bool:
        """
        Check if the antenna description is valid.

        Returns
        -------
        state : bool
            True is the antenna description looks valid.
        """
        try:
            self.validate_description(self.antenna)
            return True
        except ValueError:
            return False

    def plot(self):
        """
        Plot the antenna architecture.
        """
        if not self.is_valid():
            raise ValueError("Invalid antenna description")

        b, a, z, y, nwr, nwas, act_module_tor = self.antenna_coordinates()
        # Extract parameters from antenna dictionary
        antenna = self.antenna
        plt.figure()
        # TODO : plot all poloidal modules if more than one
        for idx_pol in range(len(y)):
            for idx_tor in range(len(z)):
                rect_pos = [z[idx_tor], y[idx_pol], b[idx_tor], a[idx_pol]]
                # Create passive/active mask
                ar_modules = np.ones(antenna["global"]["nb_mod_phi"])
                ar_pa_mask = np.array(antenna["module"]["mask"])

                # Add passive waveguides between modules
                ar_pa_mask = np.concatenate([ar_pa_mask, np.zeros(antenna["module"]["nb_pwg_btw_mod_phi"])])

                pa_mask = np.kron(ar_modules, ar_pa_mask)

                # Remove last element if there are passive waveguides between modules
                if antenna["module"]["nb_pwg_btw_mod_phi"] > 0:
                    pa_mask = pa_mask[:-1]

                # Add passive waveguides at edges
                pa_mask = np.concatenate(
                    [np.zeros(antenna["module"]["nb_pwg_edge"]), pa_mask, np.zeros(antenna["module"]["nb_pwg_edge"])]
                )

                if pa_mask[idx_tor] == 0:  # passive wg
                    plt.gca().add_patch(
                        plt.Rectangle(
                            (rect_pos[0], rect_pos[1]),
                            rect_pos[2],
                            rect_pos[3],
                            facecolor=[0.8, 0.8, 0.8],
                            edgecolor="k",
                        )
                    )
                elif pa_mask[idx_tor] == 1:  # active wg
                    plt.gca().add_patch(
                        plt.Rectangle(
                            (rect_pos[0], rect_pos[1]), rect_pos[2], rect_pos[3], fill=False, facecolor=[0.8, 0.8, 0.8]
                        )
                    )

        plt.axis("equal")
        plt.xlabel("z [m]")
        plt.ylabel("y [m]")
        plt.title(f"ALOHA antenna: {antenna['name']}\n (view from the plasma)")
        plt.show()

    def antenna_coordinates(self):
        """
        Generate waveguide positions and dimensions from the antenna description.

        Returns
        -------
        b : list
            Waveguide widths (in toroidal direction).
        a : list
            Waveguide heights (in poloidal direction).
        z : list
            Waveguide toroidal positions.
        y : list
            Waveguide poloidal positions.
        nb_wg_per_row : int
            Number of waveguides per row.
        nb_wg_total : int
            Total number of waveguides per antenna.
        act_module_tor : list
            Mask for passive/active waveguides.

        """
        antenna = self.antenna
        mod = antenna["global"]
        wg = antenna["module"]

        # (total) number of waveguides per row
        # = (nb wg in a module) + 2*(nb ext wg) + (nb of wg between modules)
        nb_wg_per_row = (
            mod["nb_mod_phi"] * wg["nb_wg_phi"]
            + 2 * wg["nb_pwg_edge"]
            + (mod["nb_mod_phi"] - 1) * wg["nb_pwg_btw_mod_phi"]
        )

        # (total) number of waveguides per column
        nb_wg_per_col = mod["nb_mod_theta"] * wg["nb_wg_theta"]

        # total number of waveguides
        nb_wg_total = nb_wg_per_row * nb_wg_per_col

        # waveguide height - supposed constant for all the waveguides of the antenna
        a = np.full(nb_wg_total, wg["wg_size_theta"])

        # b
        # Make the array b which contains all the waveguide width of a row of waveguides
        b_module = np.where(wg["mask"], wg["awg_size_phi"], wg["pwg_size_phi"])
        b_edge = np.full(wg["nb_pwg_edge"], wg["pwg_size_edge_phi"])
        b_inter = np.full(wg["nb_pwg_btw_mod_phi"], wg["pwg_size_phi"])

        b = np.concatenate(
            [b_edge, np.tile(np.concatenate([b_module, b_inter]), mod["nb_mod_phi"] - 1), b_module, b_edge]
        )

        # e
        # Make the array e which contains all the waveguide septum width of a row of waveguides
        ne_phi = (
            wg["nb_pwg_btw_mod_phi"] * (mod["nb_mod_phi"] - 1)
            + wg["nb_pwg_edge"] * 2
            + mod["nb_mod_phi"] * wg["nb_wg_phi"]
            - 1
        )
        e = np.tile(wg["pwg_size_phi"], ne_phi)

        # z
        # Make the array z which contains all the waveguide positions in the toroidal direction
        z = np.zeros(nb_wg_per_row)
        for ind in range(1, nb_wg_per_row):
            z[ind] = z[ind - 1] + b[ind - 1] + e[ind - 1]

        # y
        # Make the array y which contains all the waveguide positions in the poloidal direction
        h = np.concatenate(
            [
                np.tile(
                    np.concatenate(
                        [np.full(wg["nb_wg_theta"], wg["awg_size_phi"]), np.full(1, mod["spacing_btw_mod_phi"])]
                    ),
                    mod["nb_mod_theta"] - 1,
                ),
                np.full(wg["nb_wg_theta"], wg["awg_size_phi"]),
            ]
        )
        y = np.zeros(nb_wg_per_col)
        for ind in range(1, nb_wg_per_col):
            y[ind] = y[ind - 1] + h[ind - 1] + a[ind]

        # index of active waveguides in a module
        act_module_tor = wg["mask"]

        return b, a, z, y, nb_wg_per_row, nb_wg_total, act_module_tor
