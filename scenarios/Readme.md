The following schema conversion is used from the
deprecated Matlab scenarios to the new Python implementation.

matlab structure variable -> Python nested dictionary variable

if "-> None" means matlab variable is deprecated and no more used in the Python dict
if "None ->" means matlab variable didn't exist. A default values is indicated for the Python dict

# General
options.bool_save -> None
options.scenario_filename -> None
options.comment -> comment

# Antenna description
antenna.architecture -> antenna.file
antenna.freq -> antenna.excitation.f
options.bool_mesure -> antenna.excitation.experimental
options.bool_homeMadeExcitation -> None
antenna.a_ampl -> antenna.excitation.power
antenna.a_phase -> antenna.excitation.power (in degree)
options.TSport -> antenna.excitation.port
options.choc -> antenna.excitation.pulse
options.tps_1 -> antenna.excitation.avg_times[0]
options.tps_2 -> antenna.excitation.avg_times[1]

# Plasma description
None -> plasma.solver (default: "spectral_1D")
options.version_code -> None
version_plasma_1D -> plasma.spectral_1D.profile (3 -> 'linear; 6 -> 'bilinear')
Nmh -> None
Nme -> plasma.spectral_1D.nb_evanescent_modes

options.bool_lignes_identiques -> plasma.spectral_1D.bilinear.identical_profiles
plasma.ne0 -> plasma.spectral_1D.bilinear.ne0
plasma.lambda_n(1) -> plasma.spectral_1D.bilinear.lambda_n[0]
plasma.lambda_n(2) -> plasma.spectral_1D.bilinear.lambda_n[1]
plasma.d_couche -> plasma.spectral_1D.bilinear.plasma_layer_length
options.B0 -> plasma.spectral_1D.bilinear.B0
plasma.d_vide -> plasma.spectral_1D.bilinear.vacuum_layer_length
options.type_swan_aloha -> plasma.spectral_1D.bilinear.infinite_waveguide (0 -> False, if 1 -> True)

# Options

options.bool_display_density_profile -> None
options.bool_compute_spectrum -> None
options.bool_display_spectrum -> None
options.bool_compute_directivity -> None
options.bool_display_directivity -> None
options.definition_directivite -> None
options.bool_compute_total_field -> None
options.bool_display_total_field -> None
options.bool_compute_plasma_field -> None
options.bool_display_plasma_field -> None
options.aloha_path -> None
options.bool_debug -> options.debug

plasma.dne(:,1) -> None
plasma.dne(:,2) -> None
options.modes -> None

options.nz_min -> plasma.spectral_1D.bilinear.nz_min
options.nz_max -> plasma.spectral_1D.bilinear.nz_max
options.dnz -> plasma.spectral_1D.bilinear.dnz
options.ny_min -> plasma.spectral_1D.bilinear.ny_min
options.ny_max -> plasma.spectral_1D.bilinear.nz_max
options.dny -> plasma.spectral_1D.bilinear.dny

options.ny_nb -> None
options.nz_nb -> None

options.z_coord_min -> plasma.spectral_1D.bilinear.z_min
options.z_coord_max -> plasma.spectral_1D.bilinear.z_max
options.nbre_z_coord -> plasma.spectral_1D.bilinear.nb_z
options.x_coord_max -> plasma.spectral_1D.bilinear.x_max
options.nbre_x_coord -> plasma.spectral_1D.bilinear.nb_x

options.pas_nz_fig_plasma -> None
options.fig_Ez_ou_EzHy -> None
options.lig_fig_plasma -> None
results.S_plasma -> None
antenna_lh.setup -> None
