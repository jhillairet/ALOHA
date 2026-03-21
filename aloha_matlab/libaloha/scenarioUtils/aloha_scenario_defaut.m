function scenario = aloha_scenario_defaut(varargin)
% Create a defaut scenario for the ALOHA code
%  
% scenario = aloha_scenario_defaut(varargin)
%  
% INPUT
%  
% OUPUT
%  - scenario [structure] : ALOHA scenario structure
%  
% AUTHOR:JH
% LAST CHANGES:
%  - 05/09/2008: creation
%  - %  


plasma.d_couche = 2e-2;
plasma.d_vide   = 0;
plasma.ne0      = 2e17;
plasma.lambda_n = [2e-2 2e-2];
plasma.dne      = [1e19 1e19];
plasma.B0       = 2.95;

antenna.a_ampl  = 1;
antenna.a_phase = 0;
antenna.architecture = 'antenne_PAM_TOPLHA';
antenna.freq    = 3.7e9;

options.modes= [1 2];
options.version= 3;

options.comment= 'test';
options.scenario_filename = 'test.mat';

options.aloha_path= '/home/sccp/gchf/JH218595/codes/aloha_v2';

options.dnz= 0.0200;
options.fig_Ez_ou_EzHy= 1;
options.lig_fig_plasma= 1;
options.nbre_x_coord= 8;
options.nbre_z_coord= 60;
options.nz_max= 10;
options.nz_min= -10;
options.pas_nz_fig_plasma= 0.1000;
options.definition_directivite= 1;
options.type_swan_aloha= 1;
options.x_coord_max= 0.0500;
options.z_coord_max= 0.0750;
options.z_coord_min= -0.0150;

options.bool_compute_directivity= 0;
options.bool_compute_plasma_field= 0;
options.bool_compute_spectrum= 1;
options.bool_compute_total_field= 0;
options.bool_debug= 0;
options.bool_display_density_profile= 0;
options.bool_display_directivity= 0;
options.bool_display_plasma_field= 0;
options.bool_display_spectrum= 0;
options.bool_display_total_field= 0;
options.bool_fichier_Ez_x_variable= 0;
options.bool_lignes_identiques= 1;
options.bool_sauvegarde= 1;

options.bool_mesure= 0;
options.choc = 0000;
options.tps_1 = 1;
options.tps_2 = 2;
options.excitation = 'C2';

scenario = aloha_setfield([], plasma, antenna, options);
