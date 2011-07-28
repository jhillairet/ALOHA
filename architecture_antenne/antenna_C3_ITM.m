function antenna_lh = antenna_ITM
% ALOHA description of an LH antenna
% 
% Tore Supra C3 antenna (half antenna)
% ITM CPO antenn_lh description
% 
% INPUT : none
%  
% OUTPUT
%  - antenna_lh <structure> : Matlab CPO antenna_lh
% 
% 

% This file contains the description of an LH antenna according to the ITM
% defition, i.e. using the CPO antenna_lh format. 
% 
% An LH antenna is defined as the following :
%  - general description
%  - modules description
%  - waveguides description (dimensions, etc.)
% All the parameters are described in the (poloidal,toroidal)=(theta,phi) frame.

%% --------------------------------------
%% General description of the antenna
%% --------------------------------------
%  Antenna name
antenna_lh.name = 'Tore Supra C3 (half)';

%  Frequency [Hz]
antenna_lh.frequency = 3.7e9;

%  Power [W]
antenna_lh.power  = []; % not defined here in ALOHA

% Main parallel refractive index of the launched spectrum, for multi-junction antennas. 
antenna_lh.n_par = 2.02;


%% ----------------------------------------
%% Modules description
%% ----------------------------------------
% Number of modules per antenna in the poloidal direction. 
modules.nma_theta = 1;
%  Number of modules per antenna in the toroidal direction. 
modules.nma_phi = 8;

%  Position index of the module in the poloidal direction (from low theta to high theta, 
%  i.e. from bottom to top if the antenna is on LFS). 
modules.ima_theta = ones(1,modules.nma_phi*modules.nma_theta); % numbering in ALOHA goes from top to bottom as view from the plasma.

%  Position index of the module in the toroidal direction (from low phi to high phi, 
%  counter-clockwise when seen from above).
modules.ima_phi = [modules.nma_phi:-1:1]; % numbering in ALOHA goes from left to right as view from the plasma.

%  Spacing between poloidally neighboring modules [m]
modules.sm_theta = 0;

%% ----------------------------------------
%% Waveguides description
%% ----------------------------------------
%  Number of waveguides per module in the poloidal direction. (passive and active)
waveguides.nwm_theta = 3;

% Number of waveguides per module in the toroidal direction. (passive and active)
waveguides.nwm_phi = 6;

% Mask of passive and active waveguides for an internal module
% 1 for active -- 0 for passive.
waveguides.mask = [1 1 1 1 1 1];

% Number of passive waveguide between modules in the toroidal direction. 
waveguides.npwbm_phi = 1;

% Number of passive waveguides on each antenna edge in the toroidal direction. 
waveguides.npwe_phi = 1;

% Spacing between poloidally neighboring waveguides [m]
waveguides.sw_theta = 12e-3;

% Height of waveguides in the poloidal direction [m]
waveguides.hw_theta = 70e-3;

% Width of active waveguides [m]
waveguides.bwa = 8e-3;     

% Width of internal passive waveguides [m]
waveguides.biwp = 6.5e-3;

% Width of edge passive waveguides [m]
waveguides.bewp = 6.5e-3;

% Thickness between waveguides in the toroidal direction [m]
% Reminder : length(e_phi) = nma_phi*nwm_phi + (nma_phi - 1)*npwbm_phi + 2*npwe_phi - 1
e = 2e-3; % between active waveguides
e_pwg = 3e-3; % between passive waveguides

e_aw = [e e e e e];
e_mod = [e_pwg e_aw e_pwg];

waveguides.e_phi = repmat(e_mod, 1, 8);

% Short circuit length for passive waveguides [m]
% Reminder : length(scl) = nma_phi*npwm_phi + (nma_phi - 1)*npwbm_phi + 2*npwe_phi
nscl = waveguides.npwbm_phi*(modules.nma_phi-1) + ...
       waveguides.npwe_phi*2 + ...
       sum(not(waveguides.mask))*modules.nma_phi;
waveguides.scl = repmat(1/4, 1, nscl);     

%% --------------------------------
%% Modules Scattering parameters
%% --------------------------------
% matrice S des modules ds des fichiers .m (NB : la matrice est rangee sur une seule colonne)
modules.Sparameters.pathFrom = pwd;
modules.Sparameters.pathTo = 'S_HFSS/matrices_HFSS_C3';  % path the S-matrix files (.m files) 
% modules C3
% La modelisation HFSS de l'antenne dÃ©bute au voisinage du CM.
% Pour prendre en compte le dephasage entre la mesure de phase et le debut de la modelisation HFSS
% il faut calculer le dephasage necessaire : 
%  phase incident + phase rallonge + phase s12 HFSS = phase en bout
% soit :
%  phase rallonge = dphi - phase s12 HFSS
% ou dphi est la correction (connue, cf. doc Annika) utilisee pour calculer la phase en bout
% a partir de la phase incidente.


% modules HAUT
% -----------
% module    correction dphi     phase S12 HFSS    
% 24h (= 8H)    17               168.6118
% 23h           11              -154.7044
% 22h           1               -129.9639
% 21h           15              -118.7455
% 14h           16              -113.2875
% 13h           22              -117.8509
% 12h           17              -132.5663
% 11h (= 1B)    15              -164.3341        
%  
%  modules.Sparameters.SFileNames = ['S_C3_24h';'S_C3_23h';'S_C3_22h';'S_C3_21h';'S_C3_14h';'S_C3_13h';'S_C3_12h';'S_C3_11h'];
%  phase_rallonge = (pi/180)*[17-168.6; 11+154.7; 1+130; 15+118.7; 16+113.3; 22+117.9; 17+132.6; 15+164.3];

modules.Sparameters.SFileNames = ['S_C3_24b';'S_C3_23b';'S_C3_22b';'S_C3_21b';'S_C3_14b';'S_C3_13b';'S_C3_12b';'S_C3_11b'];


%% Phase deembedding
% These parameters are the phase correction in order to take into account
% the transmission line length between phase measurement and S-matrix description.
% This is only usefull when using input data from experiments.
% modules.Sparameters.phase_deembedded = zeros(nma_phi,1);
modules.Sparameters.phase_deembedded = (pi/180)*[0-166.44; -5.5+160.9; 8+135.3; 14+122; 8.5+116.1; 20+131.2; 10+139.8; 20+166.9];

%% --------------------------------
%% Other antenna_lh CPO parameters
%% --------------------------------
%% Not defined here in ALOHA
%  Plasma edge characteristics in front of the antenna.
antenna_lh.plasmaedge = [];

%  Amplitude of the TE10 mode injected in the module [W], Matrix (nantenna_lh,max_nmodules). Time-dependent
modules.amplitude = sqrt(4.0320e6/16)*ones(modules.nma_phi,1);

%  Phase of the TE10 mode injected in the module [rd], Matrix (nantenna_lh, max_nmodules). Time-dependent
modules.phase = -90*(pi/180)*(0:modules.nma_phi-1)';

%% Not used at all in ALOHA - 
%  Reference global antenna position. Vectors (nantenna_lh). Time-dependent
antenna_lh.position = []; 
%  Beam characteristics
antenna_lh.beam = [];


%% ------------------------------  
%% DO NOT EDIT UNDER THIS LINE
%% ------------------------------
% architecture name of the current antenna = its filename
antenna_lh.archName = mfilename;

% Make the array b which contains all the waveguide width 
% of a row of waveguides
% Not mandory for ITM CPO antenna_lh, but usefull for ALOHA
b_module = waveguides.mask.*waveguides.bwa + not(waveguides.mask).*waveguides.biwp; % waveguide width inside a module
b_edge = repmat(waveguides.bewp, 1, waveguides.npwe_phi);  % passive wg width on each side
b_inter= repmat(waveguides.biwp, 1, waveguides.npwbm_phi); % passive wg width between modules

waveguides.b = [b_edge, kron(ones(1,modules.nma_phi-1),[b_module, b_inter]),b_module, b_edge];

%  Detailed description of LH antennas.
modules.waveguides = waveguides;
setup.modules = modules;
antenna_lh.setup = setup;



