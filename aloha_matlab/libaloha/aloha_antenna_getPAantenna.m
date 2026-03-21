function antenna_lh = aloha_antenna_getPAantenna(freq, nb_modules, b_a, b_p, e)
% ALOHA
% 
% antenna_lh = aloha_antenna_getPAantenna(Nb_modules) 
% 
% Generate a Passive-Active simple antenna, made of Nb_modules: 
% * Each module consists in a single active waveguide.
% * A passive waveguide is inserted between each module (ie.between each wg)
% * A passive waveguide is inserted at both side of the antenna
% 
% NB: the excitation of the antenna is not set in this function: 
% it has to be made after calling this function !
% (ie: set the antenna_lh.setup.modules.amplitude and phase fields) 
% 
% INPUT
%  - Nb_modules : number of modules (active waveguides)
% 
% OUTPUT
%  - antenna_lh <structure>: ITM CPO ALOHA LH antenna description
%
% AUTHOR: JH
% June 2011

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
antenna_lh.name = ['One row of passive and active waveguides.', num2str(nb_modules), ' actives wg.'];

%  Frequency [Hz]
antenna_lh.frequency = freq;

%  Power [W]
antenna_lh.power  = []; % not defined here in ALOHA

% Main parallel refractive index of the launched spectrum, for multi-junction antennas. 
antenna_lh.n_par = []; % depends on the dimensions of the wg...

%% ----------------------------------------
%% Modules description
%% ----------------------------------------
% Number of modules per antenna in the poloidal direction. 
modules.nma_theta = 1;
%  Number of modules per antenna in the toroidal direction. 
modules.nma_phi = nb_modules;

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
waveguides.nwm_theta = 1;

% Number of waveguides per module in the toroidal direction. (passive and active)
waveguides.nwm_phi = 1;

% Mask of passive and active waveguides for an internal module
% 1 for active -- 0 for passive.
waveguides.mask = [1];

% Number of passive waveguide between modules in the toroidal direction. 
waveguides.npwbm_phi = 1;

% Number of passive waveguides on each antenna edge in the toroidal direction. 
waveguides.npwe_phi = 1;

% Spacing between poloidally neighboring waveguides [m]
waveguides.sw_theta = 12e-3;

% Height of waveguides in the poloidal direction [m]
waveguides.hw_theta = 70e-3;

% Width of active waveguides [m]
waveguides.bwa = b_a;     

% Width of internal passive waveguides [m]
waveguides.biwp = b_p;

% Width of edge passive waveguides [m]
waveguides.bewp = b_p;

% Thickness between waveguides in the toroidal direction [m]
% Reminder : length(e_phi) = nma_phi*nwm_phi + (nma_phi - 1)*npwbm_phi + 2*npwe_phi - 1
ep = e;
ne_phi = waveguides.npwbm_phi*(modules.nma_phi-1) + ...
	 waveguides.npwe_phi*2 + ...
	 modules.nma_phi*waveguides.nwm_phi - 1;
waveguides.e_phi = repmat(ep, 1, ne_phi);

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
modules.Sparameters.pathTo = 'S_HFSS/matrices_HFSS_elem';  % path the S-matrix files (.m files) 

filenames = ['S_elem'];

for ind_mod_fich = [2:modules.nma_phi]
    filenames = [filenames;'S_elem'];
end
modules.Sparameters.SFileNames = filenames;

%% Phase deembedding
% These parameters are the phase correction in order to take into account
% the transmission line length between phase measurement and S-matrix description.
% This is only usefull when using input data from experiments.
modules.Sparameters.phase_deembedded = zeros(modules.nma_phi,1);



%% --------------------------------
%% Other antenna_lh CPO parameters
%% --------------------------------
%% Not defined here in ALOHA
%  Plasma edge characteristics in front of the antenna.
antenna_lh.plasmaedge = [];

%  Amplitude of the TE10 mode injected in the module [W], Matrix (nantenna_lh,max_nmodules). Time-dependent
%  modules.amplitude = zeros(1, modules.nma_theta*modules.nma_phi)';
modules.amplitude = sqrt(1/modules.nma_phi)*ones(modules.nma_phi,1); % 1 W 

%  Phase of the TE10 mode injected in the module [rd], Matrix (nantenna_lh, max_nmodules). Time-dependent
%  modules.phase = zeros(1, modules.nma_theta*modules.nma_phi)';
modules.phase = (270*pi/180)*(0:modules.nma_phi-1)'; % positive n// peak

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

% % Make the array b which contains all the waveguide width 
% % of a row of waveguides
% % Not mandory for ITM CPO antenna_lh, but usefull for ALOHA
% b_module = waveguides.mask.*waveguides.bwa + not(waveguides.mask).*waveguides.biwp; % waveguide width inside a module
% b_edge = repmat(waveguides.bewp, 1, waveguides.npwe_phi);  % passive wg width on each side
% b_inter= repmat(waveguides.biwp, 1, waveguides.npwbm_phi); % passive wg width between modules
%
% waveguides.b = [b_edge, kron(ones(1,modules.nma_phi-1),[b_module, b_inter]),b_module, b_edge];

%  Detailed description of LH antennas.
modules.waveguides = waveguides;
setup.modules = modules;
antenna_lh.setup = setup;



