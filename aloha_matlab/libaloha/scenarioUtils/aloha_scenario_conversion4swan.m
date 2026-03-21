function swan_input=aloha_scenario_conversion4swan(scenario)
% ALOHA
% 
% INPUT
%  - scenario : ALOHA scenario
% 
% OUTPUT
%  - swan_input : string corresponding to the input format of SWAN
%  SWAN may be run by command line 'swan <swan_input >swan_output'
% 
% Export an ALOHA scenario in order to create the SWAN input file
% 
% AUTHOR: JH
% LAST UPDATE:
%  - 09/06/2010: creation

if exist(scenario.antenna.architecture)
    eval(scenario.antenna.architecture);    
else
    error(['Architecture : ', scenario.antenna.architecture,' doesn''t found! Check for typo or matlab PATH']);
end

%% parameters conversion ALOHA -> SWAN
filename = 'TALOHA'; % 6 char only !
freq = scenario.antenna.freq/1e9;% in GHz

dne0 = scenario.plasma.ne0./scenario.plasma.lambda_n(1);

ne0 = scenario.plasma.ne0*(1e-2)^3 % m^-3 -> cm^-3
dne = scenario.plasma.dne(1)*(1e-2)^4% m^-4 -> cm^-4

[b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(scenario.antenna.architecture);
option.gap = 1; % is there a gap in front of the antenna ?
option.nbWG = nbre_guides; % number of secundary waveguides
option.nbTM = scenario.options.modes(2); % number of evanescent TM modes
ee = repmat(e, 1, option.nbWG);
wg_dim = [b', ee']*(1e+2);% m -> cm

option.general(1) = 1; % are numerical results associated to spectrum printed ?
option.general(2) = 1; % is spectrum plotted ? 
option.general(3) = 1; % is electric field plotted ?
option.general(4) = -1; % are experimental values of electric field entered ?
option.general(5) = 1; % is spectrum stored for ray tracing calculation ?

% max n// value (defaut:8) and number of point(defaut:500)
option.n_par = '  8.00500 -1';

option.nbPWG = 0; % number of passive waveguides [TODO]
option.nbMJtype = 1; % number of multijunction type used [TODO]
option.nbMJ = option.nbWG;% number of multijunction unit [TODO]
option.nbFC = 1;% number of feeding case [TODO]

option.MJtype = repmat(1,1,option.nbMJ); % indicates that the MJ are of the same type

option.nbSWGpMJ = 1; % number of secundary WG in one MJ
option.SMJknow = -1; % precise if intrinsic scattering matrix of a MJ is known

%% generate SWAN input file
swan_input=[filename, sprintf(' %0.3e %0.3e %0.3e  %i %i %i', [freq, ne0, dne, option.gap, option.nbWG, option.nbTM])]; % line#1
swan_input=strvcat(swan_input, num2str(wg_dim, ' %0.3e %0.3e' )); % waveguide dimensions
swan_input=strvcat(swan_input, ' 1.000e+00 0.000e+00 0.000E-00'); % vacuum gap set to 0.0 (3rd nb)
swan_input=strvcat(swan_input, sprintf('  %i  %i  %i  %i  %i',option.general));
swan_input=strvcat(swan_input, sprintf('%s',option.n_par));
swan_input=strvcat(swan_input, ' 0.00  0  1'); % graphical output [obsolete]
swan_input=strvcat(swan_input, ' '); % blank line
swan_input=strvcat(swan_input, sprintf(' %i  %i  %i  %i',option.nbPWG, option.nbMJtype, option.nbMJ, option.nbFC));
% passive waveguide section [TODO]
swan_input=strvcat(swan_input, sprintf(' %i ', option.MJtype));
swan_input=strvcat(swan_input, sprintf(' %i  %i',option.nbSWGpMJ, option.SMJknow));
% scattering matrix [TODO]
% here : S matrix of a perfect WG
swan_input=strvcat(swan_input, '   0.000   0.000');
swan_input=strvcat(swan_input, '   1.000   0.000');
swan_input=strvcat(swan_input, '   0.000   0.000');   
% excitation
swan_input=strvcat(swan_input, num2str([scenario.antenna.a_ampl, scenario.antenna.a_phase*180/pi], '  %0.3e %0.3e'));
%  swan_input=strvcat(swan_input, 'test')