function [struct_name] = aloha_scenario_convert2cpo(scenario, name)
% ALOHA function : convert an ALOHA scenario to ITM CPOs matlab structures.
% 
% INPUT
%  - scenario : ALOHA scenario
%  - name ; 'antenna_lh' or 'launchs'
% 
% OUPUT
%  - struct_name : matlab CPO structure :
%        antenna_lh (Lower Hybrid antennas CPO) 
%       or 
%        launchs (RF waves launch conditions CPO)
% 
% J.Hillairet
% 
% 

switch name
    case 'antenna_lh'
    % -------------------------------------
    % creation of the antenna_lh structure
    % -------------------------------------
    antenna_lh.name = scenario.antenna.architecture; % antenna name
    antenna_lh.frequency = scenario.antenna.freq; % frequency [Hz]
    antenna_lh.power = sum(scenario.antenna.a_ampl.^2); % total input power [W]
    antenna_lh.mode = +1; % slow wave
    
    % position
    position.r = []; % major radius [m]
    position.z = []; % altitude [m]
    position.phi = []; % toroidal angle [rad]
    antenna_lh.position = position;
    
    % Modules description
    
    
    % setup
    % TODO : presently, antennas are described internally in ALOHA.
    % This will change in the future release of ALOHA, then the following
    % fields will be completed :
    % module geometric description
    modules.nma_theta = [];
    modules.nma_phi = [];
    modules.sm_theta = [];
    % waveguide geometric description
    waveguides.nwm_theta = [];
    waveguides.nwm_phi = [];
    waveguides.mask = [];
    waveguides.npwbm_phi = [];
    waveguides.npwe_phi = [];
    waveguides.sw_theta = [];
    waveguides.hw_theta = [];
    waveguides.bwa = [];
    waveguides.biwp = [];
    waveguides.bewp = [];
    waveguides.e_phi = [];
    waveguides.scl = [];
    modules.waveguides = waveguides;
    
    % amplitude and phase inputs of each modules
    modules.amplitude = scenario.antenna.a_ampl;
    modules.phase = scenario.antenna.a_phase;
    
    
    setup.modules = modules; %
    antenna_lh.setup = setup;

    
    
    
    % plasma edge
    plasmaedge.nmode = sum(scenario.options.modes); % number of modes used for antenna/plasma coupling
    % TODO : presently, ALOHA describes the plasma in front of the antenna analytically.
    % ALOHA supposes that the plasma is described by a linear increase of the density :
    % ne(x) = ne0 + grad_ne*x
    % where x is the radial direction (x=0 is the antenna's mouth, x>0 the plasma)
    plasmaedge.npoints = 10 ; % nb of points in the distance grid
    plasmaedge.distance = linspace(0,50e-2,plasmaedge.npoints) ; % grid for e- density
    plasmaedge.density = scenario.plasma.ne0 ...
        + scenario.plasma.ne0./scenario.plasma.lambda_n(1).*plasmaedge.distance; % e- density in) front of the antenna [m^-3]
    
    antenna_lh.plasmaedge = plasmaedge;
    
    % beam
    beam.spot = [];
    beam.phaseellipse = [];
    antenna_lh.beam = beam;

    struct_name = antenna_lh;

case 'launchs'
    % -------------------------------------
    % creation of the launchs structure
    % -------------------------------------
    launchs.name = scenario.antenna.architecture;
    launchs.frequency = scenario.antenna.freq;
    launchs.type = 'LH';
    % TODO : fill the launchs structure fields ?
    launchs.mode = [];
    launchs.datainfo = [];
    launchs.position = [];
    launchs.beam = [];
    launchs.codeparam = [];
    launchs.time = [];
    
    if isfield(scenario.results, 'dP')
        launchs.spectrum.nn_phi = scenario.options.nbre_ny;
        launchs.spectrum.nn_theta = scenario.options.nbre_nz;
        launchs.spectrum.n_phi = scenario.results.ny;
        launchs.spectrum.n_theta = scenario.results.nz;
        launchs.spectrum.power = scenario.results.dP;
    else
        disp('No spectrum in the scenario !');
        launchs.spectrum = [];
    end
    struct_name = launchs;

otherwise
    error('bad structure name! ');
end