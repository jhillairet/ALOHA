function sc = aloha_ondemand_scenario(pulsenb, t_start, t_end, varargin)
% ALOHA ON-DEMAND SCENARIO GENERATOR
% 
% This script generates an ALOHA scenario to be run with aloha_scenario() fror a specific TS pulse
%  
% INPUT
%   - pulsenb : Tore Supra pulse number
%   - t_start : start time (s)
%   - t_end : end time (s)
%   
%  IMPORTANT : Other optionnal input arguments :
%  
%   - [optionnal] ne : edge electron density [m^-3]
%   - [optionnal] lambda : electron edge scrape-off length (lambda = ne/grad_ne) [m] :
%                   1 value for one gradient
%                   OR
%                   2 value in an array [lambda0,lambda1] for 2 gradients
%   - [optionnal] TSport : Tore-Supra port 'Q6A' or 'Q6B'
% 
% OUTPUT
%  ALOHA scenario (without results)
%  
% AUTHOR: J.Hillairet
% LAST UPDATE:
%  - 09/10/2009 : creation

% We use the defaut ALOHA scenario example and we adapt it after
% We search for the scenario_example.m on the Matlab PATH. If it not exist
% we look into the scenario directory. This assumes that this directory exist...
if exist('scenario_example') == 2
    sc=scenario_example;
else 
    actual_dir=pwd;
    cd([aloha_utils_getRootPath,'/scenario']);
    sc=scenario_example;
    cd(actual_dir);
end

% check for optionnal input arguments
if nargin == 3
    disp(['Default density : ne0=', num2str(sc.plasma.ne0/1e17), ' x 1e17 m^-3']);
    disp(['Default scrape-off length : lambda0=', num2str(sc.plasma.lambda_n(1)), ' m']);
elseif nargin == 4
    sc.plasma.ne0 = varargin{1};
    disp(['Default scrape-off length : lambda0=', num2str(sc.plasma.lambda_n(1)), ' m']);
elseif nargin == 5
    sc.plasma.ne0 = varargin{1};
    if length(varargin{2}) ==1
        sc.plasma.version = 3;
        sc.plasma.lambda_n(1) = varargin{2};
    else
        sc.plasma.version = 6;
        sc.plasma.lambda_n = varargin{2};
    end
elseif nargin == 6
    sc.plasma.ne0 = varargin{1};
    if length(varargin{2}) ==1
        sc.plasma.version = 3;
        sc.plasma.lambda_n(1) = varargin{2};
    else
        sc.plasma.version = 6;
        sc.plasma.lambda_n = varargin{2};
    end
    TSport = varargin{3}; 
else
    error('Bad number of input arguments. See help.')
end


% Ask for the antenna port to choose for excitation
if not(exist('TSport'))
    TSport=upper(input('Which port  : ''Q6A'' [C2 before 2009, C3 after] or  ''Q6B'' [C3 before 2009, C4 after] ? : ', 's'));
end
% check user input
if not(strcmp('Q6A', TSport)) && not(strcmp('Q6B', TSport))
    error('Bad port input : should be either Q6A or Q6B only.');
end

% Adaptation of the scenario 
sc.options.TSport = TSport;   
sc.options.choc = pulsenb;
sc.options.tps_1 = t_start;
sc.options.tps_2 = t_end;
sc.options.bool_mesure = true; % activate also the phase correction

% the antenna type depends of the pulse number :
% before #43540, C2 was on port Q6A and C3 on port Q6B
% after  #43540, C3 was on port Q6A and C4 on port Q6B
% #43540 was on end november 2008. 
switch TSport
    case 'Q6A'
        if pulsenb < 43540
            sc.antenna.architecture = 'antenne_C2';
        else
            sc.antenna.architecture = 'antenne_C3';
        end
    case 'Q6B'
        if pulsenb < 43540
            sc.antenna.architecture = 'antenne_C3';
        else
            sc.antenna.architecture = 'antenne_C4';
        end
end

% JH 01/12/2011
% ITM compatibility
% This line avoid a bug, which could lead to creating a scenario working
% with two disctincts antennas... this bug was due to the matlab PATH
% 
% clear the ITM antenna description and reload it with the antenna described upstair
sc = aloha_setAntenna(sc, sc.antenna.architecture);

% In any case, the frequency in the same : 3.7GHz;
sc.antenna.freq = 3.7e9;

% Now we get the feeding information from the TS database
[sc.antenna.a_ampl, sc.antenna.a_phase] = aloha_antenna_excitation(pulsenb, t_start, t_end, TSport);
