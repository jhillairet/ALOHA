function runaloha(pulsenb, t_start, t_end, varargin)
% ALOHA ON DEMAND
% 
% This script run the LH coupling code ALOHA for a specific TS pulse
% 
% The results of the calculation, such as reflection coefficients of antenna spectra
% are automaticaly saved into a matlab structure file (.mat), which name is 
%   TSXXXXX_YY_ZZ_antenna.mat
% where XXXXX is the Tore Supra pulse number
% YY is the time start, ZZ the time end (without decimal)
% 
% INPUT
%   - pulsenb : Tore Supra pulse number
%   - t_start : start time (s)
%   - t_end : end time (s)
%   
%  Other optinnal input arguments :
%  
%   - [optionnal] ne : edge electron density [m^-3]
%   - [optionnal] lambda : electron edge scrape-off length (lambda = ne/grad_ne) [m]
%   - [optionnal] 
% 
% OUTPUT
%  none
%  
% AUTHOR: J.Hillairet
% LAST UPDATE:
%  - 09/10/2009 : creation

% In order to run ALOHA, we need an input scenario
% This input scenario set up all the options, such
% as the precision dnz of the spectrum.
% 
% We use the defaut ALOHA scenario example and we adapt it after
actual_dir=pwd;
cd('../scenario');
sc=scenario_example;
cd(actual_dir);


% Ask for the antenna port to choose for excitation
if not(exist('TSport'))
    TSport=upper(input('Which port  : ''Q6A'' [C2 before 2009, C3 after] or  ''Q6B'' [C3 before 2009, C4 after] ? : '));
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
sc.options.bool_mesure = true; % activate the phase correction

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

% In any case, the frequency in the same : 3.7GHz;
sc.antenna.freq = 3.7e9;

% Now we get the feeding information from the TS database
[sc.antenna.a_ampl, sc.antenna.a_phase] = aloha_antenna_excitation(pulsenb, t_start, t_end, TSport);

% Asking for plasma properties such as electron density and gradient

% run ALOHA
sc=aloha_scenario(sc);

% save results into an ALOHA scenario file
fileName = ['ALOHA_TS', num2str(pulsenb), '_', ...
            num2str(round(t_start)), '_', num2str(round(t_end)), '_', ...
            sc.options.TSport, '_', sc.antenna.architecture];
aloha_scenario_save(sc, fileName);

% convert and send ALOHA results into the TS database

disp('Write into the Tore Supra database ?');
pause

% save results to be imported into the DB into a matlab binary file
aloha_ntnpar(sc); % TODO

%  fileName = [specdat, int2str(pulsenb)];
%  save fileName ['tdata phcons phmes nparmax effic npar spectre1 spectre2 direc1 direc2 cert scom'];
%  
%  mm='n';
%  mm=input('Ecrire dans la base de donnees (o/n) ? ','s');
%  
%  if ((mm=='o') |  (mm=='O')), ntnpar; end
%  
%  end
%  
%  
%  file