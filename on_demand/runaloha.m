function runaloha(pulsenb, t_start, t_end, varargin)
% ALOHA ON-DEMAND
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

% Generate the input ALOHA scenario from TS database 
sc = aloha_ondemand_scenario(pulsenb, t_start, t_end, varargin);

% run ALOHA
sc = aloha_scenario(sc);

% save results into an ALOHA scenario file
fileName = ['ALOHA_TS', num2str(pulsenb), '_', ...
            num2str(round(t_start)), '_', num2str(round(t_end)), '_', ...
            sc.options.TSport, '_', sc.antenna.architecture];
aloha_scenario_save(sc, fileName);

% convert and send ALOHA results into the TS database
% save results to be imported into the DB into a matlab binary file
disp('Write into a file compatible to the Tore Supra database');
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