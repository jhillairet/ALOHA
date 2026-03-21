function antenna_lh = aloha_getAntenna(scenario)
% ALOHA
% 
% Returnes the ITM antenna_lh CPO from a given scenario.
%
% INPUT
%  - scenario <structure>: ALOHA scenario
%
% OUTPUT
%  -  antenna_lh <structure> : Matlab antenna_lh CPO
%
% AUTHOR: JH

% TODO : fix for array of scenarios
antenna_lh = scenario.antenna_lh;
