function aloha_scenario_conversion4Annika(scenario, filename)
% Convert a scenario result into a .mat field, for Vladimir Fuch
% 
% INPUTS:
%   scenario
%   filename
% OUTPUT : none
% 
% 

dP_nz = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
nz   = squeeze(aloha_scenario_get(scenario, 'nz'));

save(filename, 'dP_nz', 'nz', '-V4');