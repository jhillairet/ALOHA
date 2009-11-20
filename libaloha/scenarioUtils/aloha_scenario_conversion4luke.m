function aloha_scenario_conversion4luke(scenario, filename)
% Convert a scenario result into a .mat field, for LUKE
% 
% INPUTS:
%   scenario
%   filename
% OUTPUT : none
% 
% 
antenna.dP = squeeze(aloha_scenario_get(scenario, 'dP'));
antenna.dP_nz = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
antenna.nz = squeeze(aloha_scenario_get(scenario, 'nz'));
antenna.ny = squeeze(aloha_scenario_get(scenario, 'ny'));
antenna.dnz = squeeze(aloha_scenario_get(scenario, 'dnz'));
antenna.dny = squeeze(aloha_scenario_get(scenario, 'dny'));

disp('REDUCTION FACTOR : 4 (blocks) * 12 (rows)');
antenna.reduction_factor = 4*12;

plasma.ne0 = squeeze(aloha_scenario_get(scenario, 'ne0'));
plasma.lambda_n = squeeze(aloha_scenario_get(scenario, 'lambda_n'));


save(filename, 'plasma', 'antenna');