function aloha_scenario_conversion4vladimir(scenario, filename)
% Convert a scenario result into a .mat field, for Vladimir Fuch
% 
% INPUTS:
%   scenario
%   filename
% OUTPUT : none
% 
% 

dP_z = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
nz   = squeeze(aloha_scenario_get(scenario, 'nz'));
z    = squeeze(aloha_scenario_get(scenario, 'abs_z'));
Ez   = squeeze(aloha_scenario_get(scenario, 'Efield'));

% if the field into the plasma has been computed
% add it to the output structure
if isfield(scenario.results, 'Ez_x_z')
    Ez_x_z = scenario.results.Ez_x_z;
    Hy_x_z = scenario.results.Hy_x_z;
    x_coord = scenario.results.x_coord;
    z_coord = scenario.results.z_coord;
    nbre_x_coord = scenario.options.nbre_x_coord;
    nbre_z_coord = scenario.options.nbre_z_coord;
    save(filename, 'dP_z', 'nz', 'z', 'Ez', ...
        'Ez_x_z', 'Hy_x_z', 'x_coord', 'z_coord', 'nbre_x_coord', 'nbre_z_coord');
else
    save(filename, 'dP_z', 'nz', 'z', 'Ez','-v7.3');
end