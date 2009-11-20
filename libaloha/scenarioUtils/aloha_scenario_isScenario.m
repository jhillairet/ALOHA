function res = aloha_scenario_isScenario(scenario)
% Determine si une variable correspond a un scenario ALOHA
%  
%  res = aloha_scenario_isScenario(scenario)
% 
% INPUT
%  - scenario: variable a tester
%
% OUTPUT
%  - res [boolean] : true or false
%  
% AUTHOR:JH
% LAST CHANGES:
%  - 05/09/2008: creation



    res = false;

    if isstruct(scenario)
        if isfield(scenario, 'plasma') & ...
           isfield(scenario, 'antenna') & ...
           isfield(scenario, 'options')
            res = true;
        end
    end