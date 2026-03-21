
function scenario12 = aloha_scenario_append(scenario1, scenario2)
% Append a scenario structure to an other
%  
% EXAMPLE: 
% scenario12 = aloha_scenario_append(scenario1, scenario2)
%
% 
% INPUT:
%  - scenario1 [structure scenario]
%  - scenario2 [structure scenario]
%  
% OUTPUT:
%  - scenario12 [structure scenario]
%  
% 
% AUTHOR: JH
% LAST CHANGES:
%  - 8/8/2008 (Beijing Olympic Game): creation
% 

    scenario12 = [scenario1; scenario2];
