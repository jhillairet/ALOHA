function aloha_scenario_constructor(scenario_name)
%% ALOHA
%
% aloha_scenario_constructor(scenario_name)
%
% ALOHA scenario constructor
% 
% Generate a default scenario, which then has to be modified in order to be run by ALOHA.
% The function generate a matlab script .m which name corresponds to the scenario_name input parameter
% This file can then be used to run ALOHA.
%
% The purpose of this function is to generate a scenario witch is the latest scenario syntax of ALOHA.
%
% INPUT
%  - scenario_name [string] : scenario name. Spaces characters (' ') are automatically replaced by '_'
%                             do NOT use accentued characters (such as é, ç, etc...)
% 
% OUTPUT 
%  - none
%  
% AUTHOR: JH
% LAST UPDATE: 20/09/2011
%

% retrieve the current path from which the function is called 
cur_path = pwd;

% generate the filename of the scenario
% replace space char by underscores
scenario_filename = [regexprep(scenario_name, ' ', '_'),'.m'];

% copy the latest default example ALOHA scenario to the current directory
aloha_path = aloha_utils_getRootPath;
source = [aloha_path,'/scenario/scenario_example.m'];
destination = [cur_path,'/',scenario_filename];
[status,message,messageid] = copyfile(source, destination);
aloha_message(message);

