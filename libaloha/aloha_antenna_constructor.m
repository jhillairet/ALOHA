function aloha_antenna_constructor(antenna_name)
%% ALOHA
%
% aloha_antenna_constructor(antenna_name)
%
% ALOHA antenna constructor
% 
% Generate a default ALOHA antenna file, which then has to be modified/filled in order to be run by ALOHA.
% The function generate a matlab script .m which name corresponds to the antenna_name input parameter
% This file can then be used to run ALOHA.
%
% The purpose of this function is to generate a antenna description with match the latest syntax of ALOHA.
%
% INPUT
%  - antenna_name [string] : antenna name. Spaces characters (' ') are automatically replaced by '_'
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
antenna_filename = [regexprep(antenna_name, ' ', '_'),'.m'];

% copy the latest default example ALOHA scenario to the current directory
aloha_path = aloha_utils_getRootPath;
source = [aloha_path,'/architecture_antenne/antenna_example.m'];
destination = [cur_path,'/',antenna_filename];
[status,message,messageid] = copyfile(source, destination);
aloha_message(message);

