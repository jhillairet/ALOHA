function aloha_batch_constructor(batch_name)
%% ALOHA
%
% aloha_batch_constructor(batch_name)
%
% ALOHA batch constructor
% 
% Generate a default batch, which then has to be modified in order to be run by ALOHA.
% The function generate a matlab script .m which name corresponds to the batch_name input parameter
% This file can then be used to run ALOHA.
%
% The purpose of this function is to generate a batch file witch is the latest batch syntax of ALOHA.
%
% INPUT
%  - batch_name [string] : batch name. Spaces characters (' ') are automatically replaced by '_'
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

% generate the filename of the batch
% replace space char by underscores
batch_filename = [regexprep(batch_name, ' ', '_'),'.m'];

% copy the latest default example ALOHA batch to the current directory
aloha_path = aloha_utils_getRootPath;
source = [aloha_path,'/scenario/batch_example.m'];
destination = [cur_path,'/',batch_filename];
[status,message,messageid] = copyfile(source, destination);
aloha_message(message);

