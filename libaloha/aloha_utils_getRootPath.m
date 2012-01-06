function aloha_path = aloha_utils_getRootPath
% Get the ALOHA root path
%  
% INPUT: none
% 
% OUTPUT
%  - aloha_path [string] : ALOHA root path
%  
% AUTHOR:JH
% LAST UPDATE
%  - 20/10/2008: creation

% Hypothesis : the ALOHA root path is given from the path of 
% the function 'aloha_message' which should be located in the 
% sub-directory 'libaloha'.
% From this path, we get the ALOHA root 
aloha_message_path = which('aloha_message');

% if the function is unkown, 
% user is asked to add the libaloha to the matlab path
if isempty(aloha_message_path)
    error('the function aloha_message is unkown ! This is not normal !');
end

% the function aloha_message is already in the path
% We chop the path in order to recover the ALOHA root path
[pathstr, dummy, dummy] = fileparts(aloha_message_path);
[aloha_path, dummy, dummy] = fileparts(pathstr);
    
