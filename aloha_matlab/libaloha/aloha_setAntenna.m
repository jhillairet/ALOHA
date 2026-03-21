function scenario = aloha_setAntenna(scenario, antenna_name)
% ALOHA
% 
% Defines the ITM antenna_lh CPO for a scenario (or array of scenarios).
%
% INPUT
%  - scenario <structure>: ALOHA scenario(s) 
%  - antenna_name [string] : antenna name (ITM antenna description)
%
% OUTPUT
%  - scenario <structure>: ALOHA scenario(s) with the antenna defined in it.
%
% AUTHOR: JH
% LAST UPDATE
% - 01/12/2011. Corrected a bug on the matlab PATH and add support to array of scenarios
% 

% If the 'architecture_antenne' directory is not in the Matlab PATH, then add it
if not(exist('architecture_antenne'))
    addpath(genpath([aloha_utils_getRootPath, '/architecture_antenne']));
end


% In case of array of scenarios 
for idscen=1:length(scenario)

% Load antenna if the antenna description 
% comply to the ITM-way of definition, ie if is is a function
% (in contrast to a simple script)
% Otherwise, this is a old fashioned antenna description
% or the file does not exist !

try
   % Try to load the antenna description as a Matlab function
   % The antenna description is a Matlab function only with ITM antenna description 
   scenario(idscen).antenna_lh = eval(antenna_name);
catch ME
   % So if the previous command has failed, this may be because
   % the antenna description is not a Matlab function, but a classic script, 
   % as it is for old fashioned antenna: this means that it is a old antenna description
    if strcmp(ME.identifier, 'MATLAB:scriptNotAFunction')
        disp(aloha_message('Old fashioned antenna description... passing...'));
        scenario(idscen).antenna_lh = [];
        scenario(idscen).antenna_lh.setup = [];
    elseif strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        % or simply because the file has not been found (in the PATH)
        error('The antenna description has not been found in the Matlab PATH. Sure it exists?');
    else
        % or I don't know why !?
        error('Unknown error!!')
   end 
end

end % for
