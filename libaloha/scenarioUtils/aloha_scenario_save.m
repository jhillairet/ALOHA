function status = aloha_scenario_save(scenario, fileName)
% Save a scenario structure into a matlab file .mat (force V6 format for backward compatibility)
%   
% INPUT:
%  - scenario [structure scenario]
%  - fileName [string]
%  
% OUTPUT:
%  - status [int] : 0 OK ; 1 problem
% 
% NB: for bakward compatibility, results are save in 'V6' mode.
% 
% AUTHOR: JH
% LAST CHANGE:
%  - 03/09/2008: ouput status
%  - 08/08/2008: creation
%  


    matlab_version = ver;
    if str2num(matlab_version(1).Version(1)) > 5
        save(fileName, 'scenario', '-v6');
    else
        save(fileName, 'scenario');
    end

    status = 0;