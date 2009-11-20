function scenario = aloha_scenario_load(varargin)
% Load a scenario structure(s) from a matlab file .mat
% 
% scenario = aloha_scenario_load(fileName)
% Or,
% scenario = aloha_scenario_load(fileName1, fileName2, ...) 
%    
% INPUT
%  - fileName [string]
%  Or,
%  - fileName1 [string] 
%  - fileName2 [string]
%  - ... 
%  
%  
% OUTPUT
%  - scenario [structure(s) scenario] loaded from fileName 
%  Or 
%  - scenario [structure(s) scenario] concatened loaded scenarios from fileName1, fileName2, ...
% 
% NB: if the file doesn't exist or is unreadable, a empty variable is returned.
%  
% AUHTOR: JH
% LAST CHANGE:
%  - 13/10/2008: many filename are supported : the resulting scenario is concatened
%  - 08/08/2008: creation 
%  


scenario = [];

for idx=1:nargin
    fileName = varargin{idx};

    if exist(fileName, 'file')
        results = load(fileName);
        scenario = aloha_scenario_append(scenario, results.scenario);
    end
end