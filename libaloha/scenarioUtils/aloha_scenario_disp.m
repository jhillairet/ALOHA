function aloha_scenario_disp(scenario, level)
% Display the arborescence of an ALOHA scenario structure
% 
% aloha_scenario_disp(scenario, level)
% 
% INPUT
%  - scenario : matlab structure
%  - level [optionnal] : recursion level to dig into
% 
% AUTHOR: Algorithm found in http://code.activestate.com/recipes/576489/
% LAST CHANGES:
%  - 13/10/2008: if the scenario input is a vector of scenarios, only the 1st one is diplayed & char display modif
%  - 25/09/2008 : creation


if nargin < 2
    level = 0;
end

fn = fieldnames(scenario(1));
for n = 1:length(fn)
    tabs = '';
    for m = 1:level
        if n == length(fn) % dernier element 
            % use of unicode char U+2514 and U+2500
            tabs = [repmat('  ', 1, 2*level), '   └── '];
        else % element quelconque
            % use of unicode char U+251C and U+2500
            tabs = [repmat('  ', 1, 2*level), '   ├── '];
        end
    end
    disp([tabs fn{n}])
    fn2 = getfield(scenario(1), fn{n});
    if isstruct(fn2)
        aloha_scenario_disp(fn2, level+1);
    end
end

