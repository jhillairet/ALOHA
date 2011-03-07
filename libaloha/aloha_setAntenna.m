function scenario = aloha_setAntenna(scenario, antenna_name)
% ALOHA
% 
% Defines the ITM antenna_lh CPO for a scenario.
%
% INPUT
%  - scenario <structure>: ALOHA scenario
%  - antenna_name [string] : antenna name (ITM antenna description)
%
% OUTPUT
%  - scenario <structure>: ALOHA scenario with the antenna defined in it.
%
% AUTHOR: JH

% TODO : should we add array of scenarios support ?

% Only load antenna if the antenna description 
% comply to the ITM-way of definition, ie if is is a function
% (in contrast to a simple script)
try   
   scenario.antenna_lh = eval(antenna_name);
catch
   disp(aloha_message('Old fashioned antenna description... passing.'));
   scenario.antenna_lh.setup = [];
end
