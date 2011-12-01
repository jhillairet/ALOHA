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
catch ME
   disp(aloha_message(ME.message));
   disp(aloha_message('Antenna description file not found or old fashioned antenna description... '));
   scenario.antenna_lh.setup = [];
end
