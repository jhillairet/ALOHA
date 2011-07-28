function bool = aloha_isAntennaITM(scenario)
% ALOHA
%
% Says if the passed scenario contains an ITM antenna description
%
% INPUT
% - scenario : ALOHA scenario
%
% OUTPUT
% - bool [boolean] : true for an ITM description, fakse otherwise
%
% AUTHOR: JH
try
    if isfield(scenario.antenna_lh.setup, 'modules')
        bool = true;
    else
        bool = false;
    end
catch % there was an error in the previous line : this may indicate that this is an old antenna description
    bool = false;
end
