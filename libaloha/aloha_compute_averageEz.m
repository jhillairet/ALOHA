function [varargout] = aloha_compute_averageEz(scenario)
% Ez_mean = aloha_compute_averageEz(scenario)
% or
% [Ez_mean, Ez_max] = aloha_compute_averageEz(scenario)
% or
% [Ez_mean, Ez_max, Ez_min] = aloha_compute_averageEz(scenario)
%  
%  
% Compute the average electric field at the mouth of the grill.
% This function didn't take into account the 0-field zones which 
% correspond to the septum of the mouth.
%  
% INPUT
%  - scenario : ALOHA scenario <struct Nx1>
% 
% OUTPUT
%  - Ez_mean : average electric field (nb_WGline x N)
%  - Ez_max : maximum electric field 
%  - Ez_min : minimum electric field
%  
% AUTHOR: J.Hillairet
% LAST UPDATE : 
%  - 23/02/2010 : creation

for idx_sc=1:length(scenario)
    % check is the E field have been calculated
    if not(isfield(scenario(idx_sc).results, 'Efield'))
        error('The electric field had not been calculated !');
    end
    
    % retrieve E field from ALOHA scenario
    [Ez]=aloha_scenario_get(scenario(idx_sc), 'Efield');
    
    % supress the zero E field points
    Ez = squeeze(abs(Ez)); 
    idx = find(Ez == 0); 
    Ez(idx)=[];
    
    % output
    Ez_mean(idx_sc) = mean(Ez);
    Ez_max(idx_sc)  = max(Ez); 
    Ez_min(idx_sc)  = min(Ez);

end

% optionnal output arguments
if nargout <= 1
    varargout{1} = Ez_mean;
elseif nargout == 2
    varargout{1} = Ez_mean;
    varargout{2} = Ez_max;
elseif nargout == 3
    varargout{1} = Ez_mean;
    varargout{2} = Ez_max;
    varargout{3} = Ez_min;
end