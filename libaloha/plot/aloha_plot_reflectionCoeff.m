function aloha_plot_reflectionCoeff(scenario)
% Plot the reflection coefficients from a scenario structure. 
% 
% If more than 1 scenario are present into the structure, all
% the RC are plotted on the same figure. 
% 
%  
% INPUT 
%  - scenario [struct(n)] : structure scenario (supporte plusieurs scenario a la fois)

% EXAMPLE
%  aloha_plot_reflectionCoeff(scenario)
%  
% AUHTOR:JH
% LAST CHANGES:
%  - 06/11/2008: minor changes & use of aloha functions
%  - 10/09/2008: creation
%  

% on recupere tous les coefficients de reflexion calcules
RC = aloha_scenario_get(scenario, 'CoeffRefPuiss');

% verification de la validite du scenario
if isempty(RC)
    error('Results don''t exist in the scenario structure(s)');
end

% plot
aloha_plot_figure(figure, 'ALOHA : module reflection coefficients');
    plot(1:size(RC,2), RC);
    grid on;
    xlabel('module #');
    ylabel('Reflection coefficient (%)');
