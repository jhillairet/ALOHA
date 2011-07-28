function aloha_plot_reflectionCoeff(scenario, varargin)
% Plot the reflection coefficients from a scenario structure. 
% 
% If more than 1 scenario are present into the structure, all
% the RC are plotted on the same figure. 
% 
% By default, the reflection coefficient for each module is plotted.
% 
% If a vector is set as optionnal input, then the module average is plot vs this vector. 
% Thus, the optionnal vector length should match the length of the scenario.
%  
% INPUT 
%  - scenario [struct(n)] : structure scenario (supporte plusieurs scenario a la fois)
%  - [optionnal] x : vector (1xlength(scenario))
%  - [optionnal] label_x : label of the vector. Must be set if an optionnal vector is set.

% EXAMPLE
%  aloha_plot_reflectionCoeff(scenario [,x, label_x])
%  
% AUHTOR:JH
% LAST CHANGES:
%  - 02/06/2010: add possibility to specify a variation vector (like x=ne0)
%  - 06/11/2008: minor changes & use of aloha functions
%  - 10/09/2008: creation
%  

% on recupere tous les coefficients de reflexion calcules
RC = aloha_scenario_get(scenario, 'CoeffRefPuiss');

% verification de la validite du scenario
if isempty(RC)
    error('Results don''t exist in the scenario structure(s)');
end

if nargin == 1

% plot
aloha_plot_figure(figure, 'ALOHA : module reflection coefficients');
    plot(1:size(RC,2), RC);
    grid on;
    xlabel('module #');
    ylabel('Reflection coefficient (%)');
    % x label is only integer !
    set(gca, 'XTick', 1:size(RC,2));

elseif nargin >= 2
    x = varargin{1};

    if length(x) ~= length(scenario)
        error('The length of the optionnal vector does not match the scenario length');
    end

    aloha_plot_figure(figure, 'ALOHA : average module reflection coefficients');
        plot(x, mean(RC,2));
        ylabel('mean RC [%]');
        grid on;

        if nargin == 3
            label_x = varargin{2};
            xlabel(label_x);
        end
end