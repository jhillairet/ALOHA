function aloha_plot_spectra(scenarios, varargin)  
% Plot 1D (ie nz) spectra for different scenarios
%  
% EXAMPLE
%  aloha_plot_spectrum(scenario)
%  
% INPUT
%  - scenario [struct(1)] : structure scenario
%  [optionnal] - legend : string vector which contains the legend of the plot
%  
% AUTHOR: JH
% LAST CHANGES
%  - 02/10/2009 add a legend optionnal input and also works with aloha-2D
%  - 03/07/2009 creation
%    

h=aloha_plot_figure(figure, 'ALOHA : nz spectra');
colors = get(gca,'ColorOrder');
hold on;

for id_scen=1:length(scenarios) % for many scenarios
    scenario=scenarios(id_scen);

    current_dP_nz = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
    current_nz = squeeze(aloha_scenario_get(scenario, 'nz'));

    plot(current_nz, real(current_dP_nz), 'Color', colors(id_scen,:));
    xlabel('n_{z}');
    ylabel('Power density [W]');
    title('Power density spectrum');
    grid on;

end % id_scen

% add the legend if the optionnal input has been given
if nargin == 2
    legend(varargin{1});
end

hold off;



