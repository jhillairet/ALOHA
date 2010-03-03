function aloha_plot_spectra(scenarios, varargin)  
% Plot 1D (ie nz) spectra for different scenarios
%  
% EXAMPLE
%  aloha_plot_spectrum(scenario)
%  
% INPUT
%  - scenario [struct(1)] : structure scenario
%  [optionnal] - 'legend' : string vector which contains the legend of the plot
%  
% AUTHOR: JH
% LAST CHANGES
%  - 02/10/2009 add a legend optionnal input and also works with aloha-2D
%  - 03/07/2009 creation
%    

h=aloha_plot_figure(figure, 'ALOHA : nz spectra');
colors = [get(gca,'ColorOrder');get(gca,'ColorOrder');get(gca,'ColorOrder')];
hold on;

% default optionnal input argument value
bool_normalization = false;

% parsing optionnal input argument
if nargin > 1
    for idx=1:2:(nargin-1)
        arg_name = varargin{idx};
        arg_value= varargin{idx+1};
        switch lower(arg_name)
            case 'normalization'
                bool_normalization = arg_value;
        end
    end
end


for id_scen=1:length(scenarios) % for many scenarios
    scenario=scenarios(id_scen);

    current_dP_nz = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
    current_nz = squeeze(aloha_scenario_get(scenario, 'nz'));

    if bool_normalization
        norma = max(real(current_dP_nz));
    else
        norma = 1;
    end

    plot(current_nz, real(current_dP_nz)./norma, 'Color', colors(id_scen,:));
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



