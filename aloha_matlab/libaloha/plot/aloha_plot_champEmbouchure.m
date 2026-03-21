function aloha_plot_champEmbouchure(scenario)
%   aloha_plot_champEmbouchure(scenario)
%   
%   Plot the electric field at the mouth of the antenna.
%   
%   
switch scenario.options.version_code
    case '1D'
        aloha_plot_champEmbouchure1D(scenario)

    case '2D'
        aloha_plot_champEmbouchure2D(scenario)

    otherwise
        error('Bad version_code parameter')
end