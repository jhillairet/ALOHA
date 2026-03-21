function h=aloha_plot_densityProfile(scenarios)
%  Plot the density profile modelled in front of the antenna
%  in function of the normalized radial position
% 
%  INPUT ARGUMENTS :
%   - scenario
%  
%  OUPUT: 
%   - h : figure handler
%  
%  AUTHOR : JH
%  LAST UPDATES : 
%  - 31/07/2008 : creation
% 

    % normalized positions
    X_ANTENNA= 3;
    X_PLASMA = 1;
    % large radius of the Tokamak
    R = 3;

for idx=1:length(scenarios)
    scenario=scenarios(idx);

    % retrieve scenario plasma configuration
    ne0      = aloha_scenario_get(scenario, 'ne0');
    lambda_n = aloha_scenario_get(scenario, 'lambda_n');
    d_couche = aloha_scenario_get(scenario, 'd_couche');
    d_vide   = aloha_scenario_get(scenario, 'd_vide');

    % generate the profile line
    switch(aloha_scenario_get(scenario, 'version'))
        case 3 % ne0, lambda_n0
            lambda_n0 = lambda_n(1);

            X(idx,:) = [X_PLASMA, X_ANTENNA];
            Y(idx,:) = [ne0*(1+abs(X(end)-X(1))/lambda_n0), ne0];

        case 6 % d_vide, ne0, lambda_n0, d_couche, lambda_n1
            lambda_n0 = lambda_n(1);
            lambda_n1 = lambda_n(2);
            % gradients
            dne0 = ne0./lambda_n0;
            dne1 = (1+d_couche./lambda_n0).*ne0./lambda_n1;

            ne1 = ne0 + dne0*d_couche;
            ne2 = ne1 + dne1*abs(X_PLASMA-(X_ANTENNA-d_couche-d_vide));
            % create the 
            X(idx,:) = [X_PLASMA, X_ANTENNA-d_couche-d_vide, X_ANTENNA-d_vide,X_ANTENNA-d_vide, X_ANTENNA];
            Y(idx,:) = [ne2, ne1, ne0,0, 0];
            

        otherwise   
            error('invalid number of arguments. see help.');
    end
end

  
    % plot
    aloha_plot_figure(figure)
    line((X_ANTENNA-X'),Y', 'Marker', 'O');
    line([0,0], [0, 20e18], 'LineStyle','--', 'color', 'k');
    set(gca, 'XLim', [-1e-3, 7e-3], 'XDir', 'reverse');
    set(gca, 'YLim', [0 20e18]);

    
    grid on;
    xlabel('x [mm]');
    ylabel('Electr. dens. n_e [m^{-3}]');
    title('Electronic density in front of the antenna');

  


