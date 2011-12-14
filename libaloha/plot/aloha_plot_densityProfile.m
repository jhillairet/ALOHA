function h=aloha_plot_densityProfile(scenario)
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
    X_PLASMA = 2.5;
    % large radius of the Tokamak
    R = 3;

    

    % retrieve scenario plasma configuration
    ne0      = aloha_scenario_get(scenario,'ne0');
    lambda_n = aloha_scenario_get(scenario, 'lambda_n');
    d_couche = aloha_scenario_get(scenario, 'd_couche');
    d_vide   = aloha_scenario_get(scenario, 'd_vide');

    % generate the profile line
    switch(aloha_scenario_get(scenario, 'version'))
        case 3 % ne0, lambda_n0
            ne0 = aloha_scenario_get()
            lambda_n0 = lambda_n(1);

            X = [X_PLASMA, X_ANTENNA];
            Y = [ne0*(1+abs(X(end)-X(1))/lambda_n0), ne0];

        case 6 % d_vide, ne0, lambda_n0, d_couche, lambda_n1
            lambda_n0 = lambda_n(1);
            lambda_n1 = lambda_n(2);
            % gradients
            dne0 = ne0./lambda_n0;
            dne1 = (1+d_couche./lambda_n0).*ne0./lambda_n1;

            ne1 = ne0 + dne0*d_couche;
            ne2 = ne1 + dne1*abs(X_PLASMA-(X_ANTENNA-d_couche-d_vide));
            % create the 
            X = [X_PLASMA, X_ANTENNA-d_couche-d_vide, X_ANTENNA-d_vide,X_ANTENNA-d_vide, X_ANTENNA];
            Y = [ne2, ne1, ne0,0, 0];
            

        otherwise   
            error('invalid number of arguments. see help.');
    end

    % plot
    h=figure;
    line(X,Y, 'Marker', 'O');
    grid on;
    xlabel('\rho/R');
    ylabel('Electr. dens. n_e [m^{-3}]');
    title('Electronic density in front of the antenna');




