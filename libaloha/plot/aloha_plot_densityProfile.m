function h=aloha_plot_densityProfile(varargin)
%  Plot the density profile modelled in front of the antenna
%  in function of the normalized radial position
% 
%  INPUT ARGUMENTS :
%   
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

    % generate the profile line
    switch(nargin)
        case 2 % ne0, lambda_n0
            ne0 = varargin{1};
            lambda_n0 = varargin{2};

            X = [X_PLASMA, X_ANTENNA];
            Y = [ne0*(1+abs(X(end)-X(1))/lambda_n0), ne0];

        case 4 % ne0, lambda_n0, d_couche, lambda_n1
            ne0 = varargin{1};
            lambda_n0 = varargin{2};
            d_couche = varargin{3};
            lambda_n1 = varargin{4};

            X = [X_PLASMA, X_ANTENNA-d_couche, X_ANTENNA];
            Y = [ne0*(1+d_couche/lambda_n0) * (1+abs(X_PLASMA-(X_ANTENNA-d_couche))/lambda_n1), ...
                 ne0*(1+d_couche/lambda_n0), ...
                 ne0];

        case 5 % d_vide, ne0, lambda_n0, d_couche, lambda_n1

        otherwise   
            error('invalid number of arguments. see help.');
    end

    % plot
    h=figure;
    line(X,Y);
    grid on;
    xlabel('\rho/R');
    ylabel('Electr. dens. n_e [m^{-3}]');
    title('Electronic density in front of the antenna');




