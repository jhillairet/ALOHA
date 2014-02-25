function phi_unique = aloha_plot_angle(scenario)
% Get the phase 

%% get antenna coordinates
[b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(scenario.antenna.architecture);

%% get E field and angle
z_Efield = scenario.results.abs_z(1,:);
phi = ((angle(scenario.results.Efield(1,:))));

%% get the phase angle for each waveguide
for idx=1:length(z)
    % get the phase of the field
    % when abscisse corresponds to a waveguide position
    phi_unique(idx) = phi(find(z_Efield == z(idx)));
end

wg_idx = 1:length(phi_unique);

%% plot results
%% absolute value of the Efield phase at the mouth
aloha_plot_figure(figure)
    bar(wg_idx,180/pi*(phi_unique));
    grid on;
    xlabel('Waveguide #');
    ylabel('\phi_{wg} [deg]')

%% same thing unwrapped
aloha_plot_figure(figure)
    bar(wg_idx,180/pi*unwrap(phi_unique));
    grid on;
    xlabel('Waveguide #');
    ylabel('\phi_{wg} (unwrapped) [deg]')

% Calcul de regression lineaire pour determiner la pente 
% (ie l'increment de phase d'un guide a l'autre)
P = polyfit(wg_idx, unwrap(phi_unique), 1);
% surimpression du fit
line([wg_idx(1) wg_idx(end)], 180/pi*polyval(P, [wg_idx(1) wg_idx(end)]) , 'Color', 'r');
title(['Linear fitting : Y=', num2str(180/pi*P(1)),'[deg] X' ])


%% Calculation of the phase shift between waveguides
delta_phi = 180/pi*(phi_unique(2:end)-phi_unique(1:end-1));
% add 360° when phase shift is negative in order to have only positive values in [0 360]
delta_phi(delta_phi<0) = delta_phi(delta_phi<0)+360;

% calculating average phase shift
delta_phi_avg = mean(delta_phi);

aloha_plot_figure(figure)
    bar(wg_idx(1:end-1),delta_phi);
    grid on;
    xlabel('Waveguide #');
    ylabel('\Delta\phi_{wg} (\phi_{p+1}-\phi_p) [deg]')
    line([wg_idx(1) wg_idx(end)], [1 1]*delta_phi_avg, 'Color', 'r');
    title(['Average delta phi=',num2str(delta_phi_avg),' deg']);