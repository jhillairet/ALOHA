function aloha_plot_champEmbouchure1D(scenario)
% ALOHA
% 
% aloha_plot_champEmbouchure2D(scenario)
%
% Plot the electric field (parallel component) at the mouth of the antenna
% 
% INPUT
%  - scenario (1x1) : ALOHA scenario (2D plasma)
% 
% OUTPUT : none
% 
% AUTHORS: DV,JH
% 
% 

%%%%%%% AMPLITUDE in 2D

h=aloha_plot_figure;
hold on;

for idx_wg = 1:scenario.results.nbre_guides
    
    Ey_g = squeeze(scenario.results.Ey(idx_wg,:,:));
    Ez_g = squeeze(scenario.results.Ez(idx_wg,:,:));
    z_g = scenario.results.z(idx_wg,:);
    y_g = scenario.results.y(idx_wg,:);
    [ZG,YG]=ndgrid(z_g,y_g);
    norm_E = sqrt(abs(Ey_g).^2+abs(Ez_g).^2);

    pcolor(ZG, YG, norm_E.')
      
end

hold off;

shading interp;
xlabel('z [m]');
ylabel('y [m]');
title('Norm Electric field amplitude (V/m)');
colorbar;


%%%%%%%%% PHASE in 2D

h=aloha_plot_figure;
hold on;

for idx_wg = 1:scenario.results.nbre_guides

    Ez_g = squeeze(scenario.results.Ez(idx_wg,:,:));
    z_g = scenario.results.z(idx_wg,:);
    y_g = scenario.results.y(idx_wg,:);
    [ZG,YG]=ndgrid(z_g,y_g);
    norm_E = sqrt(abs(Ey_g).^2+abs(Ez_g).^2);

    pcolor(ZG, YG, 180/pi*angle(Ez_g))
      
end

hold off;

shading interp;
xlabel('z [m]');
ylabel('y [m]');
title('Parallel Electric field (E_z) phase (deg)');
colorbar;
    

%%%%%% AMPLITUDE in 1D
aloha_plot_figure
hold on;
 
for idx_wg = 1:scenario.results.nbre_guides
    Ey_g = squeeze(scenario.results.Ey(idx_wg,:,:));
    Ez_g = squeeze(scenario.results.Ez(idx_wg,:,:));
    z_g = scenario.results.z(idx_wg,:);
    y_g = scenario.results.y(idx_wg,:);

    % look for the middle of the waveguide (in y direction)
    % i.e. where the E field should be maximum
    y_middle = (max(y_g)-min(y_g))/2;
    % take the closest vector index from geometric middle
    [dummy, idx_y_middle] = min(abs(y_g-y_middle)); 
    % take fields for this value
    Ey_g = Ey_g(idx_y_middle,:);
    Ez_g = Ez_g(idx_y_middle,:);
    norm_E = sqrt(abs(Ey_g).^2+abs(Ez_g).^2);

    % plot the field
    plot(z_g, abs(Ez_g));

end

hold off;
xlabel('z [m]');
ylabel('|E_z| [V/m]');
grid on;
title('Parallel electric field amplitude');
