% represente le champ electrique dans une ligne de guides avec les guides en superposition derriere.
function aloha_plot_Efield_withAntenna(scenario)
% plot the absolute value of the Ez electric field in a row of waveguide
% with the waveguide surimposed on the figure
% 
% [Efield_mean=]aloha_plot_Efield_withAntenna(scenario)
% 
% INPUT
%  scenario : ALOHA scenario. (Only the first one is used.)
%  

% take the first scenario if there is many
if length(scenario) > 1
    scenario = scenario(1);
end

% if the averaged field does not exist
% we compute it
if not(isfield(scenario.results, 'Ez_average'))
    scenario = aloha_compute_averageEz_waveguides(scenario);
end

% get the dimensions and coordinates of the waveguides of the antenna
[b,h,z,y,nb_g_total_ligne,nbre_guides]=aloha_utils_getAntennaCoordinates(scenario.antenna.architecture);


MAX_Efield = 1.2*max(abs(scenario.results.Ez_average));

PassWG_idx = scenario.results.PassWG_idx;

aloha_plot_figure(figure)
    % draw rectangles for waveguides
    for idx_wg=1:nb_g_total_ligne
        if ismember(idx_wg,PassWG_idx)
            % passive waveguide
            rectangle('Position', [z(idx_wg) 0 b(idx_wg) MAX_Efield], 'FaceColor', [.8 .8 .8]);
        else
            % active waveguide
            rectangle('Position', [z(idx_wg) 0 b(idx_wg) MAX_Efield]);
        end
    end
    % surimpose the Efield
    hold on
    plot(scenario.results.abs_z(1,:), abs(scenario.results.Efield(1,:)))
    hold off

% surimpose the average field in a waveguide
for idx_wg=1:nb_g_total_ligne
    line([z(idx_wg) z(idx_wg)+b(idx_wg)], scenario.results.Ez_average(idx_wg)*[1 1], 'Color', 'r');
end

