function aloha_plot_antenna(architecture)
% plot a graphical represention of an antenna architecture (as view from inside the torus).
% 
% aloha_plot_antenna(architecture)
% 
% INPUT
%  - architecture <string> : 
%  NB : This name must correspond to an existing filename into the antenna architectures directory.
% 
% OUPUT: none
%  
% AUTHOR : JH
% LAST UPDATE:
%  - 04/03/2009 : fill in gray passive wg
%  - 28/11/2008 : reverse the z axis in order to represent the antenna as view from the plasma
%  - 10/2008 : creation 

if exist(architecture)  
    eval(architecture);
else
    error(['Antenna architecture ', architecture, 'doesn''t exist !']);
end

[b, h, z, y, nb_g_total_ligne, nbre_guides,act_module_tor] = aloha_utils_getAntennaCoordinates(architecture); 

h=aloha_plot_figure(figure, [architecture, ' geometry']);

% view from the plasma, the z axes goes from left to right
% and y axes goes from top to bottom
% in this direct frame (x,y,z), x goes toward the plasma... 
set(gca, 'YDir', 'reverse');

% plot all waveguides for all poloidal line 
hold on;
for idx_pol = 1:nb_g_pol%nb_g_module_pol
    for idx_tor = 1:length(z)
        rect_pos = [z(idx_tor) y(idx_pol) b(idx_tor) a];
        
        % if the wg is a passive wg, 
        % then fill the rectangle in gray
        if (exist('b_g_pass') & (b(idx_tor) == b_g_pass))
          rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
        elseif (exist('b_g_pass_ext') & (b(idx_tor) == b_g_pass_ext))
          rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
        % for all other cases, just plot the contour
        % of the wg
        else 
            rectangle('Position', rect_pos);
        end
    end
end
hold off;

% plot legends & options 
axis equal;
xlabel('z [m]');
ylabel('y [m]');
title(['Antenna architecture :',  aloha_utils_str4fig(architecture), sprintf('\n'), '(as view from the plasma)']);
