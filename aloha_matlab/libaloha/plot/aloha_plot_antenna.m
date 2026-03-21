function aloha_plot_antenna(varargin)
% plot a graphical represention of an antenna architecture (as view from inside the torus).
% 
% aloha_plot_antenna(architecture)
%   or
% aloha_plot_antenna(scenario)
% 
% INPUT
%  - architecture <string> : architecture_name
%  NB : This name must correspond to an existing filename into the antenna architectures directory.
%   OR 
%  - scenario : ALOHA scenario
%
% OUPUT: none
%  
% AUTHOR : JH
% LAST UPDATE:
%  - 03/2011 : Add a first ITM support 
%  - 04/03/2009 : fill in gray passive wg
%  - 28/11/2008 : reverse the z axis in order to represent the antenna as view from the plasma
%  - 10/2008 : creation 

%% Test input arguments
if isstr(varargin{1})
    architecture = varargin{1};
    if exist(architecture) 
          try % Old fashionned antenna
            eval(architecture);
            [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(architecture);
            scenario = struct();
        catch % ITM antenna
              disp(aloha_message('assuming ITM antenna description')); 
              scenario=aloha_setAntenna(struct(), architecture);
              aloha_utils_ITM2oldAntenna; % convert ITM data to old fashioned ALOHA global parameters
            [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
        end
    else
         error(['Antenna architecture ', architecture, 'doesn''t exist !']);
    end
elseif isstruct(varargin{1});
    scenario = varargin{1};
    [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
    architecture = squeeze(aloha_scenario_get(scenario, 'archName'))';
    aloha_utils_ITM2oldAntenna; % convert ITM data to old fashioned ALOHA global parameters
else
    error('Bad input argument. See help of this function or check the input argument');
end



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
        
        % new way
        if aloha_isAntennaITM(scenario)
            % we create an array of passive/active mask to determine 
            % if a waveguide is either active of passive
            ar_modules = ones(1, scenario.antenna_lh.setup.modules.nma_phi); % [1 1 ... 1] (1xnmodules)
            ar_pa_mask = scenario.antenna_lh.setup.modules.waveguides.mask;
            % adding passive waveguide between modules in the P/A mask
            ar_pa_mask = [ar_pa_mask, zeros(1, scenario.antenna_lh.setup.modules.waveguides.npwbm_phi)];
            
            pa_mask = kron(ar_modules, ar_pa_mask);
            % removing the last one, which came from the passive waveguide
            % between modules ) if there's any)
            if scenario.antenna_lh.setup.modules.waveguides.npwbm_phi > 0
                pa_mask(end)=[];
            end
            % adding passive waveguides at the edges in the P/A mask
            pa_mask = [zeros(1, scenario.antenna_lh.setup.modules.waveguides.npwe_phi), pa_mask, zeros(1, scenario.antenna_lh.setup.modules.waveguides.npwe_phi)];
            
            if pa_mask(idx_tor) == 0
                rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
            elseif pa_mask(idx_tor) == 1
                rectangle('Position', rect_pos);
            end
%              % if the current wg is a passive wg
%              if b(idx_tor) == aloha_scenario_get(scenario, 'bewp');
%               rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
%              elseif b(idx_tor) == aloha_scenario_get(scenario, 'biwp');
%               rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
%              % for all other cases (ie active wg), just plot the contour
%              % of the wg
%              else 
%                  rectangle('Position', rect_pos);
%              end
        else % old way
            % if the wg is a passive wg, 
            % then fill the rectangle in gray
            if (exist('b_g_pass') & (b(idx_tor) == b_g_pass))
              rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
            elseif (exist('b_g_pass_ext') & (b(idx_tor) == b_g_pass_ext))
              rectangle('Position', rect_pos, 'FaceColor', [.8,.8,.8]);
            % for all other cases (ie active wg), just plot the contour
            % of the wg
            else 
                rectangle('Position', rect_pos);
            end
        end
    end
end
hold off;

% plot legends & options 
axis equal;
xlabel('z [m]');
ylabel('y [m]');
title(['Antenna architecture :',  aloha_utils_str4fig(architecture), sprintf('\n'), '(as view from the plasma)']);
