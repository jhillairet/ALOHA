function aloha_ondemand_plotPower(pulsenb, port, varargin)
% aloha_ondemand_plotPower(pulsenb, port)
% 
% ALOHA on-demand : trace power and phase of a TS pulse
% 
% INPUT
%  - pulsenb: TS pulse number
%  - port (string): TS port 'Q6A' or 'Q6B'
%  - t_start: [optionnal] 
%  - t_stop : [optionnal]
%  
% OUTPUT
%  none
% 
% Author: JH
%  

%% Select the appropriate signal to plot depending on the TS port
switch upper(port)
    case 'Q6A' % C2 or C3
        sig_power = 'GPINJC1';
        sig_phase = 'GPHIC1';

    case 'Q6B' % C3 or C4
        sig_power = 'GPINJC2';
        sig_phase = 'GPHIC2';

    otherwise
        error('bad port definition. See help');
end

%% test optionnal arguments if provided
if nargin > 2
    if nargin == 3
        error('Please provide t_stop too. See help.')
    elseif nargin == 4
        t_start = varargin{1};
        t_stop = varargin{2};
    end
end

%% retrieving the signal
[power, t_power] = tsbase(pulsenb, sig_power); % in kW
[phase, t_phase] = tsbase(pulsenb, sig_phase);

%% define some constants
upper_modules = [1:2:15];
lower_modules = [2:2:16];

MAX_POWER = 300;% kW
T_MIN = min(t_power(:,1)) + 5; % s
T_MAX = max(t_power(:,1)) - 5; % s
FIG_WIDTH  = 600;
FIG_HEIGHT = 800;

%% plotting the power signals
aloha_plot_figure('Power per modules - lower (blue) - upper (red)');
    set(gcf, 'Position', [20 100 FIG_WIDTH FIG_HEIGHT])
    for idx=1:length(upper_modules)
        subplot(4,2,idx)
            plot(t_power(:,lower_modules(idx)), power(:,lower_modules(idx)), 'b', ...
                 t_power(:,upper_modules(idx)), power(:,upper_modules(idx)), 'r');
            set(gca, 'YLim', [0 MAX_POWER]);
            set(gca, 'XLim', [T_MIN, T_MAX]);
            title(['Module: #',num2str(lower_modules(idx)), '(blue)', ...
                          ' #',num2str(upper_modules(idx)), '(red)']);
            ylabel('Power [kW]');
            grid on;
            if idx==length(upper_modules)-1 | idx==length(upper_modules) % last two
                xlabel('t [s]');
            end
            if exist('t_start')
                rectangle('Position', [t_start 0 t_stop-t_start MAX_POWER]);
            end
    end
    htb = annotation(gcf, 'textbox',[0.4 0.95 0.4 0.03])
    set(htb, 'String', ['#', num2str(pulsenb),' ; ', port,' ; Blue: bottom ; Red: top']);

%% plotting the phase signals
aloha_plot_figure('Phase per modules - lower (blue) - upper (red)');
    set(gcf, 'Position', [20+FIG_WIDTH 100 FIG_WIDTH FIG_HEIGHT])
    for idx=1:length(upper_modules)
        subplot(4,2,idx)
            plot(t_phase(:,lower_modules(idx)), phase(:,lower_modules(idx)), 'b', ...
                 t_phase(:,upper_modules(idx)), phase(:,upper_modules(idx)), 'r');
            set(gca, 'YLim', [0 360]);
            set(gca, 'XLim', [T_MIN, T_MAX]);
            title(['Module: #',num2str(lower_modules(idx)), '(blue)', ...
                          ' #',num2str(upper_modules(idx)), '(red)']);
            ylabel('Phase [deg]');
            grid on;
            if idx==length(upper_modules)-1 | idx==length(upper_modules) % last two
                xlabel('t [s]');
            end
            if exist('t_start')
                rectangle('Position', [t_start 0 t_stop-t_start 360]);
            end
    end