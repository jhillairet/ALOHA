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
        sig_RC    = 'GCREFC1';

    case 'Q6B' % C3 or C4
        sig_RC    = 'GCREFC2';
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
[refpow, t_refpow] = tsbase(pulsenb, sig_RC);

%% define some constants
upper_modules = [1:2:15];
lower_modules = [2:2:16];
%% cut time
idx_time = (t_refpow >= t_start) & (t_refpow <= t_stop);

t_refpow = t_refpow(idx_time(:,1),:);
refpow = refpow(idx_time(:,1),:);

T_MIN = min(t_refpow (:,1)); % s
T_MAX = max(t_refpow (:,1)); % s

FIG_WIDTH  = 600;
FIG_HEIGHT = 800;

%% plotting the reflexion coefficient
aloha_plot_figure('RC per modules - lower (blue) - upper (red)');
    set(gcf, 'Position', [20+FIG_WIDTH 100 FIG_WIDTH FIG_HEIGHT])
    for idx=1:length(upper_modules)
        subplot(4,2,idx)
            plot(t_refpow(:,lower_modules(idx)), refpow(:,lower_modules(idx)), 'b', ...
                 t_refpow(:,upper_modules(idx)), refpow(:,upper_modules(idx)), 'r');
            set(gca, 'YLim', [0 20]);
            set(gca, 'XLim', [T_MIN, T_MAX]);
            title(['Module: #',num2str(lower_modules(idx)), '(blue)', ...
                          ' #',num2str(upper_modules(idx)), '(red)']);
            ylabel('RC [%]');
            grid on;
            if idx==length(upper_modules)-1 | idx==length(upper_modules) % last two
                xlabel('t [s]');
            end
    end

%% average RC
avg_RC = mean(refpow,1)';
aloha_plot_figure('Average RC per module');
    subplot(211)
        plot([1:8], avg_RC(lower_modules));
        grid on;
        xlabel('Lower module #');
        ylabel('RC [%]');
        title('Lower modules')

    subplot(212)
        plot([1:8], avg_RC(upper_modules));
        grid on
        xlabel('Upper module #');
        ylabel('RC [%]');
        title('Upper modules')