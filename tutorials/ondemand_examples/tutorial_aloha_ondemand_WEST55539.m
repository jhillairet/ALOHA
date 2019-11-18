% ALOHA tutorial : ALOHA simulation of Tore Supra shot #45525
% J.Hillairet
% October 2011
%
%
% In this shot, two phase shifts between modules are used in order
% to modify the launched n// main peak from 1.7 to 1.8. We want to
% calculate the launched spectra of during this two steps of phase. This
% results may be used in LUKE for example.
%

% Tore Supra pulse
TSpulse = 45524;

% Tore Supra antenna port
% C4 - PAM launcher
TSport = 'Q6B'; 

% -180� phasing => n//=1.7
t_start1 = 12.3; % s 
t_stop1  = 13.0; % s

% -150� phasing => n//=1.8
t_start2 = 25; % s
t_stop2  = 25.5; % s

% First of all, before starting the ALOHA processing, we must be sure
% that the measured data are OK. For example, sometime the phase
% measurement can be crappy and thus leads to strange & incorrect results...
aloha_ondemand_plotPower(TSpulse, TSport, t_start1, t_stop1);

disp('press a key to continue')
pause

% OK, the data are OK at least for the upper modules on this time interval.
% So let's create the ALOHA scenario corresponding to this interval
% The following function will automatically select the good antenna
% description (C4 here) and will get from the Tore Supra database the input
% excitation of the modules (power and phase input)
%
% We have to provide to this function the edge electron density as well as
% the linear profile gradient (scrape-off length). A hint on the edge
% density can come from the Langmuir probes. The scrape-off gradient is
% more delicate to determine with confidence. The following are general
% results (2mm then 2cm) which of course can be tuned.
sc_TS45525_12s = aloha_ondemand_scenario(TSpulse, t_start1, t_stop1, 3e17, [2e-3,2e-2], TSport);
sc_TS45525_25s = aloha_ondemand_scenario(TSpulse, t_start2, t_stop2, 3e17, [2e-3,2e-2], TSport);

% here the scenarios have been generated, but not yet processed by ALOHA.
% Let's do it !
% (here you can go get a good tea, it is good for the health)
sc_TS45525_12s = aloha_scenario(sc_TS45525_12s);
sc_TS45525_25s = aloha_scenario(sc_TS45525_25s);

% save the processed results 
aloha_scenario_save(sc_TS45525_12s, 'sc_TS45525_12s.mat');
aloha_scenario_save(sc_TS45525_25s, 'sc_TS45525_25s.mat');


% The matlab command line display the following ALOHA results :
% between 12.3-13s : 3.01% average power reflected ; peak n//=1.71
% between 26.0-27s : 2.75% average power reflected ; peak n//=1.81

% now we can post-process our data. First, let's have a look to the coupled
% spectra
aloha_plot_spectra([sc_TS45525_12s, sc_TS45525_25s]);
    legend('TS45525 - 12s','TS45525 - 25s');
    set(gca, 'XLim', [-5, +5]);
aloha_plot_export(gcf, 'TS45525_spectra.pdf');    % could be .png, .fig, etc.
    
% Eventually we can export the ALOHA result to LUKE. 
% This is easy the following functions 
% aloha_scenario_conversion4luke(sc_TS45525_12s, 'LUKE_TS45525_12s.mat');
% aloha_scenario_conversion4luke(sc_TS45525_25s, 'LUKE_TS45525_25s.mat');


