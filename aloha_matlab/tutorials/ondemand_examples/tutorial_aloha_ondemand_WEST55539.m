% ALOHA tutorial : ALOHA simulation of WEST shot #55539
% J.Hillairet
% October 2019
%
%
% In this shot, two phase shifts between modules are used in order
% to modify the launched n// main peak from 1.7 to 1.8. We want to
% calculate the launched spectra of during this two steps of phase. This
% results may be used in LUKE for example.
%

%%
% WEST pulse
TSpulse = 55539;

% WEST antenna port
% LH1 (FAM Antenna) : 'Q6A'
% LH2 (PAM Antenna) : 'Q6B'
TSport = 'Q6A'; 

% start and stop times for measurement averaging
t_start1 = 8; % s 
t_stop1  = 9; % s

%%
% First of all, before starting the ALOHA processing, we must be sure
% that the measured data are OK. For example, sometime the phase
% measurement can be crappy and thus leads to strange & incorrect results...
aloha_ondemand_plotPower(TSpulse, TSport, t_start1, t_stop1);

disp('press a key to continue')
pause

%%
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

% edge density [m^-3]
ne0 = 3e17;
% density scrape-off length [m] for the first and second plasma layers
lambda_n = [2e-3, 2e-2];

sc_WEST55538 = aloha_ondemand_scenario(TSpulse, t_start1, t_stop1, ne0, lambda_n, TSport);

%%
% here the scenarios have been generated, but not yet processed by ALOHA.
% Let's do it !
% (here you can go get a good tea, it is good for the health)
sc_WEST55538 = aloha_scenario(sc_WEST55538);

%%
% save the processed results 
aloha_scenario_save(sc_WEST55538, 'sc_WEST55538.mat');


%%
% now we can post-process our data. First, let's have a look to the coupled
% spectra
aloha_plot_spectra([sc_WEST55538]);
    legend('WEST 55539');
    set(gca, 'XLim', [-5, +5]);
    
%%
aloha_plot_export(gcf, 'WEST55539_spectra.pdf');    % could be .png, .fig, etc.
    
%% Eventually we can export the ALOHA result to LUKE. 
% This is easy the following functions 
aloha_scenario_conversion4luke(sc_WEST55538, 'LUKE_WEST55538.mat');


