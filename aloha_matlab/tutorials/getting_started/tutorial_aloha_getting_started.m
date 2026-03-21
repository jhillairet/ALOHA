% ALOHA tutorial: getting started with ALOHA
% J.Hillairet
% October 2011
%
% In this tutorial, we will create a ALOHA scenario for the Tore Supra PAM
% (C4) antenna. 
% Once the ALOHA scenario processed, we will post-process the result with
% some of the ALOHA functions.

% First of all, let's creating an ALOHA scenario from scratch. This may be
% done with the following command. This function will create a file in the
% current directory named 'My_first_aloha_scenario.m'
aloha_scenario_constructor('My first aloha scenario');

% if you run this matlab function, it will generate a matlab structure which
% is an 'ALOHA scenario'. This matlab structure contains all the
% parameters and informations needed by ALOHA to make its calculations.
% This structure is also the input variable of the aloha_scenario function.

disp('Please have a look at the My_first_aloha_scenario.m file first');
pause

% This tutorial contains two example scenarios
% - tutorial_aloha_scenarioC4_defaultExcitation.m
% - tutorial_aloha_scenarioC4_manualExcitation.m
% The C4 antenna is described internally in the ALOHA code. One can use
% the default excitation for this antenna or to define its own power and phase excitation. 
%
% First, let's run ALOHA on the default excitation scenario. 
% The default C4 ALOHA excitation  correspond to a power density of
% 25MW/m^2 at the mouth. In ALOHA, only half of the Tore Supra are modeled. 
% Since there is 8 modules (8 upper or 8 lower), this means that the default 
% power input in a module is 2.6725MW/(2*8) Watts. This value comes from
% the active waveguide surface at the mouth and the 25MW/m^2 power density.
% 
sc_defaultExcitation = aloha_scenario(tutorial_aloha_scenarioC4_defaultExcitation);
% We can read in the command line that the average reflection coefficient is 5.07% 
% and the peak n// is 1.72 

% Now, let's run ALOHA on the manual excitation. In the manual excitation
% scenario, I've put a -150° phase shift between module, and 300W per
% module input. I'm expecting a peak n// of 1.8 with this phase shift
% value. 
sc_manualExcitation = aloha_scenario(tutorial_aloha_scenarioC4_manualExcitation);
% Now we can read a peak n// of 1.82, as expected. Note also the power
% conservation lines: Pin=2400 W (ie. 8x300W) and the 
% Transmitted powers calculated with the reflection coefficient and with the
% spectrum integration are very close (the ratio is 1.017). 
% This indicates a good confidence on this result.

% We save these results in order to be able to re-use it later
aloha_scenario_save(sc_defaultExcitation, 'scenario_C4_defaultExcitation.mat');
aloha_scenario_save(sc_manualExcitation, 'scenario_C4_manualExcitation.mat');

% since we have saved our data, there no reasons being afraid of clean the
% matlab memory and to re-load our data. This totally needless, but it's a
% tutorial after all...
clear all;

sc_defaultExcitation=aloha_scenario_load('scenario_C4_defaultExcitation.mat');
sc_manualExcitation =aloha_scenario_load('scenario_C4_manualExcitation.mat');

% let's plot the coupled spectrum of the antenna
aloha_plot_spectrum(sc_defaultExcitation);

% or the electric field at the mouth
aloha_plot_champEmbouchure(sc_defaultExcitation);

% or the reflection coefficient at the output of each modules.
aloha_plot_reflectionCoeff(sc_defaultExcitation);

% or what the antenna front face look like for ALOHA
aloha_plot_antenna(sc_defaultExcitation);

% All the results are stored in the scenario structure.
% You can look at these results by typing in the command line:
sc_defaultExcitation.results

% for example, the reflection coefficient (in %) for each module is stored in the
% field :
sc_defaultExcitation.results.RC

