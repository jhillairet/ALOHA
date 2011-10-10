% ALOHA tutorial
% J.Hillairet October 2011
%
% In this tutorial we create a simple "grill" antenna (Brambilla antenna)
% which is a juxtaposition of phased waveguides.
% 
% After making the ALOHA antenna description and the plasma scenario, we
% test the effect of the waveguide phasing on the spectrum. 
%
% Then we create a "batch" file, which purpose is to make easily a lot of
% calculation. In this example, we select a 90°deg phase shift at the mouth
% and we vary the electron edge density to evaluate the effect on the
% coupling.

% First of all, we create an antenna that we will modify to fit our needs
aloha_antenna_constructor('antenna_simple_grill_8waveguides');
% this function had created a matlab file named
% "antenna_simple_grill_8waveguides.m"
% This file has been modfified to :
% "tutorial_aloha_antenna_simple_grill_8waveguides.m"
% Have a look to see the differences between the two files!
%

disp('OK to continue ?')
pause

% We have also created a Scattering matrix file named
% "tutorial_aloha_Selem.m" 
% which contains a perfect (loss free) waveguide scattering matrix. This
% matrix is the scattering matrix of all the waveguides of the antenna we
% have defined.
%

% Let's have a look to this antenna front face to see if it fits what we
% want:
aloha_plot_antenna('tutorial_aloha_antenna_simple_grill_8waveguides')

disp('OK to continue ?')
pause

% Now we create the associated ALOHA scenario
aloha_scenario_constructor('scenario_simple_grill_8waveguides');
% this function had created a matlab file named 
% "scenario_simple_grill_8waveguides.m"
% This has been modified to 
% "aloha_tutorial_scenario_simple_grill_8waveguides.m"
% Have a look to see the difference between the two files!
% 
% In particular, look at the excitation parameters. 
% For this example, each module of the antenna is in fact a unique 
% waveguide. So the excitation of the module should match the waveguide
% excitation. In the example, a 90° phase shift is applied between modules
% (and thus between waveguide), creating a peak n//=1.9
%
disp('OK to continue ?')
pause

% Now let's run ALOHA on the scenario !
sc = aloha_scenario(tutorial_aloha_scenario_simple_grill_8waveguides);

% ...A few seconds later...

% We save the scenario 
aloha_scenario_save(sc, 'scenario_simple_grill.mat');

% An plot some results
aloha_plot_spectrum(sc);
aloha_plot_reflectionCoeff(sc);

disp('OK to continue ?');
pause

% Now we want to evaluate the effect of the phase shift between waveguides
% on the coupled spectrum. To see this effect, we duplicate the initial
% scenario and we change the excitation phase law to test different phase
% shift between module (=waveguide here)
%
% Moreover, we make an array of scenario, which will make easier the
% post-processing after.
initial_scenario = tutorial_aloha_scenario_simple_grill_8waveguides;
scenarios = repmat(initial_scenario, 1, 6);

% and we change the phase shift for each of the scenarios:
scenarios(1).antenna.a_phase = 30*(pi/180)*(0:7)';
scenarios(2).antenna.a_phase = 60*(pi/180)*(0:7)';
scenarios(3).antenna.a_phase = 90*(pi/180)*(0:7)';
scenarios(4).antenna.a_phase = 120*(pi/180)*(0:7)';
scenarios(5).antenna.a_phase = 160*(pi/180)*(0:7)';
scenarios(6).antenna.a_phase = 180*(pi/180)*(0:7)';

% and we run ALOHA on these scenarios
scenarios = aloha_scenario(scenarios);
% ALOHA auto-detects that this is an array of 6 scenarios and make the
% calculations for all of the scenarios one after the other. You can see
% that in the command line which indicates:
% "########## Scenario 1/6 ##########"
% 

% ...A few seconds later... (x 6)

% We save the scenario 
% Saving and loading array of scenarios is the same that for unique
% scenario.
aloha_scenario_save(scenarios, 'scenarios_simple_grill.mat');

% Now we can compare the spectra
aloha_plot_spectra(scenarios);
    legend('30deg','60deg', '90deg', '120deg', '160deg', '180deg');

% As well as the reflection coefficient
aloha_plot_reflectionCoeff(scenarios);
    legend('30deg','60deg', '90deg', '120deg', '160deg', '180deg');
% we could also plot the average RC vs phase to see a sort of "tendancy"
% This can be done by giving a the x array as second input argument:
aloha_plot_reflectionCoeff(scenarios, [30,60,90,120,160,180]);
    xlabel('Phase shift between waveguides [deg]');
    title('Average RC VS waveguide phase shift');

disp('OK to continue ?');
pause

% From the previous figure, it seems that the antenna is better matched 
% for this plasma for 90° deg phasing. 
% Now we could want to evaluate the effect of the electron edge density 
% parameter on the coupling efficiency. For this purpose, we could do as
% previously in replicating X times the initial scenario and change the
% density for each element of the array of scenarios.
% (by changing the "scenarios(X).plasma.ne0" parameter).
% 
% However, when one wants to make a parameter scan with many point, it
% would be a bit inconvenient. So it's time to introduce the ALOHA batch,
% we purpose is to make a parameter scan easier !
% 
% Lets' create the batch file
aloha_batch_constructor('batch_simple_grill');
% this function had created the file "batch_simple_grill.m" which had been
% modified to "tutorial_aloha_batch_simple_grill". Look at the two files to
% see the differences !

disp('OK to go ? ')
pause

% Readay ? Go ! 
% We run this batch !
% It's time to go get a tea !

tutorial_aloha_batch_simple_grill;

disp('Let''s crunch the data !')
pause

% The bathc automatically saved the results in the matlab file :
% "tutorial_aloha_scenario_simple_grill_8waveguides". 
% So we just have to load the results
batch_sc = aloha_scenario_load('tutorial_aloha_scenario_simple_grill_8waveguides.mat');

% In order to plot the average RC VS ne0, we first get the array of ne0 used in the scenarios.
% When one want to get a data which is inside a scenario of an array of scenarions, 
% we have the following functions
ne0 = aloha_scenario_get(batch_sc, 'ne0');
% where you just have to provide the same (string) of the parameter you
% want to get. The function automatically find the correct field in the
% scenario sub-structures. Handy !
%
% Of course we can make the same with the RC
RC = aloha_scenario_get(batch_sc, 'RC');
% note that RC is a 17 x 8 matrice, because there is 17 ne0 points and 8 modules.
% An average, min or max on modules can be made:
RC_avg = mean(RC,2); % 17 x 1 array
RC_min = min(transpose(RC));
RC_max = max(transpose(RC));

% Let's make a beatiful figure with "confidence" bounds
aloha_plot_figure('RC vs ne0');
    plot(ne0/1e17, RC_avg, ne0/1e17, RC_min, 'k--', ne0/1e17, RC_max, 'k--');
    xlabel('n_{e0} [1e17 m^{-3}]');
    ylabel('RC [%]');
    grid on;
    legend('Average RC vs edge electron density n_{e0}');
    
% note that we could do same with the call to:
% aloha_plot_reflectionCoeff(batch_sc, ne0)


% That all for the moment !



