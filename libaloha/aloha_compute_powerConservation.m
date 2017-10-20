function scenario=aloha_compute_powerConservation(scenario)
%  Compute and check the power conservation of the computed results.
%  
%  INPUT
%   - scenario : ALOHA scenario
%  OUPUT
%   - scenario : the same ALOHA scenario with additional fields
%                in the 'results' sub-field, such as :
%  
%  NB : This function needs that the spectrum have been previouly calculated.
%  
%  AUTHOR: J.Hillairet
%  LAST CHANGES: 
%   - 04-2009 : Creation

% for easier matlab manipulation, load all the fields of the input scenario 'scenario'
% into matlab workspace
aloha_scenario_loadIntoWorkspace;

% check if the spectrum have been previously calculated
if not(exist('dP'))
  error('The spectrum must have been calculated. This is not the case here !')
end

disp(aloha_message('Check power conservation'));
% get the reflected power 
R = aloha_scenario_get(scenario, 'CoeffRefPuiss')/100;

% input power into the antenna by module
Pin = abs(aloha_scenario_get(scenario, 'a_ampl')).^2; 
disp(aloha_message([' - Incident power Pin = ', num2str(sum(Pin)), ' W [half-antenna]']));

Pout = abs(aloha_scenario_get(scenario, 'b_acces')).^2; 
disp(aloha_message([' - Reflected power Pout = ', num2str(sum(Pout)), ' W [half-antenna]']));

            
% sum of the transmitted power to the plasma by module
Ptr_plasma_RC = (1-R).*Pin;

disp(aloha_message([' - Transmitted power to plasma Ptr_plasma ([1-R]*Pin) = ', num2str(sum(Ptr_plasma_RC)), ' W [half-antenna]']));

% Spectrum coupled power integral(dP)
Ptr_plasma_spectrum = real(sum(sum(dP))*dny*dnz)   ;
disp(aloha_message([' - Spectrum coupled power Ptr_plasma (nz in [',...
                    num2str(min(nz)),',',num2str(max(nz)),'])= ', num2str(Ptr_plasma_spectrum), ...
                    ' W [half-antenna]']));

% ratio of the two coupled power to check the power conservation
% 
nu = Ptr_plasma_spectrum/sum(Ptr_plasma_RC);
disp(aloha_message(['-> Ratio nu=Ptr_plasma_spectre/[(1-R)P_inc] = ',num2str(nu)]));

% save values into the scenario
scenario.results = aloha_setfield([scenario.results], Pin, Ptr_plasma_RC, Ptr_plasma_spectrum, nu);
