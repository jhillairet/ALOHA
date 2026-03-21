function scenario = aloha_compute_RC(scenarios)
% ALOHA
%
% compute (or update) the reflexion coefficient of an ALOHA scenario
% 
% INPUT
%  - scenarios: an ALOHA scenario (or array of scenarios)
% 
% OUTPUT
%  - scenarios: an ALOHA scenario (or array of scenarios) with updated results field
%  
% AUTHOR: JH
% LAST UPDATES
%  - 11/01/2012: creation

for id_scen = 1:length(scenarios)

  scenario = scenarios(id_scen);

  % for easier matlab manipulation, load all the fields of the input scenario 'scenario'
  % into matlab workspace
  aloha_scenario_loadIntoWorkspace;


  % incident wave vector on antenna
  a_acces = a_ampl.*exp(+i*a_phase); 

  % incident and reflected wave vector on/from plasma 
  % 
  % These vector depends on the number of modes we used (1TE + nTM)
  % S_plasma is the scattering diffraction matrix computed from the fortran part of the code.
  % S_ant_ij is the scattering matrix of all modules and passive waveguides and for all modes
  % (the S matrix for a module is computed from HFSS or measured)
  % 
  a_plasma = inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21*a_acces;
  b_plasma = S_plasma*a_plasma;

  % reflexion coefficient @ the mouth of the antenna
  RC_mouth = 100*abs(b_plasma./a_plasma).^2;

  %  Then we plug output from antenna S matrix to input of grill/plasma S matrix
  %  to obtain the plasma coupled antenna scattering matrix
  S_acces = S_ant_11 + S_ant_12*S_plasma*inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21; 

  %  reflected wave vector from antenna
  b_acces = S_acces*a_acces; 

  % power reflexion coefficient @ input of a module
  CoeffRefPuiss = 100*abs(b_acces./a_acces).^2;
  RC=CoeffRefPuiss;

  %% Rajout du 16/05/2007 par Izacard Olivier :
  %'moyenne du coefficient de reflexion en puissance (en %)'
  %100*mean( abs(b_acces./a_acces).^2 )



  % update the scenario fields
  scenario(id_scen).results.a_acces = a_acces;
  scenario(id_scen).results.b_acces = b_acces;
  scenario(id_scen).results.S_acces = S_acces;

  scenario(id_scen).results.a_plasma = a_plasma;
  scenario(id_scen).results.b_plasma = b_plasma;

  scenario(id_scen).results.RC = RC;
  scenario(id_scen).results.CoeffRefPuiss = CoeffRefPuiss;
end % id_scen loop





