%  Calcule la reponse de l'antenne a une excitation (Coefficients de relfexion)
%  
%  Determine les coefficients de reflexion a partir 
%  de la matrice S plasma 


% incident wave vector on antenna
a_acces = a_ampl.*exp(i*a_phase);

% incident and reflected wave vector on/from plasma 
% 
% These vector depends on the number of modes we used (1TE + nTM)
% S_plasma is the scattering diffraction matrix computed from the fortran part of the code.
% S_ant_ij is the scattering matrix of all modules and passive waveguides and for all modes
% (the S matrix for a module is computed from HFSS or measured)
% 
a_plasma = inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21*a_acces;
b_plasma = S_plasma*a_plasma;

%  Then we plug output from antenna S matrix to input of grill/plasma S matrix
%  to obtain the plasma coupled antenna scattering matrix
S_acces = S_ant_11 + S_ant_12*S_plasma*inv(eye(length(S_plasma)) - S_ant_22*S_plasma)*S_ant_21; 

%  reflected wave vector from antenna
b_acces = S_acces*a_acces; 

% power reflexion coefficient 
CoeffRefPuiss = 100*abs(b_acces./a_acces).^2;

%% Rajout du 16/05/2007 par Izacard Olivier :
%'moyenne du coefficient de reflexion en puissance (en %)'
%100*mean( abs(b_acces./a_acces).^2 )









