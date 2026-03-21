% ideal waveguide scattering matrix for one waveguide.
% 
% This matrix consist in the matrix of a 
% perfect waveguide without transmission losses 
% nor reflexion 
S = [0,1; ...
     1,0];

S = reshape(S,1,4);
