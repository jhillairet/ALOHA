% ALOHA tutorial
% J.Hillairet
% October 2011
%
% Scattering matrix of a perfect (loss free) waveguide : 
%
S = [0, 1; ... 
     1, 0];
 
% For compatibility reasons with previous ALOHA versions and with the HFSS output format, 
% the scattering matrix must be reshaped to 1 line x N columns.
S = reshape(S,1,4);