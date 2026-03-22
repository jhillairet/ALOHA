function dP_all = aloha_spectrumAllAntenna(f, ny, dP, delta_h)
%  Compute the global spectrum of the antenna, 
%  which consists to use the spectrum of one half-antenna.
%  
%  INPUTS
%   - f  : frequency [Hz]
%   - ny : Refractive index in the poloidal direction (y) 
%   - dP : 2D spectrum of one-half of an antenna
%   - delta_h : distance between the center of each antenna
%  
%  OUPUT
%   - dP_all : 2D spectrum of the global antenna.
%  
%  AUTHOR : JH
%  
%  LAST UPDATE 
%  - 06/05/2009 add some comments and explanation in the code
%  - 29/07/2008 creation
%  


if ~exist('celerite', 'var')
    aloha_constants
end

k0 = 2*pi*f/celerite;

% Since an antenna is the superposition of same motif (here, only in poloidal direction), 
% the whole spectrum consist in the fourier transform of a convolution of 2 shifted Efields. 
% This lead to multiply the spectrum by a sum of 2 dirac :
CosDeltaH = transpose(4*cos(k0*ny*delta_h/2).^2) * ones(1,size(dP,2));

dP_all = dP .* CosDeltaH;
