function scenario = aloha_compute_spectrumFromEHfield(scenario)
% Calculate the power spectrum p(n//) from the Electric and Magnetic field
% along the launcher mouth
%
% INPUT
%  - scenario : ALOHA scenario (containing the coupling calculations)
%  
% OUTPUT
%  - scenario : ALOHA scenario with the following fields in addition in the 'results' field:
%       - nz : n parallel
%       - dP_nz : power spectrum (complex valued)
%
% The power spectrum density p is normalized in order to satisfy
% the continuous Parseval theorem, ie:
% 
% integral( 1/2|ExH*|^2, s) = integral(dP_nz, n_parallel)
%
%
% AUTHOR: JH,MP
% LAST UPDATE: 16/03/2011

% get parameters from the scenario
f = scenario.antenna.freq;
s = scenario.results.abs_z(1,:);
E = scenario.results.E_mouth(3,:); % Ez
H = scenario.results.H_mouth(2,:); % Hy


% constants
c = 299792458;
k0 = 2*pi*f/c;
lambda0 = c/f;
Y0 = 1/(120*pi);

% Spectrum length
B = 2^(nextpow2(length(E))+1);

% E and H field fft (without 1/length normalization)
Efft = fftshift(fft(E,B));
Hfft = fftshift(fft(H,B));


%% (spatial) frequency 
% spectral sampling interval
ds = s(2) - s(1); % assuming ds constant ! 
df = 1/(B*ds);
% spatial frequencies bins as defined in Matlab
K = [-B/2:B/2-1];
% (spatial) frequency range 
Fz = K*df;

%% parallel index
nz =(2*pi/k0)*Fz;
%(length(E)/B)*K*df*(2*pi)/k0;

%% power spectrum
%
% the power spectrum density is normalized in order to satisfy
% the continuous Parseval theorem, ie:
% 
% integral( |E|^2, s) = integral(p, n_parallel)
%

dP_nz = (ds)^2/lambda0 * (1/2*Efft.*conj(-1*Hfft)); % -1 before H because of the cross product sign

if scenario.options.type_swan_aloha == 0 % SWAN option renormalization
    dP_nz = dP_nz*sqrt(2);
end

scenario.results.nz = nz;
scenario.results.dP_nz = dP_nz;

