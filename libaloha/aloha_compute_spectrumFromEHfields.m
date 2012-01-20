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
f = aloha_scenario_get(scenario, 'freq');

% get the number of poloidal waveguides
nbPolWg = aloha_scenario_get(scenario, 'nwm_theta');

% calculates the spectrum from the EM fields of each poloidal row of waveguides
for id_row = 1:nbPolWg
  % retrieve EM field and spatial toroidal coordinate
  s = scenario.results.abs_z(id_row,:);
  E = scenario.results.E_mouth(3,:,id_row); % Ez
  H = scenario.results.H_mouth(2,:,id_row); % Hy
  
  % waveguide height
  a = scenario.antenna_lh.setup.modules.waveguides.hw_theta;

  % constants
  c = 299792458;
  k0 = 2*pi*f/c;
  lambda0 = c/f;
  Y0 = 1/(120*pi);

  % Spectrum length
  B = 2^20; % (nextpow2(length(E))+3)

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
  
  % Since ALOHA-1D works in 2D, there is additional renormalization per the waveguide heigh and sqrt(2) term for harmonic convention
  dP_nz_row(id_row,:) = a/sqrt(2)* (ds)^2/lambda0 * (1/2*Efft.*conj(-1*Hfft)); % -1 before H because of the cross product sign

  if scenario.options.type_swan_aloha == 0 % SWAN option renormalization
      dP_nz_row(id_row,:) = dP_nz*sqrt(2);
  end
  
  % Power conservation checking
  dnz=nz(2)-nz(1);
  disp(['Power conservation checking : transmited power for row#',num2str(id_row), '[W] :', num2str(dnz*sum(real(dP_nz_row(id_row,:))))]);

end % poloidal rows

% The total spectrum is the sum of the spectra of all rows
dP_nz = sum(dP_nz_row,1);

% Power conservation checking
disp(['Power conservation checking : total transmited power [W] :', num2str(dnz*sum(real(dP_nz)))])
  
% save results to the scenario  
scenario.results.nz = nz;
scenario.results.dP_nz = dP_nz;

