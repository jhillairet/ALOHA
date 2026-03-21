function scenario=aloha_compute_spectrum2D(scenario, varargin)
%  Compute the power spectrum of the antenna (aloha 2D)
%  
%  INPUT
%   - ALOHA scenario
%   - [optionnal, string] file (with path) containing the spectrum
%  OUPUT
%   - ALOHA scenario with suppl. fields into the sub-field 'results"
%     corresponding to the spectrum parameters (dP, dP_nz, nz, ny)
%  
if nargin == 1
    ascii_file_path = [scenario.options.aloha_path,scenario.options.chemin_binaire_fortran, '/ALOHA2D.out.spectralFields.dat'];
elseif nargin == 2
    ascii_file_path = varargin{1};
end

aloha_constants; % celerite, Y0

aloha_message(['Reading result file:', 'ALOHA2D.out.spectralFields.dat...']);
fid = fopen(ascii_file_path,'r');
    % read nb of mode and nb of guides
    C = textscan(fid, '%f %f', 1, 'headerLines', 2);
    nb_modes = C{1};
    nb_guides= C{2};
    % read nb of ny & nz
    C = textscan(fid, '%f %f', 1, 'headerLines', 2);
    ny_nb = C{1};
    nz_nb = C{2};
    % read spectral electric and magnetic field components
    C = textscan(fid, '%f %f %f %f %f %f %f %f', 'headerLines', 2);
    Ey = complex(C{1}, C{2});
    Ez = complex(C{3}, C{4});
    Hy = complex(C{5}, C{6});
    Hz = complex(C{7}, C{8});
fclose(fid);

Ey_ny_nz = reshape(Ey, ny_nb*nz_nb, nb_modes*nb_guides).';
Ez_ny_nz = reshape(Ez, ny_nb*nz_nb, nb_modes*nb_guides).';
Hy_ny_nz = reshape(Hy, ny_nb*nz_nb, nb_modes*nb_guides).';
Hz_ny_nz = reshape(Hz, ny_nb*nz_nb, nb_modes*nb_guides).';

% % reading admittance file
% admittance_file = 'ALOHA2D.out.admittance.dat';
% fid = fopen(admittance_file, 'r');
%     C = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', 'headerLines', 1);
% fclose(fid);
% 
% Y_nz = C{1} ;
% Y_ny = C{2} ;
% Y_yy = complex(C{3}, C{4});
% Y_yz = complex(C{5}, C{6});
% Y_zy = complex(C{7}, C{8});
% Y_zz = complex(C{9}, C{9});
% 
% Hy_ny_nz = Y0*ones(nb_modes*nb_guides, 1)*(Y_yy.').*Ey_ny_nz + ...
%            Y0*ones(nb_modes*nb_guides, 1)*(Y_yz.').*Ez_ny_nz;
% Hz_ny_nz = Y0*ones(nb_modes*nb_guides, 1)*(Y_zy.').*Ey_ny_nz + ...
%            Y0*ones(nb_modes*nb_guides, 1)*(Y_zz.').*Ez_ny_nz;
       
poids_E = scenario.results.rac_Zhe*(scenario.results.a_plasma + scenario.results.b_plasma);
poids_H = inv(scenario.results.rac_Zhe)*(scenario.results.a_plasma - scenario.results.b_plasma);
    
dny = scenario.options.dny;%(scenario.options.ny_max-scenario.options.ny_min)/scenario.options.ny_nb;
dnz = scenario.options.dnz;%(scenario.options.nz_max-scenario.options.nz_min)/scenario.options.nz_nb;
ny = [scenario.options.ny_min:dny:scenario.options.ny_max-dny]; % -dny to get the correct number of points
nz = [scenario.options.nz_min:dnz:scenario.options.nz_max-dnz];

Eyt_ny_nz = Ey_ny_nz.*(poids_E*ones(1, ny_nb*nz_nb));
Ezt_ny_nz = Ez_ny_nz.*(poids_E*ones(1, ny_nb*nz_nb));
Hyt_ny_nz = Hy_ny_nz.*(poids_H*ones(1, ny_nb*nz_nb));
Hzt_ny_nz = Hz_ny_nz.*(poids_H*ones(1, ny_nb*nz_nb));
   
k0 = 2*pi*scenario.antenna.freq/celerite;

dP_ligne = ((k0/(2*pi))^2)*( ...
     sum(Eyt_ny_nz).*sum(conj(Hzt_ny_nz)) ...
    -sum(Ezt_ny_nz).*sum(conj(Hyt_ny_nz)) );

dP = reshape(dP_ligne, ny_nb, nz_nb);

dP_nz = dny*sum(real(dP));

% display the main peak value 
[max_nz, id_max_nz] = max(abs(dP_nz));
nz0 = nz(id_max_nz);
disp(aloha_message(['Main n_z peak ("n_z0") : ', num2str(nz0)]));  

% save results into the scenario
scenario.results = aloha_setfield([scenario.results], ny, nz, nz0, dny, dnz, dP, dP_nz); 
