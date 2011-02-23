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
    ascii_file_path = [scenario.options.aloha_path,scenario.options.chemin_binaire_fortran, '/Spect_plasma2.dat'];
elseif nargin == 2
    ascii_file_path = varargin{1};
end

aloha_message(['Reading result file:', 'Spect_plasma2.dat...'])
fid = fopen(ascii_file_path,'r');
    % read nb of mode and nb of guides
    C = textscan(fid, '%f %f', 1, 'headerLines', 2);
    nbre_modes = C{1};
    nbre_guides= C{2};
    % read nb of ny & nz
    C = textscan(fid, '%f %f', 1, 'headerLines', 2);
    nbre_ny = C{1};
    nbre_nz = C{2};
    % read spectral electric and magnetic field components
    C = textscan(fid, '(%f,%f) (%f,%f) (%f,%f) (%f,%f)', 'headerLines', 2);
    Ey_r = C{1};Ey_i = C{2};
    Ez_r = C{3};Ez_i = C{4};
    Hy_r = C{5};Hy_i = C{6};
    Hz_r = C{7};Hz_i = C{8};
fclose(fid);

%  nb = nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1)
%  
%  fname =  [scenario.options.aloha_path,scenario.options.chemin_binaire_fortran, '/Spect_plasma.dat'];
%  
%  byteorder = 'ieee-be';            % IEEE Big-Endian format  (use "ieee-le" for Little-Endian)
%  [fid,msg] = fopen(fname, 'r', byteorder) ;
%  % read the trailing number of records...
%  [Ey_r,count] = fread(fid, nb,'float'); 
%  [Ey_i,count] = fread(fid, nb,'float'); 
%  [Ez_r,count] = fread(fid, nb,'float');
%  [Ez_i,count] = fread(fid, nb,'float');
%  [Hy_r,count] = fread(fid, nb,'float'); 
%  [Hy_i,count] = fread(fid, nb,'float'); 
%  [Hz_r,count] = fread(fid, nb,'float');
%  [Hz_i,count] = fread(fid, nb,'float');
%  fclose(fid);



% JH 06/10/2010
% Il y'a un bug dans la taille du vecteur lu. Je pense que le 
% bug est lie au fait que les tailles des vecteurs sont fixes
% dans le code fortran, ce qui genere un gag quelque part. 
% 
% Si la taille des vecteurs n'est pas celle que l'on attend,
% on bourre par des zeros :
array_size = round((nbre_ny+1)*(nbre_nz+1)*nbre_modes*nbre_guides);

%  if length(Ey_r) ~= array_size 
%      Ey_r = [Ey_r; zeros(array_size-length(Ey_r),1)];
%      Ey_i = [Ey_i; zeros(array_size-length(Ey_i),1)];
%      Ez_r = [Ez_r; zeros(array_size-length(Ez_r),1)];
%      Ez_i = [Ez_i; zeros(array_size-length(Ez_i),1)];
%      Hy_r = [Hy_r; zeros(array_size-length(Hy_r),1)];
%      Hy_i = [Hy_i; zeros(array_size-length(Hy_i),1)];
%      Hz_r = [Hz_r; zeros(array_size-length(Hz_r),1)];
%      Hz_i = [Hz_i; zeros(array_size-length(Hz_i),1)];
%  end


Ey_ny_nz=reshape(complex(Ey_r, Ey_i), (nbre_ny+1)*(nbre_nz+1), nbre_modes*nbre_guides).';
Ez_ny_nz=reshape(complex(Ez_r, Ez_i), (nbre_ny+1)*(nbre_nz+1), nbre_modes*nbre_guides).';
Hy_ny_nz=reshape(complex(Hy_r, Hy_i), (nbre_ny+1)*(nbre_nz+1), nbre_modes*nbre_guides).';
Hz_ny_nz=reshape(complex(Hz_r, Hz_i), (nbre_ny+1)*(nbre_nz+1), nbre_modes*nbre_guides).';

%  %   BINARY FILE
%  ascii_file_path = [scenario.options.aloha_path,scenario.options.chemin_binaire_fortran, '/Spect_plasma.dat'];
%  
%  fid = fopen(ascii_file_path,'r');
%      fread(fid,1,'int32');   
%      keyboard
%      Ey_ny_nz=[];
%      tamp = fread(fid,2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1),'double');
%      Ey_ny_nz = tamp(1:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1)) + i*tamp(2:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1));
%      Ey_ny_nz=reshape(Ey_ny_nz,nbre_modes*nbre_guides,(nbre_ny+1)*(nbre_nz+1));
%      
%       
%  
%      fread(fid,1,'double'); 
%      
%      Ez_ny_nz=[];
%      tamp = fread(fid,2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1),'double');
%      Ez_ny_nz = tamp(1:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1)) + i*tamp(2:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1));
%      Ez_ny_nz=reshape(Ez_ny_nz,nbre_modes*nbre_guides,(nbre_ny+1)*(nbre_nz+1));
%      
%      fread(fid,1,'double'); 
%      
%      Hy_ny_nz=[];
%      tamp = fread(fid,2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1),'double');
%      Hy_ny_nz = tamp(1:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1)) + i*tamp(2:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1));
%      Hy_ny_nz=reshape(Hy_ny_nz,nbre_modes*nbre_guides,(nbre_ny+1)*(nbre_nz+1));
%      
%      fread(fid,1,'double'); 
%      
%      Hz_ny_nz=[];
%      tamp = fread(fid,2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1),'double');
%      Hz_ny_nz = tamp(1:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1)) + i*tamp(2:2:2*nbre_modes*nbre_guides*(nbre_ny+1)*(nbre_nz+1));
%      Hz_ny_nz=reshape(Hz_ny_nz,nbre_modes*nbre_guides,(nbre_ny+1)*(nbre_nz+1));
%      
%  fclose(fid);





poids_E = scenario.results.rac_Zhe*(scenario.results.a_plasma + scenario.results.b_plasma);
poids_H = inv(scenario.results.rac_Zhe)*(scenario.results.a_plasma - scenario.results.b_plasma);
    
dny=(scenario.options.ny_max-scenario.options.ny_min)/scenario.options.nbre_ny;
dnz=(scenario.options.nz_max-scenario.options.nz_min)/scenario.options.nbre_nz;
ny = scenario.options.ny_min:dny:scenario.options.ny_max;
nz = scenario.options.nz_min:dnz:scenario.options.nz_max;

Eyt_ny_nz = Ey_ny_nz.*(poids_E*ones(1,(scenario.options.nbre_ny+1)*(scenario.options.nbre_nz+1)));
Ezt_ny_nz = Ez_ny_nz.*(poids_E*ones(1,(scenario.options.nbre_ny+1)*(scenario.options.nbre_nz+1)));
Hyt_ny_nz = Hy_ny_nz.*(poids_E*ones(1,(scenario.options.nbre_ny+1)*(scenario.options.nbre_nz+1)));
Hzt_ny_nz = Hz_ny_nz.*(poids_E*ones(1,(scenario.options.nbre_ny+1)*(scenario.options.nbre_nz+1)));
   
aloha_constants; % celerite
k0 = 2*pi*scenario.antenna.freq/celerite;

dP_ligne = ((k0/(2*pi))^2)*(sum(Eyt_ny_nz).*sum(((Hzt_ny_nz.')'))-sum(Ezt_ny_nz).*sum(((Hyt_ny_nz.')')));

dP = reshape(dP_ligne,scenario.options.nbre_ny+1, scenario.options.nbre_nz+1);

dP_nz = dny*sum(real(dP));

    % display the main peak value 
    [max_nz, id_max_nz] = max(abs(dP_nz));
    nz0 = nz(id_max_nz);
    disp(aloha_message(['Main n_z peak ("n_z0") : ', num2str(nz0)]));  

% save results into the scenario
scenario.results = aloha_setfield(scenario.results, ny, nz, nz0, dny, dnz, dP, dP_nz); 
