function scenario=aloha_compute_spectrum1D(scenario)
%  Compute the power spectrum of the antenna
%  
%  INPUT
%   - ALOHA scenario
%  OUPUT
%   - ALOHA scenario with suppl. fields into the sub-field 'results"
%     corresponding to the spectrum parameters (dP, dP_nz, nz, ny)
%  
%  NB : All the options such as : 
%    dny, dnz, ny_min, ny_max, nz_min, nz_max
%  should be defined into the 'option' sub-field of the input scenario
%  
%  LAST UPDATE:
%   - 11/01:2012: now able to calculate the spectrum of a scenario which input excitation has been changed 
%                (especially usefull for parametric scan on module excitation)  
%   - 04/2009 : integration of the code to the scenario formalism (input/output)

% J.Hillairet - 11/01/2012
% 
% When one wants to make a parameter study on the module's phase excitation
% he should not need to recalculate the plasma scattering matrix, which does not depends on the excitation.
% Thus, the a_plasma, b_plasma are recalculated below (in the aloha_compute_RC function), just in case, if we are in the situation described upper.
scenario = aloha_compute_RC(scenario);

% for easier matlab manipulation, load all the fields of the input scenario 'scenario'
% into matlab workspace
aloha_scenario_loadIntoWorkspace;

% test if the main results, such as the plasma scattering matrix, 
% had been calculated. If not, warn user.
if not(exist('S_plasma'))
  error('The plasma scattering matrix has not been calculated !');
end

% Load into workspace the geometrical parameter of the scenario architecture 
% and get the main geometrical parameters.
if aloha_isAntennaITM(scenario)
      disp(aloha_message('assuming ITM antenna description')); 
      aloha_utils_ITM2oldAntenna; % convert ITM data to old fashioned ALOHA global parameters
    [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));

else
    eval(architecture);
    [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(architecture);
end
%  aloha_computeSpectrum1D
%  
%  

    
    % now calculates the weight of electric and magnetic field from the power waves
    poids_E = rac_Zhe*(a_plasma + b_plasma);
    poids_H = inv(rac_Zhe)*(a_plasma - b_plasma);
    
    poids_E = poids_E*sqrt(2);
    poids_H = poids_H*sqrt(2);
      
   ny = ny_min:dny:ny_max;
   nz = nz_min:dnz:nz_max;
   
   % avoid pure 0 and +/-1 discontinuities by adding a small real part
   nz(find(nz == 0)) = 1e-4;
   nz(find(nz == 1)) = 1-1e-4;
   nz(find(nz == -1)) = -1+1e-4;

   nyt = kron(ones(1,length(nz)),ny);
   nzt = kron(nz,ones(1,length(ny)));
   
   aloha_constants
   
   % if the plasma in front of each rows is identical,
   % we duplicate plasma rows parameters
   if (bool_lignes_identiques)
      ne0 = ones(1,nb_g_pol)*ne0;
      dne0 = ones(1,nb_g_pol)*dne0;
      % Rajout le 16/05/2007 par Izacard Olivier
      dne1 = ones(1,nb_g_pol)*dne1;
      d_couche = ones(1,nb_g_pol)*d_couche;
      d_vide = ones(1,nb_g_pol)*d_vide;
      % MODIF JH 05/2009
      % Si les lignes sont identiques, 
      % Il n'y a pas besoin de faire le meme calcul plusieurs fois.
      % Donc on force le calcul pour qu'il se fasse qu'une seule fois.

   end
   
   k0 = 2*pi*freq/celerite;
   nc = ((2*pi*freq)^2)*me*Eps0/(qe^2);
   X0 = ne0./nc;
   D0 = k0*ne0./dne0;
    
   % for version == 6 only
   n1=ne0+dne0.*d_couche;
   X1=n1/nc;
   D1=k0*n1./dne1;

   % calcul la matrice des spectres E et H des modes

   modes_E_sp = [];
   modes_H_sp = [];
   
   version=scenario.plasma.version; % version could be a intrinsic matlab function !
   if (version == 3) | (version == 4)
     if exist('progress')
        progress('init');
     end    
     for ind = 1:nb_g_pol
    disp(['Ligne pol. ', num2str(ind)]);

     neta_0 = ((nzt.*nzt - 1).^(1/3))*((D0(ind)/X0(ind))^(2/3))*(X0(ind) - 1)*exp(-i*pi/3);
     Ys = -((airy(1,neta_0)./airy(0,neta_0))./((nzt.*nzt-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
   
       for k = 1:nb_g_total_ligne
  
          % coeff. de normalisation en mode TE des fct. de forme sin et cos
       
          coeff_h = - sqrt(2/(a*b(k)));

          m = 1;
          n = 0;

          fy = -(m*pi/(a*(k0^2)))*((1 - ((-1)^m)*exp(-i*k0.*nyt*a))./(nyt.^2 - (m*pi/(k0*a))^2)).*exp(-i*k0*nyt*y(ind));
          fz = -i*(nzt./k0).*((1 - ((-1.)^n).*exp(-i*k0*nzt*b(k)))./(nzt.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nzt*z(k));

          modes_E_sp(1+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:) = coeff_h.*fy.*fz;
          modes_H_sp(1+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:) = -Y0*Ys.*coeff_h.*fy.*fz;


          for n = 1:Nme
      
             m = 1;

             % coeff. de normalisation en mode TM des fct. de forme sin et cos

             coeff_e = -(2*n./sqrt(b(k)/a+(n.^2)*a./b(k)))./b(k);
	 
	 
	         fy = -(m*pi/(a*(k0^2)))*((1 - ((-1)^m)*exp(-i*k0.*nyt*a))./(nyt.^2 - (m*pi/(k0*a))^2)).*exp(-i*k0*nyt*y(ind));
                 fz = -i*(nzt./k0).*((1 - ((-1.)^n).*exp(-i*k0*nzt*b(k)))./(nzt.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nzt*z(k));

                 modes_E_sp(1+n+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:)= coeff_e.*fy.*fz;
	         modes_H_sp(1+n+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:)= -Y0*Ys.*coeff_e.*fy.*fz;

           end % for n = 1:Nme
		
        end % for k = 1:nb_g_total_ligne
     
     end % for ind = 1:nb_g_pol
     
     
   elseif (version == 6)
     
     for ind = 1:nb_g_pol

       neta_1 = ((nzt.*nzt - 1).^(1/3))*((D1(ind)/X1(ind))^(2/3))*(X1(ind) - 1)*exp(-i*pi/3);
       Y_L = -((airy(1,neta_1)./airy(0,neta_1))./((nzt.*nzt-1).^(2/3)))*((X1(ind)/D1(ind))^(1/3))*exp(i*pi/6);
       
       neta_0 = ((nzt.*nzt - 1).^(1/3))*((D0(ind)/X0(ind))^(2/3))*(X0(ind) - 1)*exp(-i*pi/3);
       Ya = -((airy(1,neta_0)./airy(0,neta_0))./((nzt.*nzt-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
       Yb = -((airy(3,neta_0)./airy(2,neta_0))./((nzt.*nzt-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
       
       neta_d = neta_0 + ((nzt.*nzt - 1).^(1/3))*((X0(ind)/D0(ind))^(1/3))*k0*d_couche(ind)*exp(-i*pi/3); % Modif de 16/05/2007 'd_couche(ind)' au lieu de 'd_couche'
       
       num = Y_L.*(Yb.*airy(0,neta_d)./airy(0,neta_0) - Ya.*airy(2,neta_d)./airy(2,neta_0)) + ...
       				Ya.*Yb.*(airy(3,neta_d)./airy(3,neta_0) - airy(1,neta_d)./airy(1,neta_0));
   
       denom = Y_L.*(airy(0,neta_d)./airy(0,neta_0) - airy(2,neta_d)./airy(2,neta_0)) - ...
       				Ya.*airy(1,neta_d)./airy(1,neta_0) + Yb.*airy(3,neta_d)./airy(3,neta_0);
				
       Ys = num./denom;
       
       %% admittance propagation trough the vacuum gap of width d_vide
       alph = 1;%-i*pertes; % BUG?"pertes" is not always in the matlab memory?
       gamm = k0*(nzt.^2-alph).^(1/2); 
       u = gamm*d_vide(ind);
       % Hyperbolic tangent calculation: Taylor expansion (source: wikipedia) 
       %tanh_d_vide = u - (u.^3)/3 + 2*(u.^5)/15 - 17*(u.^7)/315; 
       tanh_d_vide = tanh(u);

       y2 = i*k0./gamm;
       Ys = y2.*(Ys+y2.*tanh_d_vide)./(Ys.*tanh_d_vide+y2);

       
   
       for k = 1:nb_g_total_ligne
  
          % coeff. de normalisation en mode TE des fct. de forme sin et cos
       
          coeff_h = - sqrt(2/(a*b(k)));

          m = 1;
          n = 0;

          fy = -(m*pi/(a*(k0^2)))*((1 - ((-1)^m)*exp(-i*k0.*nyt*a))./(nyt.^2 - (m*pi/(k0*a))^2)).*exp(-i*k0*nyt*y(ind));
          fz = -i*(nzt./k0).*((1 - ((-1.)^n).*exp(-i*k0*nzt*b(k)))./(nzt.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nzt*z(k));

          modes_E_sp(1+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:) = coeff_h.*fy.*fz;
          modes_H_sp(1+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:) = -Y0*Ys.*coeff_h.*fy.*fz;


          for n = 1:Nme
      
             m = 1;

             % coeff. de normalisation en mode TM des fct. de forme sin et cos

             coeff_e = -(2*n./sqrt(b(k)/a+(n.^2)*a./b(k)))./b(k);
	 
	 
	         fy = -(m*pi/(a*(k0^2)))*((1 - ((-1)^m)*exp(-i*k0.*nyt*a))./(nyt.^2 - (m*pi/(k0*a))^2)).*exp(-i*k0*nyt*y(ind));
                 fz = -i*(nzt./k0).*((1 - ((-1.)^n).*exp(-i*k0*nzt*b(k)))./(nzt.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nzt*z(k));

                 modes_E_sp(1+n+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:)= coeff_e.*fy.*fz;
	         modes_H_sp(1+n+(k-1)*(Nmh+Nme)+nb_g_total_ligne*(Nme+Nmh)*(ind-1),:)= -Y0*Ys.*coeff_e.*fy.*fz;

          end % for n = 1:Nme

       end % for k = 1:nb_g_total_ligne
     
     end % for ind = 1:nb_g_pol
    
   else
      error('Which version ??')
     
   end
  
   E_spc = modes_E_sp.*(poids_E*ones(1,length(nyt)));
   H_spc = modes_H_sp.*(poids_E*ones(1,length(nyt)));

   % calcul de la densite de puissance
   %dP_ligne = -((k0/(2*pi))^2)*sum(E_spc).*sum(((H_spc.')'));
   %dP = reshape(dP_ligne,length(ny),length(nz));
    
   % MAJ : J.Hillairet 19/12/07
   % optimisation calcul (gain vitesse * 2) 
    % MAJ JH Jullet 2009 :
    % Il y a un facteur 2 de trop ! 
    % On divise donc le spectre par 2.
    % cela reviend certainement ?? supprimer les racines de 2 devant les poids des champs
    dP = -1/2*(k0/2/pi)^2 .* reshape((sum(E_spc,1).*conj(sum(H_spc,1))), length(ny), length(nz));
    dP_nz = dny*sum(dP);

    % display the main peak value 
    [max_nz, id_max_nz] = max(real(dP_nz));
    nz0 = nz(id_max_nz);
    disp(aloha_message(['Main n_z peak ("n_z0") : ', num2str(nz0)]));  

    
    
    
    
    
% save results into the scenario
scenario.results = aloha_setfield([scenario.results], ny, nz, nz0, dny, dnz, dP, dP_nz); 
