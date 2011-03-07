function scenario=aloha_compute_champPlasma1D(scenario)
% compute the electrid and magnetic fields into the plasma region
% 
% 
% AUTHORS: D.Voyer/O.Izacard/J.Hillairet
% LAST CHANGE:  
% - 04/2009: scenario

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


% plasma region coordinates
z_coord = z_coord_min:(z_coord_max-z_coord_min)/nbre_z_coord:z_coord_max;
x_coord = 0:x_coord_max/nbre_x_coord:x_coord_max;
% prevent for empty x_coord, which can occurs when user whant just 1 coord on x=0
if isempty(x_coord)
  x_coord = 0;
end
  

nz_fig_plasma = -100:pas_nz_fig_plasma:100;

poids_E = rac_Zhe*(a_plasma + b_plasma);
poids_H = inv(rac_Zhe)*(a_plasma - b_plasma);
     
poids_E = poids_E*sqrt(2);
poids_H = poids_H*sqrt(2);

nz = nz_fig_plasma;
   
nz(find(nz == 0)) = 1e-4;
nz(find(nz == 1)) = 1-1e-4;
nz(find(nz == -1)) = -1+1e-4;

aloha_constants
   
if (bool_lignes_identiques)
  ne0 = ones(1,nb_g_pol)*ne0;
  dne0 = ones(1,nb_g_pol)*dne0;
  % Rajout le 16/05/2007 par Izacard Olivier
  dne1 = ones(1,nb_g_pol)*dne1;
  d_couche = ones(1,nb_g_pol)*d_couche;
end

k0 = 2*pi*freq/celerite;
nc = ((2*pi*freq)^2)*me*Eps0/(qe^2);
X0 = ne0./nc;
D0 = k0*ne0./dne0;

ind = lig_fig_plasma;% indice de la ligne
   
Ez_x_z = [];
Hy_x_z = [];

version=scenario.plasma.version; % version could be a intrinsic matlab function !     
if (version == 3) | (version == 4)  
   modes_E_sp = [];

   neta_0 = ((nz.*nz - 1).^(1/3))*((D0(ind)/X0(ind))^(2/3))*(X0(ind) - 1)*exp(-i*pi/3);
   Ys = -((airy(1,neta_0)./airy(0,neta_0))./((nz.*nz-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
   
   
   for k = 1:nb_g_total_ligne
  
      % coeff. de normalisation en mode TE des fct. de forme sin et cos
       
      coeff_h = - sqrt(2/(a*b(k)));
      m = 1;
      n = 0;
      fz = -i*(nz./k0).*((1 - ((-1.)^n).*exp(-i*k0*nz*b(k)))./(nz.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nz*z(k));

      modes_E_sp(1+(k-1)*(Nmh+Nme),:) = coeff_h.*fz;

      for n = 1:Nme
      
          m = 1;

          % coeff. de normalisation en mode TM des fct. de forme sin et cos

          coeff_e = -(2*n./sqrt(b(k)/a+(n.^2)*a./b(k)))./b(k);
          fz = -i*(nz./k0).*((1 - ((-1.)^n).*exp(-i*k0*nz*b(k)))./(nz.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nz*z(k));

          modes_E_sp(1+n+(k-1)*(Nmh+Nme),:)= coeff_e.*fz;

      end
   end
  
  
   E_spc = modes_E_sp.*(poids_E(1+nb_g_total_ligne*(Nme+Nmh)*(ind-1):nb_g_total_ligne*(Nme+Nmh)*ind)*ones(1,length(nz)));
   
   
   Ez_0 = sum(E_spc);
   
   
   for ind_x = 1:length(x_coord)
      
      neta = ((nz.*nz - 1).^(1/3))*((X0(ind)/D0(ind))^(1/3))*(k0*x_coord(ind_x) + (D0(ind)/X0(ind))*(X0(ind) - 1))*exp(-i*pi/3);
   
      tampon = (Ez_0.*airy(0,neta)./airy(0,neta_0)).'*ones(1,length(z_coord));   
      Ez_x_z = [Ez_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)];
   
      tampon = ((Y0./((nz.*nz-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6).*...
              (Ez_0.*airy(1,neta)./airy(0,neta_0))).'*ones(1,length(z_coord)); 
              
      Hy_x_z = [Hy_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)]; 


   end
    
   
elseif (version == 6)
   
   x_coord_2nd_couche=x_coord;
       
   modes_E_sp = [];
   modes_H_sp = [];


   neta_1 = ((nz.*nz - 1).^(1/3))*((D1(ind)/X1(ind))^(2/3))*(X1(ind) - 1)*exp(-i*pi/3);
   Y_L = -((airy(1,neta_1)./airy(0,neta_1))./((nz.*nz-1).^(2/3)))*((X1(ind)/D1(ind))^(1/3))*exp(i*pi/6);
       
   neta_0 = ((nz.*nz - 1).^(1/3))*((D0(ind)/X0(ind))^(2/3))*(X0(ind) - 1)*exp(-i*pi/3);
   Ya = -((airy(1,neta_0)./airy(0,neta_0))./((nz.*nz-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
   Yb = -((airy(3,neta_0)./airy(2,neta_0))./((nz.*nz-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6);
       
   neta_d = neta_0 + ((nz.*nz - 1).^(1/3))*((X0(ind)/D0(ind))^(1/3))*k0*d_couche(ind)*exp(-i*pi/3); % Modif de 16/05/2007 'd_couche(ind)' au lieu de 'd_couche'
       
   num = Y_L.*(Yb.*airy(0,neta_d)./airy(0,neta_0) - Ya.*airy(2,neta_d)./airy(2,neta_0)) + ...
            Ya.*Yb.*(airy(3,neta_d)./airy(3,neta_0) - airy(1,neta_d)./airy(1,neta_0));
   
   denom = Y_L.*(airy(0,neta_d)./airy(0,neta_0) - airy(2,neta_d)./airy(2,neta_0)) - ...
            Ya.*airy(1,neta_d)./airy(1,neta_0) + Yb.*airy(3,neta_d)./airy(3,neta_0);
            
   Ys = num./denom;
   
   for k = 1:nb_g_total_ligne
  
      % coeff. de normalisation en mode TE des fct. de forme sin et cos
       
      coeff_h = - sqrt(2/(a*b(k)));

      m = 1;
      n = 0;

      fz = -i*(nz./k0).*((1 - ((-1.)^n).*exp(-i*k0*nz*b(k)))./(nz.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nz*z(k));

      modes_E_sp(1+(k-1)*(Nmh+Nme),:) = coeff_h.*fz;
      modes_H_sp(1+(k-1)*(Nmh+Nme),:) = -Y0*Ys.*coeff_h.*fz;


      for n = 1:Nme
      
          m = 1;

          % coeff. de normalisation en mode TM des fct. de forme sin et cos

          coeff_e = -(2*n./sqrt(b(k)/a+(n.^2)*a./b(k)))./b(k);
          
          fz = -i*(nz./k0).*((1 - ((-1.)^n).*exp(-i*k0*nz*b(k)))./(nz.^2 - (n*pi/(k0*b(k)))^2)).*exp(-i*k0*nz*z(k));

          modes_E_sp(1+n+(k-1)*(Nmh+Nme),:)= coeff_e.*fz;
          modes_H_sp(1+n+(k-1)*(Nmh+Nme),:)= -Y0*Ys.*coeff_e.*fz;

      end
   end
  
   idx_g = [1+nb_g_total_ligne*(Nme+Nmh)*(ind-1):nb_g_total_ligne*(Nme+Nmh)*ind];
   
   E_spc = modes_E_sp.*(poids_E(idx_g)*ones(1,length(nz)));
   H_spc = modes_H_sp.*(poids_E(idx_g)*ones(1,length(nz)));

   % couche n1
   
   Ez_0 = sum(E_spc);
   Hy_0 = sum(H_spc);
   
   
   Va = (Hy_0 + Y0*Yb.*Ez_0)./(Y0*(Yb-Ya));
   Vb = (Hy_0 + Y0*Ya.*Ez_0)./(Y0*(Ya-Yb));
   
   x_coord = (0:d_couche(1)/6:d_couche(1)-d_couche(1)/6); % Modif de 16/05/2007 'd_couche(1)' au lieu de 'd_couche'
   x_tot = x_coord;
   
   
   for ind_x = 1:length(x_coord)
      
      neta = ((nz.*nz - 1).^(1/3))*((X0(ind)/D0(ind))^(1/3))*(k0*x_coord(ind_x) + (D0(ind)/X0(ind))*(X0(ind) - 1))*exp(-i*pi/3);
   
      tampon = (Va.*airy(0,neta)./airy(0,neta_0)+Vb.*airy(2,neta)./airy(2,neta_0)).'*ones(1,length(z_coord));   
      Ez_x_z = [Ez_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)];
   
      tampon = ((Y0./((nz.*nz-1).^(2/3)))*((X0(ind)/D0(ind))^(1/3))*exp(i*pi/6).*...
                  (Va.*airy(1,neta)./airy(0,neta_0)+Vb.*airy(3,neta)./airy(2,neta_0))).'*ones(1,length(z_coord)); 
                  
      Hy_x_z = [Hy_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)]; 
   
   end
   
   % couche n2
   
   Vc = Va.*airy(0,neta_d)./airy(0,neta_0)+Vb.*airy(2,neta_d)./airy(2,neta_0);
   
   x_coord = x_coord_2nd_couche;
   x_tot= [x_tot,x_coord+d_couche(2)];    % Modif de 16/05/2007 'd_couche(2)' au lieu de 'd_couche'
   
   for ind_x = 1:length(x_coord)
      
      neta = ((nz.*nz - 1).^(1/3))*((X1(ind)/D1(ind))^(1/3))*(k0*x_coord(ind_x) + (D1(ind)/X1(ind))*(X1(ind) - 1))*exp(-i*pi/3);
   
      tampon = (Vc.*airy(0,neta)./airy(0,neta_1)).'*ones(1,length(z_coord));   
      Ez_x_z = [Ez_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)];
   
      tampon = ((Y0./((nz.*nz-1).^(2/3)))*((X1(ind)/D1(ind))^(1/3))*exp(i*pi/6).*...
                (Vc.*airy(1,neta)./airy(0,neta_1))).'*ones(1,length(z_coord)); 
                
      Hy_x_z = [Hy_x_z;sum((exp(i*k0*nz'*z_coord).*tampon)).*k0*pas_nz_fig_plasma./(2*pi)]; 
   end
   
   x_coord=x_tot;
end
 

% save resukts to the scenario
scenario.results = aloha_setfield(scenario.results, x_coord, z_coord, Ez_x_z, Hy_x_z);
   
