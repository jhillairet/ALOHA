function scenario=aloha_compute_champEmbouchure1D(scenario)
% Compute the parallel component of electric field (Ez) in the  mouth of the antenna.
% Calcule le champ electrique parallele (Ez) dans l'embouchure de l'antenne
% 
% 1D case : the electric field correspond to the modes of a parallel-plate waveguide, ie 1 TEM and n TM modes.
%  

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
eval(architecture);
[b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(architecture);


% Modification JH 04/2009
% Precedemment, l'echantillonnage n'etait pas regulier : les champs dans les guides etaient 
% decrit sur 100 points, tandis que le champs au niveau des septums etaient decrit avec 1 seul point.
% Lorsqu'on faisait la transformee de fourier du champ total, on obtenait alors un resultat "suprenant".
% 
% La modification du code consiste donc a exprimer les champs avec un pas constant et controle.
dz = 1e-4; % ATTENTION : doit etre absolument inferieur a la plus petite decimale des dimensions b ou e


poids_E = sqrt(2)*rac_Zhe*(a_plasma + b_plasma);
%  poids_H = sqrt(2)*inv(rac_Zhe)*(a_plasma - b_plasma);


abs_z = [];
Efield = [];


% index of TM modes
n = [1:Nme]';

% for all poloidal lines
for idx_pol = 1:nb_g_pol 
    abs_z_idx_pol = [];
    fztot = [];
    % for all waveguides in a poloidal line
    for idx_tor = 1:nb_g_total_ligne
        % Coefficients modaux (pas de dependance spatiale : elle est ajoutee plus bas)
        % TE mode(s) (Nmh,1)
        coeff_h = - sqrt(2/(a*b(idx_tor)));
        % TM mode(s) (Nme,1)
        coeff_e = -(2*n./sqrt(b(idx_tor)/a+(n.^2)*a./b(idx_tor)))./b(idx_tor);
        
        % ponderation par les resultats d'ALOHA
        mode_idx=1+(idx_tor-1)*(Nme+Nmh)+(idx_pol-1)*nb_g_total_ligne*(Nme+Nmh):Nme+Nmh+(idx_tor-1)*(Nme+Nmh)+(idx_pol-1)*nb_g_total_ligne*(Nme+Nmh);
        coeff = [coeff_h;coeff_e].*...
                poids_E(mode_idx);
           
        % echantillonage de l'abscisse du guide en cours
        z_g = [0:dz:b(idx_tor)];
    
        % fonctions spatiales modales
        fhz = ones(1,length(z_g));
        n = (1:Nme)';
        fez = cos(pi*n/b(idx_tor)*z_g);

        fz = [fhz;fez];
    
        if (Nme ~= 0)
            fztot = (sum((coeff*ones(1,length(z_g))).*fz));
        else
            fztot = ((coeff*ones(1,length(z_g))).*fz);
        end
       
       
        abs_z_idx_pol = [abs_z_idx_pol, z(idx_tor) + z_g];



        % Le champ electrique transverse est nul au niveau des septums,
        % cad entre deux guides consecutifs.
        % (NB : on part de dz pour aller a e-dz afin d'eviter d'avoir plusieurs fois la meme abscisse)
        if idx_tor<nb_g_total_ligne % sauf pour le dernier
          z_e = [z(idx_tor)+b(idx_tor)+dz:dz:z(idx_tor+1)-dz];
          abs_z_idx_pol = [abs_z_idx_pol, z_e];
          fztot = [fztot, zeros(1,length(z_e))];
        end
%         %  Previous version 
%           z_norm = 0:1/100:1;
%           % mode TEM
%          fhz = ones(1,length(z_norm));
%           % modes TM
%          n = (1:Nme)';
%          fez = cos(pi*n*z_norm);
%  
%          fz = [fhz;fez];
%          if (Nme ~= 0)
%              fztot = (sum((coeff*ones(1,length(z_norm))).*fz));
%          else
%              fztot = ((coeff*ones(1,length(z_norm))).*fz);
%          end
%          
%          abscisse = (0:1/100:1.02)*b(idx_tor) + z(idx_tor);
%          fztot = [0,fztot,0];
%  
%          %abscisse = (0:1/100:1.02)*(b(idx_tor)./(b(idx_tor)+e)) + idx_tor;


        % AJOUT JH 10/2008 in order to save the results
        if (idx_tor == 1)
            Ez_idx_pol = fztot;
            Ez_avPhase_idx_pol = mean(angle(fztot))*ones(size(fztot));
        else
            Ez_idx_pol = [Ez_idx_pol, fztot];
            Ez_avPhase_idx_pol= [Ez_avPhase_idx_pol, mean(angle(fztot))*ones(size(fztot))];
        end

    end % idx_tor
   
    % variables pour sortie & sauvegarde
    abs_z(idx_pol,:) = abs_z_idx_pol;
    Efield(idx_pol,:)= Ez_idx_pol;
    Efield_average_phase(idx_pol,:)= Ez_avPhase_idx_pol;

end % idx_pol


% add the calculated results to the scenario
scenario.results.abs_z=abs_z;
scenario.results.dz=dz;
scenario.results.Efield=Efield;
scenario.results.Efield_average_phase=Efield_average_phase;


