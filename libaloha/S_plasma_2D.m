% Select the ALOHA2D binary file to use
switch version
    case 1
        chemin_binaire_fortran = '/code_2D/couplage_2D/ALOHA2D_F90';
        binary_name = 'ALOHA2D.OK';
    case 10
       chemin_binaire_fortran = '/code_2D/couplage_2D/1_couche_quad/bin';
       binary_name = 'couplage_plasma_2D';
    case 20
       chemin_binaire_fortran = '/code_2D/couplage_2D/1_couche_quad_fem/bin';
       binary_name = 'couplage_plasma_2D';
    case 60
       chemin_binaire_fortran = '/code_2D/couplage_2D/2_couches/bin';
       binary_name = 'couplage_plasma_2D';
    otherwise 
        error('Bad defined variable : ''version''');
end

% nombre d'onde
k0 = 2*pi*freq/celerite;

% nombre de modes code 2D
% ajout JH 01/2009
nbre_modes = Nmh+Nme;

Nmh = 1;	
Nme = nbre_modes-1;	

a_tot=a*ones(1,nbre_guides);
b_tot=kron(ones(1,nb_g_pol),b);
y_tot=kron(y,ones(1,nbre_guides/nb_g_pol));
z_tot=kron(ones(1,nb_g_pol),z);
    


%%%%%%%%%%%%%%%%%%%%%%%
% Changement de repertoire :
% - pour creer le fichier d'argument du code en fortran
% - lancer le code fortran
% - recuperer les resultats
% avant de revenir au repertoire actuel.
% 
% En cas d'erreur lors de l'execution du binaire, le programme affiche une alerte et continu.
chemin_retour = pwd;
% on se place dans le repertoire des binaires fortran 
% cd([aloha_utils_getRootPath, chemin_binaire_fortran]);
   
%  ecriture du fichier par_grill_2D.dat
%  
if (version == 1)
    aloha_save_ALOHA2D_inputFile('ALOHA2D.in', ...
        'wg_nb', nbre_guides, 'wg_modes_nb', nbre_modes, 'f', freq, ...
        'nz_min', nz_min, 'nz_max', nz_max, 'ny_min', ny_min, 'ny_max', ny_max, ...
        'b', b_tot, 'a', a_tot, 'y', y_tot, 'z', z_tot, ...
        'ne', ne0, 'dne', dne0, 'B0', B0);
    
elseif (version == 10 || version == 20 ) % modeles plasma 1 couche
    varCell = {nbre_guides nbre_modes B0 freq ne0 dne0 a_tot b_tot y_tot z_tot ny_min nz_min ny_max nz_max nbre_ny nbre_nz};
    
    aloha_saveFortranInputFile('par_grill_2D.dat', varCell);
%      save par_grill_2D.dat  -ascii 

elseif (version == 60) % modele plasma 2 couches
    
    save par_grill_2D.dat nbre_guides nbre_modes B0 freq ne0 [dne0 dne1] e_couche a_tot b_tot y_tot z_tot ny_min nz_min ny_max nz_max nbre_ny nbre_nz -ascii 

end



%
% execution du binaire fortran
% 
try
    binary_full_path = [aloha_utils_getRootPath,chemin_binaire_fortran,'/',binary_name];
    disp(aloha_message(['Lancement du binaire ', binary_full_path]));
    [status,result] = system(binary_full_path);

    if bool_debug
        disp(result);
    end
catch
    % cd(chemin_retour);
    error(aloha_message('?! Binary execution problem ?!'));
end

if (version == 1)
    % FORTRAN 90 ascii output
    % coupling matrix
    data=importdata('ALOHA2D.out.K.dat', ' ', 1);

    K_cpl = complex(data.data(:,1), data.data(:,2));
    
    % characteristic impedances
    data=importdata('ALOHA2D.out.Zc.dat');
    Zhe = complex(data(:,1), data(:,2));
    rac_Zhe = diag(sqrt(Zhe));
    
elseif (version == 10 || version == 20 || version == 60)

    %  % extraction des donn???es du fichier binaire K_cpl.dat
    %  fichier_data = fopen('K_cpl.dat','r');
    %  fread(fichier_data,1,'int64'); 
    %  K_cpl_tamp = fread(fichier_data,2*nbre_modes*nbre_guides*nbre_modes*nbre_guides,'double');   
    %  K_cpl = K_cpl_tamp(1:2:2*nbre_modes*nbre_guides*nbre_modes*nbre_guides) + ...
    %                  i*K_cpl_tamp(2:2:2*nbre_modes*nbre_guides*nbre_modes*nbre_guides);
    %  K_cpl=reshape(K_cpl,nbre_modes*nbre_guides,nbre_modes*nbre_guides);
    %  fread(fichier_data,1,'double');
    %  Zhe_tamp = fread(fichier_data,2*nbre_modes*nbre_guides,'double');
    %  Zhe = Zhe_tamp(1:2:2*nbre_modes*nbre_guides) + i*Zhe_tamp(2:2:2*nbre_modes*nbre_guides);
    %  fclose(fichier_data);

    %  Extraction des donnees des fichiers ascii
    [K_r, K_i]=textread('K_cpl2.dat', '(%f,%f)', 'headerlines', 2);
    [Zhe_r, Zhe_i]=textread('Zhe2.dat', '(%f,%f)', 'headerlines', 2);

    Zhe = complex(Zhe_r, Zhe_i);
    rac_Zhe = diag(sqrt(Zhe));

    K_cpl = complex(K_r, K_i);


end

% on s'attend ?? ce que la matrice soit une matrice carree
K_cpl=reshape(K_cpl, sqrt(length(K_cpl)), sqrt(length(K_cpl)));

%  K_cpl=reshape(K_cpl,nbre_modes*nbre_guides,nbre_modes*nbre_guides);


%  calcul de la matrice S_plasma
H = rac_Zhe*K_cpl*rac_Zhe;
S_plasma = inv(eye(length(Zhe))+H)*(eye(length(Zhe))-H); 
  
    
% cd(chemin_retour);

    
