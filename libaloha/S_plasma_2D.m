% Determine le nom du binaire a utiliser
%  
% NB JH 01/2009 : Pour le code 2D, le nom du binaire est toujours le meme. 
% Par contre, c'est le repertoire du binaire qui change.
binary_name = 'couplage_plasma_2D';


switch version
    case 1
       chemin_binaire_fortran = '/code_2D/couplage_2D/1_couche_quad/bin';
    case 2
       chemin_binaire_fortran = '/code_2D/couplage_2D/1_couche_quad_fem/bin';
    case 6
       chemin_binaire_fortran = '/code_2D/couplage_2D/2_couches/bin';
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
cd([aloha_utils_getRootPath, chemin_binaire_fortran]);
   
%  ecriture du fichier par_grill_2D.dat
%  
if (version == 1|2|3 ) % modeles plasma 1 couche
    varCell = {nbre_guides nbre_modes B0 freq ne0 dne0 a_tot b_tot y_tot z_tot ny_min nz_min ny_max nz_max nbre_ny nbre_nz};
    
    aloha_saveFortranInputFile('par_grill_2D.dat', varCell);
%      save par_grill_2D.dat  -ascii 

elseif (version == 6) % modele plasma 2 couches
    
    save par_grill_2D.dat nbre_guides nbre_modes B0 freq ne0 [dne0 dne1] e_couche a_tot b_tot y_tot z_tot ny_min nz_min ny_max nz_max nbre_ny nbre_nz -ascii 

end



%
% execution du binaire fortran
% 
try
    disp(aloha_message(['Lancement du binaire ', binary_name]));
    [status,result] = system(['./',binary_name]);

    if bool_debug
        disp(result);
    end
catch
    cd(chemin_retour);
    error(aloha_message('?! Binary execution problem ?!'));
end


%  % extraction des donn�es du fichier binaire K_plasma.dat
%  fichier_data = fopen('K_plasma.dat','r');
%  fread(fichier_data,1,'int64'); 
%  K_plasma_tamp = fread(fichier_data,2*nbre_modes*nbre_guides*nbre_modes*nbre_guides,'double');   
%  K_plasma = K_plasma_tamp(1:2:2*nbre_modes*nbre_guides*nbre_modes*nbre_guides) + ...
%                  i*K_plasma_tamp(2:2:2*nbre_modes*nbre_guides*nbre_modes*nbre_guides);
%  K_plasma=reshape(K_plasma,nbre_modes*nbre_guides,nbre_modes*nbre_guides);
%  fread(fichier_data,1,'double');
%  Zhe_tamp = fread(fichier_data,2*nbre_modes*nbre_guides,'double');
%  Zhe = Zhe_tamp(1:2:2*nbre_modes*nbre_guides) + i*Zhe_tamp(2:2:2*nbre_modes*nbre_guides);
%  fclose(fichier_data);

%  Extraction des donnees des fichiers ascii
[K_r, K_i]=textread('K_plasma2.dat', '(%f,%f)', 'headerlines', 2);
[Zhe_r, Zhe_i]=textread('Zhe2.dat', '(%f,%f)', 'headerlines', 2);

Zhe = complex(Zhe_r, Zhe_i);
rac_Zhe = diag(sqrt(Zhe));

K_plasma = complex(K_r, K_i);

% on s'attend à ce que la matrice soit une matrice carree
K_plasma=reshape(K_plasma, sqrt(length(K_plasma)), sqrt(length(K_plasma)));

%  K_plasma=reshape(K_plasma,nbre_modes*nbre_guides,nbre_modes*nbre_guides);


%  calcul de la matrice S_plasma
H = rac_Zhe*K_plasma*rac_Zhe;
S_plasma = inv(eye(length(Zhe))+H)*(eye(length(Zhe))-H); 
  
    
cd(chemin_retour);

    
