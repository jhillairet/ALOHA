%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;			% nbre de lignes de guides poloidales 

% mettre autant de modules que de guides
nb_modules_tor = 8;		% nbre de modules sur une ligne poloidale 

% mettre 1
nb_g_module_pol = 1;    	% nbre de guides par module dans le sens poloidal

% mettre 1
nb_g_module_tor = 1;		% nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules 
nb_g_passifs_bord = 0;		% nbre de guides passifs sur chaque bord 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3;	% espacement entre les grills juxtapos�s poloidalement 

a = 60e-3;			% hauteur des guides ds le sens poloidal  

antenne_standard = 1; % = 0 il faut d�crire l'antenne sur une ligne : tableaux b et e  
                      % = 1 param�tres scalaires : b_g_actif, b_g_pass et e


b_g_actif = 5.5e-3;		% largeur des guides actifs 
b_g_pass = 6.5e-3;		% largeur des guides passifs 

e = 1.5e-3;			% �paisseur des parois des guides dans le sens toroidal 

lcc = 1/4;			% longueur du court-circuit (en lambda guid�e); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang�e sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_retour = pwd;
chemin_aller = 'cd ../S_HFSS/matrices_HFSS_elem';



% mettre autant de fichiers que de guides

nom_fichiers = ['S_elem'];
for ind_mod_fich = 2:nb_modules_tor
    nom_fichiers = [nom_fichiers;'S_elem'];
end




% parametres 1D

T_grill = 1;                % p�riodicit� du grill / diminue les tps de calcul
D_guide_max = 100;          % d�couplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} num�ros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les int�grales
pertes = 1e-6;              % tangente delta pour �viter les singularit�s	  
