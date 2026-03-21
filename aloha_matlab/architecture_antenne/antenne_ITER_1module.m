% QUART de l'antenne ITER (soit 5MW) - 5GHz - DDD2001

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;			% nbre de lignes de guides poloidales 
nb_modules_tor = 1;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 1;    	% nbre de guides par module dans le sens poloidal 

nb_g_module_tor = 15;		% nbre de guides par module ds le sens toroidal 
pass_module_tor = [2,4,6,8,10,12,14];

nb_g_passifs_inter_modules = 1;	% nbre de guides passifs entre les modules
nb_g_passifs_bord = 1;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 10e-3;	% espacement entre les grills juxtapos??? poloidalement 

a = 58e-3;			% hauteur des guides ds le sens poloidal  

antenne_standard = 1; % = 0 il faut d???rire l'antenne sur une ligne : tableaux b, e et lcc 
                      % = 1 param???res scalaires : b_g_actif, b_g_pass e et lcc

b_g_actif = 9.25e-3;		% largeur des guides actifs 
b_g_pass = 7.25e-3;		% largeur des guides passifs

e = 3e-3;			% ???aisseur des parois des guides dans le sens toroidal 

lcc = 1/4;			% longueur du court-circuit (en lambda guid???);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang??? sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_aller = 'cd ../S_HFSS/matrices_HFSS_ITER';
chemin_retour = 'cd ../../matlab';

% modules
% nom_fichiers = ['LH4Iter_module''];   

 nom_fichiers = ['LH4Iter_module_longVersion_58mm'];   
%  nom_fichiers = ['LH4Iter_module_shortVersion'];   
%  nom_fichiers = ['LH4Iter_module_longVersion_comsol'];   


% parametres 1D

T_grill = 1;                % p???iodicit???du grill / diminue les tps de calcul
D_guide_max = 100;          % d???ouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} num???os de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les int???rales
pertes = 1e-6;              % tangente delta pour ???iter les singularit???	  