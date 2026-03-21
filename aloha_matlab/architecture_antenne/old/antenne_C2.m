% moitié de C2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 2;			% nbre de lignes de guides poloidales 
nb_modules_tor = 8;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 2;    	% nbre de guides par module dans le sens poloidal

nb_g_module_tor = 4;		% nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules 
nb_g_passifs_bord = 1;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3;	% espacement entre les grills juxtaposés poloidalement 

a = 76e-3;	% 70e-3		% hauteur des guides ds le sens poloidal 

b_g_actif = 8.5e-3;		% largeur des guides actifs 
b_g_pass = 8.5e-3;	%6.5e-3	% largeur des guides passifs

e = 2e-3;			% épaisseur des parois des guides dans le sens toroidal

lcc = 1/4;			% longueur du court-circuit (en lambda guidée);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangée sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_aller = 'cd ../S_HFSS/matrices_HFSS_C2';
chemin_retour = 'cd ../../matlab';

% modules C3
nom_fichiers = ['S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2'];	



% parametres 1D

T_grill = 100;              % périodicité du grill / diminue les tps de calcul
D_guide_max = 100;          % découplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numéros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les intégrales
pertes = 1e-6;              % tangente delta pour éviter les singularités	

