% moitié de C2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 2;			% nbre de lignes de guides poloidales 
nb_modules_tor = 8;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 2;    	% nbre de guides par module dans le sens poloidal

nb_g_module_tor = 4;		% nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules 
nb_g_passifs_bord = 2;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3;	% espacement entre les grills juxtaposés poloidalement 

a = 72e-3;			% hauteur des guides ds le sens poloidal 

antenne_standard = 0; % = 0 il faut décrire l'antenne sur une ligne : tableaux b e et lcc  
                      % = 1 paramètres scalaires : b_g_actif, b_g_pass e et lcc

% largeur des guides 
b =9*ones(1,36)*1e-3;
% épaisseur entre les guides 
e = [2,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),6.4,2*ones(1,3),4.38,2*ones(1,3),2]*1e-3;
% longueur du court-circuit (en lambda guidée);  
lcc = [1/2,1/4,1/4,1/2];			


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangée sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_aller = 'cd ../S_HFSS/matrices_HFSS_JET';
chemin_retour = 'cd ../../matlab';

% modules C3
nom_fichiers = ['S_JET';'S_JET';'S_JET';'S_JET';'S_JET';'S_JET';'S_JET';'S_JET'];	



% parametres 1D

T_grill = 100;              % périodicité du grill / diminue les tps de calcul
D_guide_max = 100;          % découplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numéros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les intégrales
pertes = 1e-6;              % tangente delta pour éviter les singularités	

