% Antenne FAM for Iter - 5GHz 
% Author : J.Belo
% Date : 12/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;			% nbre de lignes de guides poloidales 
nb_modules_tor = 1;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 1;    	% nbre de guides par module dans le sens poloidal 

nb_g_module_tor = 8;		% nbre de guides par module ds le sens toroidal 
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules
nb_g_passifs_bord = 0;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 1;	% espacement entre les grills juxtapose poloidalement 

a = 50e-3;			% hauteur des guides ds le sens poloidal  

antenne_standard = 1; % = 0 il faut derire l'antenne sur une ligne : tableaux b, e et lcc 
                      % = 1 parameres scalaires : b_g_actif, b_g_pass e et lcc

b_g_actif = 6.99e-3;		% largeur des guides actifs 
b_g_pass = 1;		% largeur des guides passifs

e = 3e-3;			% eaisseur des parois des guides dans le sens toroidal 

lcc = 1/4;			% longueur du court-circuit (en lambda guide);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est range sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_aller = 'cd ../S_HFSS/matrices_HFSS_ITER';
chemin_retour = 'cd ../../matlab';

% modules
nom_fichiers = ['FAM_8ActiveWG_120_short_module'];   


% parametres 1D

T_grill = 1;                % peiodicitedu grill / diminue les tps de calcul
D_guide_max = 100;          % deouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeos de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les interales
pertes = 1e-6;              % tangente delta pour eiter les singularite	  