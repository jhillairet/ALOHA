% JH 06/10/2009
% COMPASS-D LHCD SYSTEM
% This antenna consists in 8 independant waveguides
% fed by 2*300kW 1.3 GHz source
% Port size is 170 x 140 mm
% Antenne size is 170 x 120 mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;			% nbre de lignes de guides poloidales 
nb_modules_tor = 8;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 1;    	% nbre de guides par module dans le sens poloidal

nb_g_module_tor = 1;		% nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules 
nb_g_passifs_bord = 0;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 0;	% espacement entre les grills juxtapos�s poloidalement 

a = 16.50e-2;			% hauteur des guides ds le sens poloidal 

antenne_standard = 1; % = 0 il faut d�crire l'antenne sur une ligne : tableaux b, e  et lcc
                      % = 1 param�tres scalaires : b_g_actif, b_g_pass, e et lcc

b_g_actif = 1.48e-2;		% largeur des guides actifs 
b_g_pass = 0;		% largeur des guides passifs

e = 0.2e-2;			% �paisseur des parois des guides dans le sens toroidal

lcc = 1/4;			% longueur du court-circuit (en lambda guid�e);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang�e sur une seule colonne)

chemin_retour =  pwd;
chemin_aller = 'S_HFSS/matrices_HFSS_Simple_Antenna';

nom_fichiers =  repmat(['S_1active_wg1'], 8, 1);

% parametres 1D

T_grill = 1;                % p�riodicit� du grill / diminue les tps de calcul
D_guide_max = 100;          % d�couplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} num�ros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les int�grales
pertes = 1e-6;              % tangente delta pour �viter les singularit�s	

