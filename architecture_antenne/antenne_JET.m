% 2 rows of the JET LHCD Antenna (8x4=32 active waveguides per row)

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

espacement_g_pol = 12e-3;	% espacement entre les grills juxtaposes poloidalement 

a = 72e-3;			% hauteur des guides ds le sens poloidal 

antenne_standard = 0; % = 0 il faut decrire l'antenne sur une ligne : tableaux b e et lcc  
                      % = 1 parametres scalaires : b_g_actif, b_g_pass e et lcc

% largeur des guides 
b_g_pass_ext = 9e-3;
b =[b_g_pass_ext,b_g_pass_ext, [repmat([9 9 9 9],1,nb_modules_tor)]*1e-3, b_g_pass_ext,b_g_pass_ext];

% epaisseur entre les guides 
e =[2, 2, repmat([2,2,2,2],1,nb_modules_tor),2]*1e-3;
%  e = 2e-3
% longueur du court-circuit (en lambda guidee);  
lcc = (1/4)*ones(1,2*nb_g_passifs_bord+nb_g_passifs_inter_modules*(nb_modules_tor-1));			


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang�e sur une seule colonne)

chemin_aller = 'S_HFSS/matrices_HFSS_JET';
chemin_retour = pwd;

% modules C3
nom_fichiers = repmat('S_JET_damien',nb_modules_tor, 1);


% parametres 1D

T_grill = 100;              % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-6;  	    % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites	

