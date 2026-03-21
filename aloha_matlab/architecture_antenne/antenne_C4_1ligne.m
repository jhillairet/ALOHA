% moitié de C4 sur 1 ligne poloidale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;			% nbre de lignes de guides poloidales 
nb_modules_tor = 8;		% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 1;    	% nbre de guides par module dans le sens poloidal 

nb_g_module_tor = 3;		% nbre de guides par module ds le sens toroidal 
pass_module_tor = [2];

nb_g_passifs_inter_modules = 1;	% nbre de guides passifs entre les modules
nb_g_passifs_bord = 1;		%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3;	% espacement entre les grills juxtaposés poloidalement 

a = 76e-3;			% hauteur des guides ds le sens poloidal  

antenne_standard = 0; % = 0 il faut décrire l'antenne sur une ligne : tableaux b, e et lcc 
                      % = 1 paramètres scalaires : b_g_actif, b_g_pass e et lcc

%b_g_actif = 14.65e-3;		% largeur des guides actifs 
%b_g_pass = 13e-3;		% largeur des guides passifs

%e = 3.825e-3;			% épaisseur des parois des guides dans le sens toroidal 

% largeur des guides 
mod = [14.65,13,14.65];  % un module interne : voie Active + Passive + Active (en mm)
mod_int = 13;               % Une voie Passive entre 2 modules internes (en mm)
mod_ext = 11;               % Une voie Passive à l'extérieur de l'antenne (en mm)
b =[mod_ext,mod,mod_int,mod,mod_int,mod,mod_int,mod,mod_int,mod,mod_int,mod,mod_int,mod,mod_int,mod,mod_ext]*1e-3;
% épaisseur entre les guides 
ep = 3.825;
ep_mod = [ep,ep*ones(1,2),ep];
e = [ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod]*1e-3;


lcc = 1/4;			% longueur du court-circuit (en lambda guidée);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangée sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';

chemin_aller = 'cd ../S_HFSS/matrices_HFSS_C4';
chemin_retour = 'cd ../../matlab';

% modules C4
%nom_fichiers = ['S_C4_24b';'S_C4_23b';'S_C4_22b';'S_C4_21b';'S_C4_14b';'S_C4_13b';'S_C4_12b';'S_C4_11b'];   
% nom_fichiers = ['S_C4_24h';'S_C4_23h';'S_C4_22h';'S_C4_21h';'S_C4_14h';'S_C4_13h';'S_C4_12h';'S_C4_11h'];
nom_fichiers = ['S_jorge'];

% parametres 1D

T_grill = 2;                % périodicité du grill / diminue les tps de calcul
D_guide_max = 100;          % découplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numéros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-4;  	    % erreur sur les intégrales
pertes = 1e-6;              % tangente delta pour éviter les singularités	  