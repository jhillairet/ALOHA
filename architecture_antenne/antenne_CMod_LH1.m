% JH 02/2009
% Antenne LH du Tokamak ALCATOR C-Mod
% Grill classique composé de 22 (24 - 2 HS) guides d'ondes phasés
% (sur 4 lignes poloidale : ici, on n'en considere qu'une seule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;% nbre de lignes de guides poloidales 
nb_modules_tor = 22;% nbre de modules sur une ligne poloidale 

nb_g_module_pol = 1;    % nbre de guides par module dans le sens poloidal

nb_g_module_tor = 1;% nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 0;% nbre de guides passifs entre les modules 
nb_g_passifs_bord = 0;%nbre de guides passifs sur chaque bord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3;% espacement entre les grills juxtaposes poloidalement 

a = 60e-3;% hauteur des guides ds le sens poloidal 

antenne_standard = 1; % = 0 il faut decrire l'antenne sur une ligne : tableaux b, e  et lcc
                      % = 1 parametres scalaires : b_g_actif, b_g_pass, e et lcc

b_g_actif = 5.5e-3;% largeur des guides actifs 
b_g_pass = 5.5e-3;% largeur des guides passifs

e = 1.5e-3;% epaisseur des parois des guides dans le sens toroidal

lcc = 1/4;% longueur du court-circuit (en lambda guidee);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangee sur une seule colonne)

chemin_retour =  pwd;
chemin_aller = 'S_HFSS/matrices_HFSS_Simple_Antenna';

nom_fichiers =  repmat(['S_1active_wg1'], 22, 1);

% parametres 1D

T_grill = 1;                % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 / diminue encore les tps de calcul
      
erreur_rel = 1e-6;      % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites

