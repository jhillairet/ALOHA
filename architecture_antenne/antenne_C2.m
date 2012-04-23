% moiti� de C2

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

espacement_g_pol = 12e-3;	% espacement entre les grills juxtapos�s poloidalement 

a = 76e-3;			% hauteur des guides ds le sens poloidal 

antenne_standard = 0; % = 0 il faut d�crire l'antenne sur une ligne : tableaux b e et lcc 
                      % = 1 param�tres scalaires : b_g_actif, b_g_pass e et lcc

b_g_actif = 8.5e-3;       % largeur des guides actifs 
b_g_pass = 6.5e-3;      % largeur des guides passifs

% largeur des guides 
b =[b_g_pass, b_g_actif*ones(1,32), b_g_pass];
% �paisseur entre les guides 
e = [2,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),4.38,2*ones(1,3),2]*1e-3;
% longueur du court-circuit (en lambda guid�e);  
lcc = (1/4)*ones(1,2);			


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang�e sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';
chemin_retour = pwd;
%  chemin_aller = 'S_HFSS/matrices_HFSS_C2';

% modules C2
%  names of the modules (as viewed from the plasma)
%  
%  upper : 1B 2B 1B 2B 1B 2B 1B 2B
%  lower : 1H 2H 1H 2H 1H 2H 1H 2H

% Matrice S theorique (D.Voyer)
%  nom_fichiers = ['S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2';'S_C2'];	

% Ajout J.Hillairet 22/1/2008 
% Matrice S theorique du module "0" provenant du document TS9789 du 11/4/89 de Ph.Bibet et J.Achard
%  nom_fichiers = ['S_C2_0';'S_C2_0';'S_C2_0';'S_C2_0';'S_C2_0';'S_C2_0';'S_C2_0';'S_C2_0'];	

% Matrice S experimentale du module "0" provenant du document TS9789 du 11/4/89 de Ph.Bibet et J.Achard
%  nom_fichiers = ['S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp';'S_C2_0_exp'];	

%  %  Matrices S provenant du document TS9789 du 11/4/89 de Ph.Bibet et J.Achard
%  Lower modules (after antenna rotation)
%  nom_fichiers = ['S_C2_1B';'S_C2_2B';'S_C2_1B';'S_C2_2B';'S_C2_1B';'S_C2_2B';'S_C2_1B';'S_C2_2B'];   
%  %  Upper modules (after antenna rotation)
%  nom_fichiers = ['S_C2_1H';'S_C2_2H';'S_C2_1H';'S_C2_2H';'S_C2_1H';'S_C2_2H';'S_C2_1H';'S_C2_2H'];	

%  % Simulation des matrices S par J.Belo en Janvier 2008 
%  nom_fichiers = [repmat('JBelo_01_2008_',8,1), ...
%                  ['C2_Module_2H';'C2_Module_1H';'C2_Module_2H';'C2_Module_1H';'C2_Module_2H';'C2_Module_1H';'C2_Module_2H';'C2_Module_1H']];   


%% Simulation des matrices S par J.Belo en Juillet 2008
%  chemin_aller = 'S_HFSS/matrices_HFSS_C2/JBelo_07_2008';
%                   
%  module #1, le plus a droite vu du plasma : C2_Module_2B_R_8_H
%  ...
%  module #8, le plus a gauche vu du plasma : C2_Module_1B_R_1_H
%  
% modules HAUT de C2 (apres retournement)
%  nom_fichiers = [repmat('JBelo_07_2008_C2_Module_',8,1), ...
%                  ['2B_R_8_H';'1B_R_7_H';'2B_R_6_H';'1B_R_5_H';'2B_R_4_H';'1B_R_3_H';'2B_R_2_H';'1B_R_1_H']]; 
%  
%  nom_fichiers = [repmat('JBelo_07_2008_C2_Module_',8,1), ...
%                  ['1B_R_1_H';'2B_R_2_H';'1B_R_3_H';'2B_R_4_H';'1B_R_5_H';'2B_R_6_H';'1B_R_7_H';'2B_R_8_H']]; 

%  % modules BAS de C2 (apres retournement) 
%  nom_fichiers = [repmat('JBelo_07_2008_C2_Module_',8,1), ...
%                  ['2H_R_8_B';'1H_R_7_B';'2H_R_6_B';'1H_R_5_B';'2H_R_4_B';'1H_R_3_B';'2H_R_2_B';'1H_R_1_B']]; 
%  
%  nom_fichiers = [repmat('JBelo_07_2008_C2_Module_',8,1), ...
%                  ['1H_R_1_B';'2H_R_2_B';'1H_R_3_B';'2H_R_4_B';'1H_R_5_B';'2H_R_6_B';'1H_R_7_B';'2H_R_8_B']]; 

%%  dephasage mesure entre les fenetres HF et l'entree de la jonction hybride
%  
%  JH 4/11/2008 : pour determiner le dephasage necessaire entre les fenetres HF 
%  et l'entree de la jonction hybride, j'ai utilise les dernieres matrices S realisees
%  par J.Belo (07/2008). 
%  
%  Connaissant la correction de dephasage (delta_phi) utilisee sur TS (relative au module 8B)
%  (A.Ekedahl : rhyb/hybeta/Coef_Hyb20081006.dat) 
%  (rappel : phase_en_bout = phase_inc + delta_phi)
%  
%  NB: les mesures sont relatives au module 8B
%  
% C2 : numerotation TS des modules situes physiquement en HAUT du coupleur :
%  1H    5.00 % module #1, le plus a droite vu du plasma
%  2H  -49.00
%  3H   17.00
%  4H  -64.00
%  5H    6.00
%  6H  -49.00
%  7H    7.00
%  8H  -63.00 % module #8, le plus a gauche vu du plasma
%
% C2 : numerotation TS des modules situes physiquement en BAS du coupleur :
%  1B  -57.00 % module #1, le plus a droite vu du plasma
%  2B   14.00
%  3B  -56.00
%  4B    0.00
%  5B  -61.00
%  6B   -3.00
%  7B  -65.00
%  8B    0.00 % module #8, le plus a gauche vu du plasma
%  
%  Connaissant la phase de la première voie de la matrice S J.Belo (07/2008): angle(S(1,2)) 
%  ATTENTION : la nomenclature des modules ne correspond pas a la nomenclature TS !!
%  
%  C2 : modules situes physiquement en HAUT du coupleur :
%  2B_R_8_H  139.8274 % module #1, le plus a droite vu du plasma
%  1B_R_7_H  99.0504
%  2B_R_6_H  145.4388
%  1B_R_5_H  104.7924
%  2B_R_4_H  145.4938
%  1B_R_3_H  100.9410
%  2B_R_2_H  141.6767
%  1B_R_1_H  89.7361 % module #8, le plus a gauche vu du plasma
%  
%  C2 ; Modules situes physiquement en BAS du coupleur :
%  2H_R_8_B  111.3714 % module #1, le plus a droite vu du plasma
%  1H_R_7_B  153.3568
%  2H_R_6_B  117.0598
%  1H_R_5_B  160.7192
%  2H_R_4_B  117.0220
%  1H_R_3_B  153.2958
%  2H_R_2_B  111.3019
%  1H_R_1_B  147.6672 % module #8, le plus a gauche vu du plasma
%    
%  j'en deduis les corrections necessaires : phase_rallongue = delta_phi - angle(S12)

%  modules HAUT : ['2B_R_8_H';'1B_R_7_H';'2B_R_6_H';'1B_R_5_H';'2B_R_4_H';'1B_R_3_H';'2B_R_2_H';'1B_R_1_H']
%  phase_rallonge = pi/180*[5-139.8274; -49-99.0504; 17-145.4388; -64-104.7924; 6-145.4938; -49-100.9410; 7-141.6767; -63-89.7361];

% modules BAS : ['2H_R_8_B';'1H_R_7_B';'2H_R_6_B';'1H_R_5_B';'2H_R_4_B';'1H_R_3_B';'2H_R_2_B';'1H_R_1_B']
%phase_rallonge = pi/180*[-57-111.3714; 14-153.3568; -56-117.0598; 0-160.7192; -61-117.0220; -3-153.2958; -65-111.3019; 0-147.6672];


%% Simulation des matrices S par J.Hillairet en Mars 2012
chemin_aller = 'S_HFSS/matrices_HFSS_C2/JHillairet_03_2012';
% upper modules
nom_fichiers = [repmat('C2_Module_',8,1), ...
                ['1B';'2B';'1B';'2B';'1B';'2B';'1B';'2B']]; 
%  % lower modules
%  nom_fichiers = [repmat('C2_Module_',8,1), ...
%                  ['1H';'2H';'1H';'2H';'1H';'2H';'1H';'2H']]; 

phase_rallonge = [0,0,0,0,0,0,0,0];

% parametres 1D

T_grill = 100;              % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} num�ros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
  
erreur_rel = 1e-4;  	    % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites	

