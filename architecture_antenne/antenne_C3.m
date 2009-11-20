% half-C3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 3; % nbre de lignes de guides poloidales 
nb_modules_tor = 8; % nbre de modules sur une ligne poloidale 

nb_g_module_pol = 3; % nbre de guides par module dans le sens poloidal

nb_g_module_tor = 6; % nbre de guides par module ds le sens toroidal
pass_module_tor = [];

nb_g_passifs_inter_modules = 1; % nbre de guides passifs entre les modules 
nb_g_passifs_bord = 1; %nbre de guides passifs sur chaque bord

% distance entre la derniere ligne de guide d'une demi-antenne 
% et la premiere ligne de guide de l'autre demi-antenne : 8.2cm
% D'ou : distance entre les centres des demi-antennes :
% 82e-3 + 3*a + 2*espacement_g_pol = 0.3160 m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 12e-3; % espacement entre les grills juxtaposes poloidalement 

a = 70e-3; % hauteur des guides ds le sens poloidal 

antenne_standard = 0; % = 0 il faut d�crire l'antenne sur une ligne : tableaux b, e et lcc 
                      % = 1 param�tres scalaires : b_g_actif, b_g_pass e et lcc

b_g_actif = 8e-3;		% largeur des guides actifs 
b_g_pass = 6.5e-3;		% largeur des guides passifs

e = 2e-3; % epaisseur des parois des guides dans le sens toroidal
e_p = 3e-3; % epaisseur des parois des voies passives

lcc = 1/4; % longueur du court-circuit (en lambda guid�e);  

% MAJ JH 06/05/2009
% En realite la largeur des parois des voies passives n'est pas 
% celle des parois entre guides actifs : 3 mm au lieu de 2 mm. 
% Pour prendre en compte les dimensions realistes de l'antenne :
b = [b_g_pass, repmat([b_g_actif*ones(1,nb_g_module_tor), b_g_pass], 1, nb_modules_tor)];
e = repmat([e_p, e*ones(1,nb_g_module_tor-1), e_p], 1, nb_modules_tor);
lcc = 1/4*ones(1,nb_g_passifs_inter_modules*(nb_modules_tor-1) + 2*nb_g_passifs_bord);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rang�e sur une seule colonne)


chemin_retour = pwd;
chemin_aller = 'S_HFSS/matrices_HFSS_C3';
%  chemin_retour =  'cd ../../matlab';

% modules C3
% La modelisation HFSS de l'antenne débute au voisinage du CM.
% Pour prendre en compte le dephasage entre la mesure de phase et le debut de la modelisation HFSS
% il faut calculer le dephasage necessaire : 
%  phase incident + phase rallonge + phase s12 HFSS = phase en bout
% soit :
%  phase rallonge = dphi - phase s12 HFSS
% ou dphi est la correction (connue, cf. doc Annika) utilisee pour calculer la phase en bout
% a partir de la phase incidente.

% modules BAS
% -----------
% module    correction dphi     phase S12 HFSS    
% 24b (= 8B)  0                   166.4453
% 23b        -5.5                -160.9112
% 22b         8                  -135.2895
% 21b        14                  -122.0280
% 14b         8.5                -116.1493
% 13b        20                  -131.2212
% 12b        10                  -139.8052
% 11b (= 1B) 20                  -166.9180
%  nom_fichiers = ['S_C3_24b';'S_C3_23b';'S_C3_22b';'S_C3_21b';'S_C3_14b';'S_C3_13b';'S_C3_12b';'S_C3_11b'];   
%  phase_rallonge = (pi/180)*[0-166.44; -5.5+160.9; 8+135.3; 14+122; 8.5+116.1; 20+131.2; 10+139.8; 20+166.9]:


% modules HAUT
% -----------
% module    correction dphi     phase S12 HFSS    
% 24h (= 8H)    17               168.6118
% 23h           11              -154.7044
% 22h           1               -129.9639
% 21h           15              -118.7455
% 14h           16              -113.2875
% 13h           22              -117.8509
% 12h           17              -132.5663
% 11h (= 1B)    15              -164.3341        

nom_fichiers = ['S_C3_24h';'S_C3_23h';'S_C3_22h';'S_C3_21h';'S_C3_14h';'S_C3_13h';'S_C3_12h';'S_C3_11h'];
phase_rallonge = (pi/180)*[17-168.6; 11+154.7; 1+130; 15+118.7; 16+113.3; 22+117.9; 17+132.6; 15+164.3];



% parametres 1D

T_grill = 7;                % p�riodicit� du grill / diminue les tps de calcul
D_guide_max = 100;          % d�couplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} num�ros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
  
erreur_rel = 1e-6;  	    % erreur sur les int�grales
pertes = 1e-6;              % tangente delta pour �viter les singularit�s	

