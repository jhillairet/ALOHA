% moitie de C4 -half C4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

% nbre de lignes de guides poloidales - number of rows
nb_g_pol = 3;

% nbre de modules sur une ligne poloidale - number of modules in a row
nb_modules_tor = 8;

% nbre de guides par module dans le sens poloidal - number of waveguides per module in one column
nb_g_module_pol = 3;

% nbre de guides par module ds le sens toroidal - number of waveguides per module in one row
nb_g_module_tor = 3;

% nbre total de guides passifs par module -total number of passive waveguides per module in one row 
pass_module_tor = [2];  

% nbre de guides passifs entre les modules - number of passive weveguide between modules per row
nb_g_passifs_inter_modules = 1;

%nbre de guides passifs sur chaque bord - number of passive waveguides at each edge of the row 	
nb_g_passifs_bord = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

% espacement entre les grills juxtaposes poloidalement- spacing between rows
espacement_g_pol = 12e-3;	

% hauteur des guides ds le sens poloidal - height of waveguides 
a = 76e-3;

%  antenna_standard :
% = 0 il faut decrire l'antenne sur une ligne : 
%   tableaux b, e et lcc - array for b, e, lcc (width of WG, septum width, passive WHG length
% = 1 parametres scalaires : 
%   b_g_actif, b_g_pass e et lcc - scalar (1 value for all waveguides)
antenne_standard = 0;
 
% largeur des guides actifs 
b_g_actif = 14.65e-3;
% largeur des guides passifs
% MAJ JH 27/06/2008 : l'epaisseur des voies passives est en realite de 12 mm sur le C4 au lieu de 13 mm.
b_g_pass = 12e-3;
% Une voie Passive e l'exterieur de l'antenne 
% MAJ JH 27/06/2008 : l'epaisseur des voies passives ext est en realite de 10 mm sur le C4 au lieu de 11.
b_g_pass_ext = 10*1e-3; %11;               
    
% largeur des guides 
% un module interne : voie Active + Passive + Active 
mod = [b_g_actif, b_g_pass, b_g_actif];  

b =[b_g_pass_ext,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass,mod,b_g_pass_ext];


% Epaisseur des parois des guides dans le sens toroidal 
% MAJ JH 27/06/2008 : l'epaisseur de matière est plus importante de 0.5mm du a la taille relle des VP 
ep = 4.325e-3; %3.825e-3; % epaisseur entre les guides 

ep_mod = [ep,ep*ones(1,2),ep];
e = [ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod,ep_mod];


%lcc = 1/4;			% longueur du court-circuit (en lambda guidee);   passive WG length normalized to lambdag
lcc = 1/4*ones(1,17);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangee sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';
chemin_retour = pwd;
chemin_aller = 'S_HFSS/matrices_HFSS_C4';  % path to get the HFSS matrix 
%  chemin_retour = 'cd ../../matlab';                % back to the main directory


%  MAJ J.Hillairet le 14/12/07
%  Matrices S Antenne C4 calculees par J.Belo (HFSS) le 13/12/07
%  Commentaire : j'ai conserve les noms originels des fichiers transmis. 
%  La correspondance entre ces fichiers et les anciens est detaillee 
%  dans le document .doc situee dans le repertoire des .m
bool_newSmatrix = true;

% modules C4
if bool_newSmatrix
    nom_fichiers = strvcat('Module_bas_plaque_long7_ext2', ...
                    'Module_bas_plaque_long6_long7', ...     
                    'Module_bas_plaque_long5_long6', ...                
                    'Module_bas_plaque_long4_long5', ...                
                    'Module_bas_plaque_long3_long4', ...                           
                    'Module_bas_plaque_long2_long3', ...                
                    'Module_bas_plaque_long1_long2', ...                
                    'Module_bas_plaque_ext1_long1');
%nom_fichiers = ['Module_haut_plaque_long7_ext2'; ...
%               'Module_haut_plaque_long6_long7'; ...     
%               'Module_haut_plaque_long5_long6'; ...                
%               'Module_haut_plaque_long4_long5'; ...                
%               'Module_haut_plaque_long3_long4'; ...                           
%               'Module_haut_plaque_long2_long3'; ...                
%               'Module_haut_plaque_long1_long2'; ...                
%               'Module_haut_plaque_ext1_long1'];

else
	nom_fichiers = ['S_C4_24b';'S_C4_23b';'S_C4_22b';'S_C4_21b';'S_C4_14b';'S_C4_13b';'S_C4_12b';'S_C4_11b'];   
	% nom_fichiers = ['S_C4_24h';'S_C4_23h';'S_C4_22h';'S_C4_21h';'S_C4_14h';'S_C4_13h';'S_C4_12h';'S_C4_11h'];
end

phase_rallonge = zeros(8,1);
phase_rallonge = -pi/180.*[-56.5;-44.7;-36.5;0;0;-36.5;-44.7;-56.5];

% parametres 1D

T_grill = 2;                % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-4;  	    % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites	  
