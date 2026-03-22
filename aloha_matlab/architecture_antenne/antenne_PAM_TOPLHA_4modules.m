% moitie de C4 -half C4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% architecture

nb_g_pol = 1;       % nbre de lignes de guides poloidales - number of rows
nb_modules_tor = 4;	% nbre de modules sur une ligne poloidale - number of modules in a row

nb_g_module_pol = 1; % nbre de guides par module dans le sens poloidal - number of waveguides per module in one column

nb_g_module_tor = 7; % nbre de guides par module ds le sens toroidal - number of waveguides per module in one row
pass_module_tor = [2 4 6]; % INDEX des voies passives d'un module 

nb_g_passifs_inter_modules = 0;	% nbre de guides passifs entre les modules - number of passive weveguide between modules per row
nb_g_passifs_bord = 1;	%nbre de guides passifs sur chaque bord - number of passive waveguides at each edge of the row 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions

espacement_g_pol = 0;	% espacement entre les grills juxtaposes poloidalement- spacing between rows

a = 76e-3;	% hauteur des guides ds le sens poloidal - height of waveguides 

antenne_standard = 0; % = 0 il faut decrire l'antenne sur une ligne : tableaux b, e et lcc - array for b, e, lcc (width of WG, septum width, passive WHG length
                      % = 1 parametres scalaires : b_g_actif, b_g_pass e et lcc - scalar (1 value for all waveguides)

%b_g_actif = 14.65e-3;		% largeur des guides actifs 
%b_g_pass = 13e-3;		% largeur des guides passifs

%e = 3.825e-3;			% epaisseur des parois des guides dans le sens toroidal 

% largeur des guides 
mod = [14.65, 13, 14.65, 13, 14.65, 13, 14.65];  % un module interne : voie Active + Passive + Active (en mm)
          
VP_int = 13; % Une voie Passive interne entre 2 modules (en mm)
VP_ext = 13; % Une voie Passive e l'exterieur de l'antenne (en mm)
b = [VP_ext, repmat(mod,1,4), VP_ext]*1e-3;

ep = 3.825e-3; % epaisseur entre les guides 

e = ep*ones(1,length(b)-1);


%lcc = 1/4;	 % longueur du court-circuit (en lambda guidee);   passive WG length normalized to lambdag
lcc = 1/4*ones(1, length(pass_module_tor)*nb_modules_tor+2*nb_g_passifs_bord);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrice S des modules ds des fichiers .m (NB : la matrice est rangee sur une seule colonne)

% chemin_aller = '';
% chemin_retour = '';
chemin_retour =  pwd;
chemin_aller = 'S_HFSS/matrices_HFSS_PAM_TOPLHA';  % path to get the HFSS matrix 
%  chemin_retour = 'cd ../../matlab';                % back to the main directory

% modules       
nom_fichiers = ['PAM_TOPLHA';'PAM_TOPLHA';'PAM_TOPLHA';'PAM_TOPLHA'];   

%  dephasage mesure entre les fenetres HF et l'entree du CM pour les 8 modules
phase_rallonge = zeros(1, nb_modules_tor);






% parametres 1D

T_grill = 2;                % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-4;  	    % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites	  
