% ALOHA initialisation script
% 
% Plusieurs operations indispensables au code ALOHA sont effectues ici :
% - Charge les fonctions necessaires au fonctionnement d'ALOHA dans le PATH de matlab.
% - Determine la constant ALOHA_PATH qui contient le chemin absolu vers le dossier racine d'ALOHA
% - Charge en memoire les constantes physiques fondamentales.
%  
% AUTHOR: JH
% LAST UPDATE
%  - 2008: creation

ALOHA_ROOT = aloha_utils_getRootPath;

% charge dans le PATH les fonctions utiles pour aloha
addpath(genpath([ALOHA_ROOT, '/libaloha']));
addpath(genpath([ALOHA_ROOT, '/architecture_antenne']));
addpath(genpath([ALOHA_ROOT, '/S_HFSS']));

% TODO : Determine la constant ALOHA_PATH qui contient le chemin absolu vers le dossier racine d'ALOHA

%  load physical constants and version
disp(aloha_message('Load usual physical constants.'));
aloha_constants;

disp(aloha_message('Define some ALOHA constants.'));
% parametres 1D
T_grill = 2;                % periodicite du grill / diminue les tps de calcul
D_guide_max = 100;          % decouplage : pr |i-j| > D_guide_max, Kij = 0 avec {i,j} numeros de deux guides => S  
                            % puis pr |i-j| > D_guide_max -5, Sij = 0 	/ diminue encore les tps de calcul
			    	  
erreur_rel = 1e-4;  	    % erreur sur les integrales
pertes = 1e-6;              % tangente delta pour eviter les singularites


% Cree le fichier .mat contenant Ez pour differentes valeurs de x [TODO]
bool_fichier_Ez_x_variable = false;


disp(aloha_message(' --------- computation start ---------'));
