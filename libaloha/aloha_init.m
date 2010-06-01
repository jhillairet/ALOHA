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
aloha_constants;
disp(aloha_message('Load usual physical constants.'))

disp(aloha_message(' --------- computation start ---------'));