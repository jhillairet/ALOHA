% ALOHA (main file for script mode)
%
% #####################################################################
% #                                                                   # 
% #                     -=  ALOHA =-                                  #
% #                                                                   #
% #              (Advanced LOwer Hybrid Antenna)                      #
% #                                                                   #
% #####################################################################
%  
%  Description du code : cf. README.m 
%  Installation du code : cf. INSTALL.m
%  
%  NB : au prealable, le repertoire contenant les librairies
%  d'ALOHA doit etre charge dans le PATH matlab. Pour cela, 
%  on doit effectuer dans matlab : 
%  >> addpath(genpath([chemin_du_code_ALOHA_V2/libaloha'])); 
%  Cette commande peut etre rajoutee dans le fichier startup.m de 
%  matlab pour effectuer cette operation a chaque demarrage de matlab
%  (cf help matlab)
%  
%  UTILISATION:
%  Tous les parametres necessaires pour lancer une simulation sont 
%  detailles dans ce fichier (mis a part les parametres intrinseques
%  des antennes qui sont specifies dans le dossier 'architecture_antenne').
%  L'execution de ce script Matlab doit permettre le lancement des calculs
%  par ALOHA sur un scenario.
%  
%  Une utilisation plus avancee d'ALOHA est possible, en concatenant les
%  scenarios. Cf. aide ALOHA_1D
%  
%  LAST CHANGES:
%   - 10/09/2008: ALOHA_1D et scenario : OK.

clear all; 

% #####################
%  Simulation comment.
%  This comment is saved into the scenario.
comment = ['test'];

% #####################
%  Save results in a matlab file (.mat) option
%  
bool_sauvegarde = false;
% matlab result filename
scenario_filename = 'sim_test.mat';

% ==================================================================
%  ANTENNA DESCRIPTION
% ==================================================================

%% #####################
% antenna architecture:
%  'antenne_C2', 'antenne_C3', 'antenne_C4' ...
%  Ces valeurs correspondent aux noms des fichiers disponibles dans
%  le dossier achitecture_antenne.
% 
architecture = 'antenna_C4_ITM';      
%architecture = 'antenne_C4';% Old fashioned antenna description. [OBSOLETE]

%% #####################
% Antenna excitation
%  

% source frequency [Hz]
%
% [OBSOLETE] : for compatibility with previous version
% NB : Following the new antenna definition (ITM)
% the source frequency should now be defined in the antenna description. 
freq = 3.7e9; 

% Use Tore Supra data ? true/false
% 
% - true : on dispose de mesures reelles 
% - false : on utilise les valeurs predefinies (dans aloha_antenna_excitation)
% 
bool_mesure = false;	

% Use a predefenite excitation (true) 
% or a user excitation (false) [only if bool_mesure = false]
bool_homeMadeExcitation = true;

% User-made excitation.
% [only if bool_homeMadeExcitation = true & bool_mesure = false]
% 
% Warning : the length of the array must be correct 
% in respect to the number of module of your antenna !
a_ampl = sqrt(1)*ones(8,1); 
% a_ampl([1,2,8]) = 0; % example of disfunctionnal klystrons : number 1,2 & 8
a_phase = (270*pi/180)*(0:7)';


% [bool_mesure=true only]
% Quel queusot doit-on utiliser pour fournir les informations d'excitation des antennes ?
% Which Tore Supra port should we use for feeding antenna input ?
% 
% 'Q6A' [C2 before 2009, C3 after] ou 'Q6B' [C3 before 2009, C4 after]
TSport = 'Q6A';   

%% numero du choc Tore Supra [bool_mesure = true]
choc = 39286; 

%% Temps debut et fin de lecture [s]. [bool_mesure = true]
% Les valeurs mesurees sont moyennees entre ces deux valeurs.
tps_1 = 4; 
tps_2 = 7;






% ==================================================================
%   PLASMA DESCRIPTION
% ==================================================================
%% #####################
% Modele d'ALOHA : '1D' ou '2D'
version_code = '1D';

%% #####################
% Nombre de modes utilise pour decrire le champ dans les guides
% (Modele electromagnetique 1D)
% 
Nmh = 1;    % nbre de modes TE
Nme = 2;    % nbre de modes TM 

%% #####################
% Choix de la version de modelisation du profil 
% de densite electronique devant l'antenne.
% --------- pour la version 1D ------------
%  
% version_plasma_1D 3: profil lineaire, methode des residus 
%            + domaine d'integration limite e t = 1/(max_nz+1)
% 
% version_plasma_1D 5: profil lineaire, couche de vide 
%            + couche de plasma [Deneb only]
%  
% version_plasma_1D 6: profil lineaire, deux couches de plasma.
%            ne0 et dne0 decrivent la couche ne1
%            
% version_plasma_1D 7: plasma step [Deneb only]
version_plasma_1D = 3;

% --------- pour la version 2D ------------
% version_plasma_2D 1 : 1 couche. Integration numerique (libquadpack)
% version_plasma_2D 2 : 1 couche, Integration numerique et FEM
% version_plasma_2D 3 : 
% ...TODO...
% version_plasma_2D 6 : 2 couches, Integration numerique et FEM
version_plasma_2D = 1;



% #####################
% Utiliser des lignes poloidales identiques ?
% 
% true  : densite Haut et Bas identiques 
%       => ne0 et dne0 scalaires. ex : ne0 = 0.2
% false : densite Haut et Bas differentes 
%       => ne0 et dne0 tableaux. ex : ne0 = [0.2,0.4,0.2]
% NB : les tableaux sont de la taille du nombre 
% de lignes poloidale d'un module (pour une demi-antenne)
% 
bool_lignes_identiques = true;      

% #####################
% Densite electronique minimale devant l'antenne (x=0). [m^-3] 
% Rappels : 
%  - densite de coupure 3.7GHz : 1.71e17
%  
ne0 = 5e17;

% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la premiere couche de plasma en partant de l'antenne. [m]
% 
% Exemple : 2e-3 (2mm)
lambda_n0 = 2e-3;

% #####################
% Epaisseur de la premiere couche de plasma (version=6 uniquement). [m]
% 
% Exemple : 2e-3 (2mm)
d_couche = 2e-3;    

% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la deuxieme couche de plasma (version=6 uniquement). [m]
% 
% Exemple : 2e-2 (2cm)
lambda_n1 = 2e-2;   


% #####################
% Epaisseur de la couche de vide entre l'antenne et 
% la premiere couche de plasma (version=5). [OBSOLETE?] [m]
d_vide = 0;     

% ####################
%  Intensite du champ magnetique devant l'antenne. [T]
B0_version7 = 2.95;

% #####################
% Definition de l'impedance caracteristique, de type ALOHA ou SWAN, 
% pour la renormalisation de la matrice de scattering S_grill-plasma.
% 
% = 0 impedance de mode ds le plasma type SWAN (mode TEM)
% = 1 impedance de mode ds le plasma type ALOHA (mode TE)
% NB :
% = 0 : valable uniquement pour les coefficient de reflexion ; 
% En champ, il y a un pb de normalisation a regler [TODO]
type_swan_aloha = 1;    




% ==================================================================
%   ALOHA options (true/false)
% ==================================================================

% display electronic density profile
bool_display_density_profile = false;

% N_parallel spectrum
bool_compute_spectrum = true; % compute n// spectrum
bool_display_spectrum = true; % display n// spectrum 

% Directivity
bool_compute_directivity = true;
bool_display_directivity = true;   
    
% directivity definition :
%  * definition_directivite=1 : 
%        D = (1/P)*int{nz=[+1,+inf]}(dP)
%   This is the most usual definition : ratio of the positive part of the spectrum 
%   over total power. This results may be considered in percents.
%  
%  * definition_directivite=2 : 
%       D = delta_cd = (1-R)*nz0^2/P*( int{nz=+1,+inf}(dP) - int{nz=-inf,-1}(dP))
%  This is the 'weighted directivity' as defined in [Litaudon & Moreau, Nucl Fusion 30 (1990), 471]
%  This directivity definition is based on the Fish current drive efficiency of the LH waves.
%  It is a quantitative measure, based on the CD efficiency which decrease as 1/n//^2, 
%  and of the quality of the launched spectrum with regard to an ideal Dirac spectrum, 
%  centered on n//0.
%  
%   * definition_directivite=3 : 
%       D = (1/P)*int{+1,+inf}(dP/nz^2);
definition_directivite = 1;

% (Parallel) Electric field in the mouth of the antenna
% champ electrique dans l'embouchure de l'antenne 
bool_compute_total_field = true;
bool_display_total_field = true;

% (Parallel) Electric field in some mm into the plasma
% [Plane-Wave propagation]
bool_compute_plasma_field = false;
bool_display_plasma_field = false;


%  ==================================================================
%  ==================================================================
%  ==================================================================
%  
%  DO NOT EDIT AFTER THIS LINE
%  (unless you know what you are doing ! :)
%  
%  ==================================================================
%  ==================================================================
%  ==================================================================

% #####################
%  Debug mode option 
%  If true, display more informations.
% 
bool_debug = false;

% #####################
%  set the ALOHA root path.
%  
if exist('aloha_utils_getRootPath') == 0
    disp('!!!! ERROR : the ALOHA functions are not registered into the MATLAB PATH !!!!');
    disp('You must add the "libaloha" directory to the MATLAB PATH');
    disp('In order to do so, please change directory (cd) to the aloha root directory, ');
    disp('from where you should see the "libaloha" directory.');
    disp('Then enter in the MATLAB prompt the command : "addpath(genpath([pwd, ''/libaloha'']));"');
    disp('This add the functions located into the directory "libaloha" (and sub-dir) to the MATLAB PATH');
    error('Please correct the previous error before restarting ALOHA.');
else
    aloha_path = aloha_utils_getRootPath;
end

% #####################
% Plasma edge properties
%   

% plasma electronic density model
if version_code == '1D'
    version = version_plasma_1D;
elseif version_code == '2D'
    version = version_plasma_2D;
else 
    error('Bad defined variables : ''version_code''');
end

% deduction du gradient de densite de la premiere couche
dne0 = ne0./lambda_n0;  
% Deduction du gradient de densite de la deuxieme couche
% La densite au bout de la couche d_couche est : 
%   ne0/lambda_n1*(1+d_couche/lambda_n0)
%   
dne1 = (1+d_couche./lambda_n0).*ne0./lambda_n1;


% #####################
% Retrieve antenna excitation
% 
% It depends on the choice the user made: 
%  - use a predefinite excitation
%  - or use Tore Supra data
% 
disp(aloha_message('Prise en compte de l''excitation antenne'));
if (bool_mesure) % measured excitation
    [a_ampl, a_phase] = ...
       aloha_antenna_excitation(choc, tps_1, tps_2, TSport);
else % theoritical excitation
    if not(bool_homeMadeExcitation) % predefinite excitation
        [a_ampl, a_phase] = ...
            aloha_antenna_excitation(architecture);
    end
end



%  % ##################### TODO
%  
%  % Affiche le profile du plasma de bord.
%  %
%  % TODO
%  if bool_display_density_profile
%      switch(version)
%          case 3
%              figure; aloha_plotDensityProfile(ne0, lambda_n0);
%          case 6
%              figure; aloha_plotDensityProfile(ne0, lambda_n0, d_couche, lambda_n1);
%      end    
%      
%      if bool_debug
%          disp('Appuyez sur une touche pour continuer');
%          pause
%      end
%  end

% #####################
% spectre computing parameters
%  
nz_min = -10;   % depart en nz 
nz_max = 10;    % arrivee en nz
dnz = 0.05;     % pas en nz

ny_min = -2;    % depart en ny 
ny_max = 2;     % arrivee en ny
dny = 0.5;      % pas en ny

% NB JH 01/2009 : compatibilite avec code 2D       
nbre_ny=(ny_max-ny_min)/dny;
nbre_nz=(nz_max-nz_min)/dnz;  


% ##################### 
% trace du chp dans le plasma 
% (pour les versions 3 et 6)
%  
z_coord_min = -0.005;   % coordonnees dans le plasma
z_coord_max = 0.06;
nbre_z_coord = 3000;
x_coord_max = 0.01;
nbre_x_coord = 5;

pas_nz_fig_plasma = 0.1;% description du chp spatial a partir de la description spectrale 
                        %     !! => periodisation du chp suivant z Tz=lambda/pas_nz
fig_Ez_ou_EzHy = 1;     % description chp Ez ou vecteur de poynthing S = Ez.Hy* :  = 1 Ez  / = 0 Ez.Hy*
lig_fig_plasma = 1;     % numero de la ligne



% #####################
%  Creation du scenario
%  

%  creation de la structure 'scenario' 
%  contenant tous les parametres de la simulation
aloha_scenario_create;


% Load the antenna structure into the scenario
%
% NB : in previous version of ALOHA, this was made inside the aloha_scenario function. 
% However, in order to allow batch on antenna's dimension, the antenne dimensions have 
% been inserted into the scenario : this must be done before aloha_scenario, thus here.
scenario = aloha_setAntenna(scenario, architecture);


% #####################
%  Main Program
%  
scenario = aloha_scenario(scenario);


% #####################
%  Results saving
%  
if (bool_sauvegarde)
    disp(aloha_message(['Sauvegarde des results dans le fichier ', scenario_filename]));
    status = aloha_scenario_appendScenarioToFile(scenario, scenario_filename);
end


disp(aloha_message('--------- fin du calcul.---------'));
