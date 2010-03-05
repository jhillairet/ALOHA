function scenario = scenario_example 
% ALOHA scenario structure generation example
% 
% Ce fichier est une fonction qui genere un scenario pour ALOHA
% 
% Tous les parametres necessaires a une simulation d'ALOHA sont contenus
% dans ce fichier. Ces parametres sont detailles ci-dessous.
% 
% La fonction renvoie une structure "scenario", qui contient tous les parametres
% d'entree d'ALOHA. Ce fichier scenario est utilise comme parametre d'entree d'ALOHA.
% La fonction 'aloha_1D', une fois son calcul termine, va renvoyer cette structure, 
% ou les resultats de la simulation y sont enregistres.
% 
%  INPUT: none
%  
%  OUPUT: 
%   - scenario <structure>: ALOHA scenario
% 
% AUTHOR:JH
% LAST UPDATE:
%  - 10/2008: creation

% ##################################################################
%  RESULTS SAVING 
% ##################################################################
% Sauvegarde automatique du scenario resultat dans un fichier ?
options.bool_sauvegarde = false;
% matlab result filename
% (avoid spaces if possible)
options.scenario_filename = ''; 
% Commentaire sur le scenario 
options.comment = [''];

% ==================================================================
%  ANTENNA DESCRIPTION
% ==================================================================

%% #####################
% antenna architecture:
%  'antenne_C2', 'antenne_C3', 'antenne_C4' ...
%  Ces valeurs correspondent aux noms des fichiers disponibles dans
%  le dossier 'achitecture_antenne'.
% 
antenna.architecture = 'antenne_C2';      

%% #####################
% Antenna excitation
%  

% source frequency [Hz]
% 
antenna.freq = 3.7e9;

% Use Tore Supra data ? true/false
% 
% - true : on dispose de mesures reelles 
% - false : on utilise les valeurs predefinies (dans aloha_antenna_excitation)
% 
options.bool_mesure = false;    

% Use a predefenite excitation (true) 
% or a user excitation (false) [only if bool_mesure = false]
options.bool_homeMadeExcitation = false;

% User-made excitation.
% [only if bool_homeMadeExcitation = true & bool_mesure = false]
% 
% Warning : the length of the array must be correct 
% in respect to the number of module of your antenna !
% help aloha_antenna_excitation for some examples
antenna.a_ampl = [];
antenna.a_phase = [];


% [bool_mesure=true only]
% Quel queusot doit-on utiliser pour fournir les informations d'excitation des antennes ?
% Which Tore Supra port should we use for feeding antenna input ?
% 
% 'Q6A' [C2 before 2009, C3 after] ou 'Q6B' [C3 before 2009, C4 after]
options.TSport = 'Q6A';   

% [bool_mesure=true only]
% numero du choc Tore Supra 
% Tore Supra pulse number
options.choc = 39286; 

% [bool_mesure=true only]
% Temps debut et fin de lecture [s]. 
% Les valeurs mesurees sont moyennees entre ces deux valeurs.
% Start and stop times for averaging experimental data [s]
options.tps_1 = 4; % s
options.tps_2 = 7; % s



% ==================================================================
%   PLASMA DESCRIPTION
% ==================================================================
%% #####################
% Modele d'ALOHA : '1D' ou '2D'
options.version_code = '1D';

%% #####################
% Nombre de modes utilise pour decrire le champ dans les guides
% (Modele electromagnetique 1D)
% 
Nmh = 1;    % nbre de modes TE
Nme = 2;    % nbre de modes TM

%% #####################
% Choix de la version de modelisation du profil 
% de densite electronique devant l'antenne.
%  --------- pour la version 1D ------------
%  
% version 3: profil lineaire, methode des residus 
%            + domaine d'integration limite e t = 1/(max_nz+1)
% 
% version 5: profil lineaire, couche de vide 
%            + couche de plasma [Deneb only]
%  
% version 6: profil lineaire, deux couches de plasma.
%            ne0 et dne0 decrivent la couche ne1
%            
% version 7: plasma step [Deneb only]
version_plasma_1D = 3;

% --------- pour la version 2D ------------
% version_plasma_2D 1 : 1 couche. Integration numerique (libquadpack)
% version_plasma_2D 2 : 1 couche, Integration numerique et FEM
% version_plasma_2D 3 : 
% ...TODO...
% version_plasma_2D 6 : 2 couches, Integration numerique et FEM
version_plasma_2D = 2;

%% #####################
% Utiliser des lignes poloidales identiques ?
% 
% true  : densite Haut et Bas identiques 
%       => ne0 et dne0 scalaires. ex : ne0 = 0.2
% false : densite Haut et Bas differentes 
%       => ne0 et dne0 tableaux. ex : ne0 = [0.2,0.4,0.2]
% NB : les tableaux sont de la taille du nombre 
% de lignes poloidale d'un module (pour une demi-antenne)
% 
options.bool_lignes_identiques = true;     
 
%% #####################
% Densite electronique minimale devant l'antenne (x=0). [m^-3] 
% Rappels : 
%  - densite de coupure 3.7GHz sur TS : 1.71e17 m^-3
plasma.ne0 = 5e17;        

%% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la premiere couche de plasma en partant de l'antenne. [m]
% 
% Exemple : 2e-3 (2mm)
plasma.lambda_n(1) = 2e-3;   

%% #####################
% Epaisseur de la premiere couche de plasma [version=6 only]. [m]
% 
% Exemple : 2e-3 (2mm)
plasma.d_couche = 2e-3;   

%% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la deuxieme couche de plasma (version=6). [m]
% 
% Exemple : 2e-2 (2cm)
plasma.lambda_n(2) = 2e-2;   


%% #####################
% Epaisseur de la couche de vide entre l'antenne et 
% la premiere couche de plasma (version=5). [m]
plasma.d_vide = 0;


%% ####################
%  Intensite du champ magnetique devant l'antenne. [T]
% [OBSOLETE pour le moment]  
options.B0 = 2.95;        

%% #####################
% Definition de l'impedance caracteristique, de type ALOHA ou SWAN, 
% pour la renormalisation de la matrice de scattering S_grill-plasma.
% 
% = 0 impedance de mode ds le plasma type SWAN (mode TEM)
% = 1 impedance de mode ds le plasma type ALOHA (mode TE)
% NB :
% = 0 : valable uniquement pour les coefficient de reflexion ; 
% En champ, il y a un pb de normalisation a regler [TODO]
options.type_swan_aloha = 1;    


% ==================================================================
%   ALOHA options (true/false)
% ==================================================================
% display electronic density profile
options.bool_display_density_profile = false;

% N_parallel spectrum
options.bool_compute_spectrum = true; % compute n// spectrum
options.bool_display_spectrum = false; % display n// spectrum 

% Directivity
% NB : The spectrum is needed to compute the directivity !
options.bool_compute_directivity = true;
options.bool_display_directivity = false;   
    
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
options.definition_directivite = 1; 

% champ electrique dans l'embouchure de l'antenne
options.bool_compute_total_field = true;       
options.bool_display_total_field = false;

% champ dans le plasma
options.bool_compute_plasma_field = true;
options.bool_display_plasma_field = false;


options.bool_fichier_Ez_x_variable = false; % Cree le fichier .mat contenant Ez pour differentes valeurs de x




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
% path to the main root of the ALOHA code
if exist('aloha_utils_getRootPath') == 0
    disp('!!!! ERROR : the ALOHA functions are not registered into the MATLAB PATH !!!!');
    error('Please correct the previous error before restarting ALOHA. See INSTALL.m');
else
    options.aloha_path = aloha_utils_getRootPath;
end

% #####################
%  Debug mode option 
%  If true, display more informations.
% 
options.bool_debug = false;

% #####################
% Plasma edge properties
%   
% plasma electronic density model
if options.version_code == '1D'
    plasma.version = version_plasma_1D;
elseif options.version_code == '2D'
    plasma.version = version_plasma_2D;
else 
    error('Bad defined variables : ''version_code''');
end

% gradient de la premiere couche de plasma
plasma.dne(1) = plasma.ne0./plasma.lambda_n(1);   
% gradient de densite couche ne2 [version=6 only]
plasma.dne(2) = (1+plasma.d_couche./plasma.lambda_n(1)).*plasma.ne0./plasma.lambda_n(2);

% EM modes
options.modes = [Nmh, Nme]; 

% #####################
% Retrieve antenna excitation
% 
% It depends on the choice the user made: 
%  - use a predefinite excitation
%  - or use Tore Supra data
% 
if (options.bool_mesure)
    [antenna.a_ampl, antenna.a_phase] = ...
        aloha_antenna_excitation(options.choc, options.tps_1, options.tps_2, options.TSport);
else
    if not(options.bool_homeMadeExcitation) % predefinite excitation
        [antenna.a_ampl, antenna.a_phase] = ...
            aloha_antenna_excitation(antenna.architecture);
    end
end




%  %%%%%%%%%%%
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
% WARNING : depending of these parameter values, 
% aloha 2D may crash (seg fault). If it occurs, 
% increase the dnz or dny values or decrease nz_min/max, ny_min/max ! 
% TODO : aloha 2D spectrum/nz,ny memory bug
options.nz_min = -10;   % depart en nz 
options.nz_max = 10;        % arrivee en nz
options.dnz = 0.01;     % pas en nz

options.ny_min = -2;    % depart en ny 
options.ny_max = 2;     % arrivee en ny
options.dny = 0.1;        % pas en ny

% NB JH 01/2009 : compatibilite avec code 2D       
options.nbre_ny=(options.ny_max-options.ny_min)/options.dny;
options.nbre_nz=(options.nz_max-options.nz_min)/options.dnz;  

% ##################### TODO
% trace du chp dans le plasma 
% (pour les versions 3 et 6)
%  
options.z_coord_min = -0.015;       % coordonnees dans le plasma
options.z_coord_max = 0.075;
options.nbre_z_coord = 60;
options.x_coord_max = 0.05;
options.nbre_x_coord = 8;

options.pas_nz_fig_plasma = 0.1;    % description du chp spatial e partir de la description spectrale 
                        %     !! => periodisation du chp suivant z Tz=lambda/pas_nz
options.fig_Ez_ou_EzHy = 1;     % description chp Ez ou vecteur de poynthing S = Ez.Hy* :  = 1 Ez  / = 0 Ez.Hy*
options.lig_fig_plasma = 1;     % numero de la ligne

%%%%%%%%%%%%%%%%%%%%%%%%
results.S_plasma = [];

% #####################
%  Creation du scenario
%  
scenario = aloha_setfield([], plasma, antenna, options, results); 

