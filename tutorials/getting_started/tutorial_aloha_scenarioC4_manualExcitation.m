function scenario = scenario_example 
% ALOHA scenario example
% VERSION : 1.5 (20/09/2011)
% 
% This file is a function which create an ALOHA scenario.
%
% All the parameters needed by ALOHA to make a run are stored inside
% this file, except the antenna description. All the parameters are detailled below.
% 
% This function returns a matlab stucture "scenario" which contains all the parameters
% used in ALOHA. This structure is used an input argument of the "aloha_scenario"
% functiion. Once its calculation done, ALOHA returns a similar structure "scenrio"
% which contains the results in addition (in the field scenario.results). 
%
% 
%  INPUT: none
%  
%  OUPUT: 
%   - scenario <structure>: ALOHA scenario
% 
% AUTHOR:JH

% ##################################################################
%  RESULTS SAVING 
% ##################################################################
% Automatic saving of the scenario (once calculated) ?
options.bool_save = false;
% matlab result filename
% (avoid spaces if possible)
options.scenario_filename = ''; 
% Commentaire sur le scenario 
options.comment = ['ALOHA tutorial - C4 antenna'];

% ==================================================================
%  ANTENNA DESCRIPTION
% ==================================================================

%% #####################
% antenna architecture:
%  'antenne_C2', 'antenne_C3', 'antenne_C4' ...
%  Ces valeurs correspondent aux noms des fichiers disponibles dans
%  le dossier 'achitecture_antenne'.
% 
antenna.architecture = 'antenna_C4_ITM';      

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
options.bool_homeMadeExcitation = true;

% User-made excitation.
% [only if bool_homeMadeExcitation = true & bool_mesure = false]
% 
% Warning : the length of the array must be correct 
% in respect to the number of module of your antenna !
% help aloha_antenna_excitation for some examples
antenna.a_ampl = sqrt(300)*ones(8,1);
antenna.a_phase = -150*(pi/180)*[0:7]';


% [bool_mesure=true only]
% Quel queusot doit-on utiliser pour fournir les informations d'excitation des antennes ?
% Which Tore Supra port should we use for feeding antenna input ?
% 
% 'Q6A' [C2 before 2009, C3 after] ou 'Q6B' [C3 before 2009, C4 after]
options.TSport = 'Q6B';   

% [bool_mesure=true only]
% numero du choc Tore Supra 
% Tore Supra pulse number
options.choc = []; 

% [bool_mesure=true only]
% Temps debut et fin de lecture [s]. 
% Les valeurs mesurees sont moyennees entre ces deux valeurs.
% Start and stop times for averaging experimental data [s]
options.tps_1 = []; % s
options.tps_2 = []; % s



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
% Examples:
plasma.ne0 = 5e17; % if bool_lignes_identiques = true
%  plasma.neO = [3;4;5]*1e17; % if bool_lignes_identiques = false

%% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la premiere couche de plasma en partant de l'antenne. [m]
% 
% Example : 
plasma.lambda_n(1) = 2e-3; % 2mm if bool_lignes_identiques = true
%  plasma.lambda_n(:,1) = [2;2;2]*1e-3; % if bool_lignes_identiques = false


%% #####################
% Epaisseur de la premiere couche de plasma [version=6 only]. [m]
% 
% Example : 2e-3 (2mm)
plasma.d_couche = 2e-3;   

%% #####################
% longueur de decroissance a l'embouchure de l'antenne
% pour la deuxieme couche de plasma (version=6). [m]
% 
% Exemple :
plasma.lambda_n(2) = 2e-2; % 2cm if bool_lignes_identiques = true
%  plasma.lambda_n(:,2) = [2;2;2]*1e-2; % if bool_lignes_identiques = false
 


%% #####################
% Epaisseur de la couche de vide entre l'antenne et 
% la premiere couche de plasma (version=5). [m]
plasma.d_vide = 0;


%% ####################
%  Intensite du champ magnetique devant l'antenne. [T]
% [OBSOLETE pour le moment]  
options.B0 = 2.95;        

%% #####################
% 2D or 3D normalization choice
% 
% In SWAN, the waveguides were described as parallel plate waveguides (2D). 
% In ALOHA-1D, the waveguides are parallel plate waveguides (2D)
% In ALOHA-2D, the waveguides are rectangular waveguides (3D). 
% 
% In order to be able to compare the Electric field amplitudes given 
% between ALOHA-1D and ALOHA-2D, one must renormalize the ALOHA-1D results
% and to take into account the height of the waveguide ("a" parameter).
% 
% Moreover, if one want to compare the scattering parameters given
% by ALOHA-1D/2D and SWAN, the scattering parameters reference impedance is
% not the same between the two codes. So a impedance renormalization if necessary.
% (in SWAN, the ref impedance is the vaccuum impedance, while it is the rectangular
% waveguide impedance in ALOHA)
% 
% type_swan_aloha = 0:
%  - SWAN reference impedance in S-parameters (parallel plate wg, TEM+TM mode)
%  - Electric field does not depends of the height of the waveguide ("a")
%    ALOHA-1D E-field should be of the same order than a 2D full wave simulation
%
% type_swan_aloha = 1: (default)
%  - ALOHA reference impedance in S-parameters (rectangular waveguide, TE+TM modes)
%  - Electric field depends of the waveguide height ("a"). 
%    ALOHA-1D or ALOHA-2D E-field should be of the same order.
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

% Compute and display the Efield in the waveguides at the mouth of the antenna
options.bool_compute_total_field = true;       
options.bool_display_total_field = false;

% Compute and display the Efield inside the few mm of the plasma
% The domain of calculation (a rectangle in front of the waveguides)
% should be defined in the parameters below.  
options.bool_compute_plasma_field = true;
options.bool_display_plasma_field = false;





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
plasma.dne(:,1) = plasma.ne0./plasma.lambda_n(:,1);   
% gradient de densite couche ne2 [version=6 only]
plasma.dne(:,2) = (1+plasma.d_couche./plasma.lambda_n(:,1)).*plasma.ne0./plasma.lambda_n(:,2);

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
options.nz_min = -20;   % depart en nz 
options.nz_max = 20;        % arrivee en nz
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

%%%%%%%%%%%%%%%%%%%%%%%% making void structures in order to avoid bugs later in filling the structure...
results.S_plasma = [];
antenna_lh.setup = [];

% #####################
%  Creation du scenario
%  
scenario = aloha_setfield([], plasma, antenna, options, results); 

% Load the antenna structure into the scenario
%
% NB : in previous version of ALOHA, this was made inside the aloha_scenario function. 
% However, in order to allow batch on antenna's dimension, the antenne dimensions have 
% been inserted into the scenario : this must be done before aloha_scenario, thus here.
scenario = aloha_setAntenna(scenario, antenna.architecture);
