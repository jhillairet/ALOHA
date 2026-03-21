function [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor,varargout]=aloha_utils_getAntennaCoordinates(architecture,varargin)
% Renvoie les coordinees qui definissent les dimensions de l'antenne ; largeur, hauteur et position des guides
% dans l'espace.
% 
% NB : 
%  - on considere que le z=0 est situe à droite vue du plasma
%  - on considere que le y=0 est situe sur la plus basse ligne de guides.
%  
%  [b,h,z,y,nb_g_total_ligne,nbre_guides[,scenario]]=aloha_utils_getAntennaCoordinates(architecture[,scenario])
%  
%  INPUT : 
%   - architecture <string> : nom de l'architecture. 
%   - [optionnal] : scenario <ALOHA scenario>
%   NB : Ce nom doit correspondre a un nom de fichier existant dans le repertoire des architectures antennes
%  
%  OUPUT :
%   - b [vecteur (1,nb_g_total_ligne)] : largeurs des guides sur une ligne poloidale
%   - h : hauteur des guides sur une ligne poloidale (constante pour tous les guides d'une ligne)
%   - z [vecteur (1,nb_g_total_ligne)]: abscisse des guides 
%   - y [vecteur (1,nb_g_pol)]: ordonnée des guides 
%   - nb_g_total_ligne : nombre total de guides par ligne poloidale
%   - nbre_guides : nombre total de guides par antenne
%   - act_module_tor [vecteur (1,nb_g_actif=nb_g_module_tor-length(pass_module_tor)] : index des guides actifs dans un module
%   - [optionnal] : scenario <ALOHA scenario>
% 
% AUTHOR: Damien Voyer
% LAST UPDATE:
%  - 10/2008,JH : mise sous forme d'une fonction au lieu d'un script. 
%   

%% test if the architecture input is correct
if exist(architecture)
    eval(architecture);    
else
    error(['Architecture : ', architecture,' doesn''t exist!']);
end


%% treat optionnal input/output arguments if any
% JH 12/10/2010
% The main purpose of this code is to save the geometrical properties of 
% the antenna into the scenario (ie. a,b,e,z,y), in order to facilitate parametric studies on
% antenna dimensions.
% 
% If there is an optionnal input argument
if not(isempty(varargin)) 
    % this optionnal input argument should be a scenario
    % thus we test if this input is a matlab structure
    if isstruct(varargin{1})
        % OK, it's an ALOHA structure.
        % Import the scenario
        scenario=varargin{1};

        % Then we check if the antenna geometrical properties are
        % stored into the structure. 
        % If this is the case, these values are used 
        % instead of the ones defined in the architecture 
        % (Yes, this code is a bit crappy, but hum...)
        if isfield(scenario.antenna, {'b', 'a', 'e'})
            a=scenario.antenna.a;
            b=scenario.antenna.b;
            e=scenario.antenna.e;
        end
    else
        % optinnal input is not a matlab structure
        error('Bad optional input argument: should be an ALOHA scenario!')
    end
end



% nombre total de guides par ligne
nb_g_total_ligne = nb_g_module_tor*nb_modules_tor + 2*nb_g_passifs_bord + nb_g_passifs_inter_modules*(nb_modules_tor - 1);
%  nombre total de guides 
nbre_guides=nb_g_total_ligne*nb_g_pol;

b_module = zeros(1,nb_g_module_tor);    % largeur des guides dans le sens toroidal

act_module_tor = 1:nb_g_module_tor;
act_module_tor(:,pass_module_tor) = [];

% create the z and y arrays which contains the width and height of 
% all the waveguides of the antenna
% 
% TODO : remove the option antenne_standard and force to define antenna
% by arrays of a,b,e
if (antenne_standard == 1) 

   b_module(act_module_tor) = b_g_actif;
   b_module(pass_module_tor) = b_g_pass;

   b_bord = b_g_pass*ones(1,nb_g_passifs_bord);         % largeur des guides passifs du bord dans le sens toroidal
   b_inter_modules = b_g_pass*ones(1,nb_g_passifs_inter_modules);   % largeur des guides passifs du bord dans le sens toroidal

   b = [b_bord,kron(ones(1,nb_modules_tor-1),[b_module,b_inter_modules]),b_module,b_bord];

   z = zeros(1,nb_g_total_ligne);

   for ind = 2:nb_g_total_ligne

       z(1,ind) = z(1,ind-1) + b(ind-1) + e;
    
   end
   
else
   
   z = zeros(1,nb_g_total_ligne);
   for ind = 2:nb_g_total_ligne

       z(1,ind) = z(1,ind-1) + b(ind-1) + e(ind-1);
    
   end
    
end

h = espacement_g_pol*ones(1,nb_g_pol-1);

y = zeros(1,nb_g_pol);

for ind = 2:nb_g_pol

    y(ind) = y(ind-1) +  h(ind-1) + a;
    
end


% export output scenario is exist
if exist('scenario')
        % store geometrical properties of the antenna into the scenario
        scenario.antenna.a=a;
        scenario.antenna.b=b;
        scenario.antenna.e=e;
        scenario.antenna.z=z;
        scenario.antenna.y=y;
        scenario.antenna.nma_toro=nb_modules_tor; 
%          scenario.antenna.nma_polo= %  TODO

        varargout{1}=scenario;
end