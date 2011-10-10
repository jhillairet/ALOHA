% creation d'un "batch" de scenarios pour ALOHA.
% 
% Il s'agit dans cet exemple de calculer la reponse du plasma a une excitation pour 
% un scan en densite electronique ne0, et pour 3 different types de scenario. 
% Par exemple ici, 3 type de profil de densite sont utilises 
%  - 1 gradient de 2cm, 
%  - 1 gradient de 2mm 
%  - et un double gradient de 2mm/cm. 
% 
% NB: cet exemple suppose qu'il existe dans le PATH de matlab des fichiers
% portant les noms de  : 
%  'scenario_C2_2cm.m' 
%  'scenario_C2_2mm.m' 
%  'scenario_C2_2mm_2cm.m'
% contenant la description de ces scenarios.
% 
% AUTHOR:JH
% LAST UPDATE
%  - 10/2008: creation
clear all;

% on veut faire un scan en densite sur les valeurs suivantes :
dens = [1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 20, 30, 40]*1e17;

% Dans la boucle qui suit, on genere des scenarios de type different.
% Pour cela, on appelle a chaque iteration une fonction differente, 
% qui genere le scenario 'type'. Les noms des scenarios type sont definis ici: 
% 
% NB : il faut bien sur s'assurer que les fichiers scenario_C2_2cm.m, 
% scenario_C2_2mm, ... existent dans le path.
scen{1} = 'tutorial_aloha_scenario_simple_grill_8waveguides';


for id_scen = 1:length(scen)
    disp(['-----------', num2str(id_scen),'/',num2str(length(scen)),'------']);
    % on charge un scenario "type". 
    % NB : On s'est assure au prealable que tous les parametres autre que
    % la densite sont OK.
    sc = eval(scen{id_scen});
    % on reproduit ce scenario pour le nombre de valeur de densite:
    sc = repmat(sc, length(dens), 1);
    
    % on definit la valeur de densite ne0 dans les scenario 
    for id=1:length(dens)
        sc(id).plasma.ne0 = dens(id);
    end  
    
    % on lance ALOHA sur l'ensemble de ces memes scenarios 'type' 
    sc = aloha_scenario(sc);
    
    % on sauve le resultat d'un type de scenario
    aloha_scenario_save(sc, [scen{id_scen}, '.mat']);

end % id_scen

