
clear all;

% on veut faire un scan en densite sur les valeurs suivantes :
b_a = [10:11]*1e-3;

scen{1} = 'scenario_PA_ITM';


for id_scen = 1:length(scen)
    disp(['-----------', num2str(id_scen),'/',num2str(length(scen)),'------']);
    % on charge un scenario "type". 
    % NB : On s'est assure au prealable que tous les parametres autre que
    % la densite sont OK.
    sc = eval(scen{id_scen});
    % on reproduit ce scenario pour le nombre de valeur de densite:
    sc = repmat(sc, length(b_a), 1);
    
    % on definit la valeur de densite ne0 dans les scenario 
    for id=1:length(b_a)
        sc(id).antenna_lh.setup.modules.waveguides.bwa = b_a(id);
    end  
    
    % on lance ALOHA sur l'ensemble de ces memes scenarios 'type' 
    sc = aloha_scenario(sc);
    
    % on sauve le resultat d'un type de scenario
    aloha_scenario_save(sc, [scen{id_scen}, '.mat']);

end % id_scen

