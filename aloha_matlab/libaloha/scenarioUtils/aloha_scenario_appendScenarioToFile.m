function status = aloha_scenario_appendScenarioToFile(scenario, scenario_filename)
% Save and append a scenario on a .mat file. 
% 
% Results are saved in a matlab structure, 
% which length correspond to the number of scenarii
% computed with the same filename.
%  
% status = aloha_scenario_appendScenarioToFile(scenario, scenario_filename)
% 
%  
% INPUT
%  - scenario [structure] : scenario you wish to save and append to the .mat file
%  - scenario_filename [string] : matlab .mat file.
%  
% OUPUT
%  - status : 0 if OK, 1 if any problem
% 
% NB: if the .mat file is empty, a new .mat file is created.
% If it allready exist, then the scenario is appended to the scenario 
% which is in the .mat file, then saved with the same name.
% 
% AUTHOR: JH
% LAST CHANGES:
%  - 09/2008: added status output
%  - 09/2008: creation
%  

    % load previously saved scenario(s)
    % if the file doesn't exist, a void will be returned
    old_scenario = aloha_scenario_load(scenario_filename);
    idx_current_scenario = length(old_scenario)+1;
    
    %  Append and save the current scenario into the matlab file.
    %  
    %  si la structure actuelle est égale à la structure précédente, 
    %  pas la peine de sauvegarder ! (ie: il s'agit de la meme simulation!)
    if (idx_current_scenario > 1) & isequal(old_scenario(end).results.S_plasma, scenario.results.S_plasma)
        disp(aloha_message('!!! PAS de sauvegarde : la simulation precedente est identique !!!'));
        status = 1;
    else
        disp(aloha_message(['creation de la structure scenario #', num2str(idx_current_scenario), ',', scenario.options.comment]));
        disp(aloha_message('Sauvegarde des resultats.'));
        %  structure concatenation 
        new_scenario = aloha_scenario_append(old_scenario, scenario); 
        % save into .mat file
        status = aloha_scenario_save(new_scenario, scenario_filename);
    end