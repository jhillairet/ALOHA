      % TODO : fix 
      % when using an ITM antenna architecture description,
      % we have to translate the parameters into the old-fashioned way
      scenario.antenna.freq = scenario.antenna_lh.frequency;% TODO : correct this bug
      freq = scenario.antenna_lh.frequency; % TODO : correct this bug 
      nb_g_pol = aloha_scenario_get(scenario, 'nma_theta');
      nb_modules_tor = aloha_scenario_get(scenario, 'nma_phi');
      nb_g_module_pol = aloha_scenario_get(scenario, 'nwm_theta');
      nb_g_module_tor = aloha_scenario_get(scenario, 'nwm_phi');
      a = aloha_scenario_get(scenario, 'hw_theta');
      nb_g_passifs_bord = aloha_scenario_get(scenario, 'npwe_phi');
      pass_module_tor = reshape(find(aloha_scenario_get(scenario, 'mask') == 0),1,[]); % the reshape force the array to be horizontal
      nb_g_passifs_inter_modules = aloha_scenario_get(scenario, 'npwbm_phi');
      antenne_standard = 0;
      lcc = squeeze(aloha_scenario_get(scenario, 'scl'))';
      chemin_retour = pwd;
      chemin_aller = scenario.antenna_lh.setup.modules.Sparameters.pathTo;
      nom_fichiers = scenario.antenna_lh.setup.modules.Sparameters.SFileNames;
      phase_rallonge = scenario.antenna_lh.setup.modules.Sparameters.phase_deembedded;

% Les parametres utiles pour S_plasma1D ont ete deplace dans aloha_init
