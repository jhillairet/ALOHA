function scenarios = aloha_scenario(scenarios)
% Code de couplage ALOHA pour le calcul du couplage d'une antenne hybride avec un plasma de Tokamak.
% 
% EXAMPLE
%  scenario = aloha_scenario(scenario)
% 
% INPUT
%  - scenario <structure length=N>: structure scenario ALOHA de taille N 
% 
% OUPUT
%  - scenario <structure length=N>: structure passee en argument, 
%               mais contenant en plus les resultats de calcul 
%               effectues par ALOHA (matrices S plasma, etc...)
%               
% AUTHOR: S.Berio / O.Izacard / D.Voyer / J.Hillairet
% LAST UPDATE:
%  - 09/2008 : mise sous forme d'une fonction et 
%              traitement des donnees sous forme de scenario(s)
%  


% ##################################################################
% ALOHA initialisation
% ##################################################################
aloha_init; 
% pour tous les scenarios contenus dans la structure
for idx_sc = 1:length(scenarios)
    tic; % for computing time display
    scenario = scenarios(idx_sc);
    disp(aloha_message(['########## Scenario ', num2str(idx_sc),'/', num2str(length(scenarios)),' ##########']));
  
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Compatibility section
    % 
    % load each field of the stucture scenario into the workspace
    % (for compatibility reasons, it's simpler than to change evvery variables in the code...)
    aloha_scenario_loadIntoWorkspace
    
    % variables classiques
    if bool_lignes_identiques
      dne0=dne(1); 
      dne1=dne(2); 
      lambda_n0 = lambda_n(1);
      lambda_n1 = lambda_n(2);    
    else
      dne0=dne(:,1); 
      dne1=dne(:,2); 
      lambda_n0 = lambda_n(:,1);
      lambda_n1 = lambda_n(:,2);
    end
    Nmh = scenario.options.modes(1); 
    Nme = scenario.options.modes(2);

    % test the validity of the plasma gradient variable. 
    % they must be consistant, i.e. the dne value must be related to lambda_n
    if bool_lignes_identiques
        if (dne0 ~= ne0./lambda_n0)
            disp(aloha_message('WARNING: Non consistent plasma gradient dne0 definition ! forcing value using lambda_n as reference : '));
            dne0 = ne0./lambda_n0;
            disp(aloha_message(['WARNING: lambda_n0=', num2str(lambda_n0)]));
            disp(aloha_message(['WARNING: ==> dne0=', num2str(dne0)]));
        end
    
        if (dne1 ~= (1+d_couche./lambda_n0).*ne0./lambda_n1)     
            disp(aloha_message('WARNING: Non consistent plasma gradient dne1 definition ! forcing value using lambda_n as reference : '));
            dne1 = (1+d_couche./lambda_n0).*ne0./lambda_n1;
            disp(aloha_message(['WARNING: lambda_n1=', num2str(lambda_n1)]));
            disp(aloha_message(['WARNING: ==> dne1=', num2str(dne1)]));
        end
    else
        for idx=1:length(ne0)
            if (dne0(idx) ~= ne0(idx)./lambda_n0(idx))
                disp(aloha_message('WARNING: Non consistent plasma gradient dne0 definition ! forcing value using lambda_n as reference : '));
                dne0(idx) = ne0(idx)./lambda_n0(idx);
                disp(aloha_message(['WARNING: lambda_n0=', num2str(lambda_n0(idx))]));
                disp(aloha_message(['WARNING: ==> dne0=', num2str(dne0(idx))]));
            end
        
            if (dne1(idx) ~= (1+d_couche./lambda_n0(idx)).*ne0(idx)./lambda_n1(idx))     
                disp(aloha_message('WARNING: Non consistent plasma gradient dne1 definition ! forcing value using lambda_n as reference : '));
                dne1(idx) = (1+d_couche./lambda_n0(idx)).*ne0(idx)./lambda_n1(idx);
                disp(aloha_message(['WARNING: lambda_n1=', num2str(lambda_n1(idx))]));
                disp(aloha_message(['WARNING: ==> dne1=', num2str(dne1(idx))]));
            end
        end
    end
    % MODIF JH 23/10/2009
    % Jua found that the gradient values dne wasn't save into the scenario 
    % It has an impact on the wat the 1D spectrum is computed, because we need 
    % the dne value into this calculation. 
    % So, we must update the scenario file: 
    scenario.plasma.dne = [dne0 dne1];


    % display main plasma parameters
    if bool_lignes_identiques
        disp(aloha_message(['Main plasma parameters : ', ...
                        'ne0=', num2str(ne0'/1e17) , 'x10^17 m^-3, '...
                        'lambda1=', num2str(lambda_n0*1e3), ' mm. ' ...
                        '(lambda2=', num2str(lambda_n1*1e3), ' mm.) ']));
    else
        disp(aloha_message('Main plasma parameters : '));
        for idx=1:length(ne0)
         disp(aloha_message(['row #', num2str(idx), ' ne0=', ...
                         num2str(ne0(idx)'/1e17) , 'x10^17 m^-3, '...
                        'lambda1=', num2str(lambda_n0(idx,:)*1e3), ' mm. ' ...
                        '(lambda2=', num2str(lambda_n1'*1e3), ' mm.) ']));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % load physicals constants
    aloha_constants;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % architecture, geometrie de l'antenne
    disp(aloha_message(['Take into account antenna architecture : ', architecture]));
    % read the architecture as a function which should return a structure 
    % in the ITM antenna_lh CPO. If there is an error in this attempt, this should mean
    % that the file is not function but a script, as it was before ALOHA ITM compliance
    if aloha_isAntennaITM(scenario)
      disp(aloha_message('ITM antenna description')); 
      aloha_utils_ITM2oldAntenna; % convert ITM data to old fashioned ALOHA parameters
	
    else
      disp(aloha_message('Old-fashioned antenna description'));
      eval(architecture);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % antenna's module loading
    disp(aloha_message('Module''s excitation:'))
    disp(a_ampl.*exp(i*a_phase));

	
    %%%%%%%%%%%%%%%%%%%%%%%%
    % coordonnees des guides
    
    disp(aloha_message('Take into account waveguide coordinates'));
    %coordonnees_guides;   
    % if the field antenna_lh exist, this means that we use the ITM ALOHA version
    if aloha_isAntennaITM(scenario)
      disp(aloha_message('ITM antenna description'));
      [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
    else
      disp(aloha_message('Old-fashioned antenna description'));
      [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor,scenario]=aloha_utils_getAntennaCoordinates(architecture,scenario);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul de la matrice S du plasma
    % pour le cas 1D ou 2D
    disp(aloha_message('Compute plasma Scattering matrix'));
    switch version_code
        case '1D'
            S_plasma_1D;    

        case '2D'
            S_plasma_2D;
        
        otherwise
            error('Version_code not (correctly) defined !');
    end
    scenario.options = aloha_setfield(scenario.options, chemin_binaire_fortran); 
    
    % Save plasma scattering parameters into the scenario
    scenario.results = aloha_setfield(scenario.results, S_plasma, rac_Zhe, K_cpl); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % voies passives et voies actives
    
    disp(aloha_message('Take into account passive waveguides (if any)'));
    voies_actives_passives;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul de la matrice S antenne
    disp(aloha_message('Take into account Antenna Scattering matrices'));
    S_antenne;
     % save global antenna S parameters into the scenario
    scenario.results = aloha_setfield(scenario.results, S_ant_11, S_ant_12, S_ant_21, S_ant_22);

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculte the coupling response of the antenna 
    disp(aloha_message('Compute antenna/plasma interactions'));
    %reponse_antenne;
    scenario = aloha_compute_RC(scenario)
          
    % show the reflexion coefficient 
    disp(aloha_message('Reflexion Coefficients (RC) per module :'));
    disp(scenario.results.RC);    
    
    % save some constants into the scenario (for check purpose essentially)
    scenario.results = aloha_setfield(scenario.results, k0);

    % show and save execution time
    disp(aloha_message(['Execution time : ', num2str(toc), ' s']));
    scenario.results.computationTime = toc;

    % moyenne du coefficient de reflexion en puissance
    R = mean( scenario.results.RC ) ;
    disp(aloha_message(['Average reflexion coefficient |Gamma^2| : R = ', num2str(R), ' %']));  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    % Antenna spectrum 
    if (bool_compute_spectrum)
        try
            disp(aloha_message('Compute of the power spectrum'));  
            tic;
            switch version_code
                case '1D'
                    % compute all the spectrum parameters : ny,nz,dP,dP_nz
                    scenario=aloha_compute_spectrum1D(scenario);                  

                case '2D'
                    % retrieve the spectrum from the fortran output file
                    scenario=aloha_compute_spectrum2D(scenario); 
            end
            toc        
            if (bool_display_spectrum)
                disp(aloha_message('Display the power spectrum'));
                aloha_plot_spectrum(scenario);
            end
                        
            % check the power conservation
            scenario = aloha_compute_powerConservation(scenario);
            
        catch
            s = lasterror;
            disp(aloha_message(['ERROR:',s.message]));
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna directivity    
    % NB : Directivity require the computation of the spectrum.
    %
    if (bool_compute_directivity & bool_compute_spectrum)
        try
            disp(aloha_message('Compute the antenna directivity'));
            scenario=aloha_compute_directivity1D(scenario);
        
            if (bool_display_directivity)
                disp(aloha_message('Display the antenna directivity'));
                aloha_plot_directivity(scenario);
            end
        catch
            s = lasterror;
            disp(aloha_message(['ERROR:',s.message]));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % figure champ total
    % Computation of the electric field into the mouth of the antenna
    % 
    if (bool_compute_total_field)
        try 
            disp(aloha_message('Compute the parallel E field in the mouth'));
            switch version_code
                case '1D'
                    % JH 12/09/2011
                    % Changed aloha_compute_champEmbouchure1D to aloha_compute_champEmbouchure1D_v2
                    scenario=aloha_compute_champEmbouchure1D_v2(scenario);
                case '2D'
                    scenario=aloha_compute_champEmbouchure2D(scenario);
            end
        
            if (bool_display_total_field)
                disp(aloha_message('Display the parallel E field in the mouth'));
                aloha_plot_champEmbouchure1D(scenario);
            end
        catch
            s = lasterror;
            disp(aloha_message(['ERROR:',s.message]));
        end        
    end
    
    % Propagation of the electric field into the plasma.
    % ATTENTION: this is a plane-wave propagation algorithm, 
    % which does not take into account hot plasma effect, collisions, absorption, etc...
    % In order to give the right picture of the LH propagation, one as to run a 
    % fokker-plank code.
    %  %%%%%%%%%%%%%%%%%%%%%%%%
    % figure champ plasma
    if (bool_compute_plasma_field)
        try
            disp(aloha_message('Compute the parallel E field into the plasma'));
            switch version_code
                case '1D'
                    scenario=aloha_compute_champPlasma1D(scenario);
                case '2D'

            end

            if (bool_display_plasma_field)
               disp(aloha_message('Display the parallel E field into the plasma')); 
               aloha_plot_champPlasma1D(scenario);
            end
        catch
            s = lasterror;
            disp(aloha_message(['ERROR:',s.message]));
        end  
    end
   
    disp(aloha_message('Saving results into the scenario.'));
    scenarios(idx_sc) = scenario;
end



