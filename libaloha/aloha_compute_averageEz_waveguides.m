function scenario=aloha_compute_averageEz_waveguides(scenario)
% Compute the average value of the Ez electric field in each waveguides of the antenna.
% 
% scenario = aloha_compute_averageEz_waveguides(scenario)
% 
% INPUT
%  - scenario : ALOHA scenario. (Only the first one is used.)
%  
% OUTPUT 
%  - scenario : ALOHA scenario witch contains the following new 'results' fields :
%       * Ez_average (NbrowsxNwg) : average absolute Ez field in each waveguide of a row.
%       * PassWG_idx(NbrowsxNwg) : index of passive waveguides : 1 for passive, 0 for actives WG
%
% AUTHOR: JH
% LAST UPDATE: nov.2010
% NB : work only for 1D calculation
% 

for idx_sc = 1:length(scenario)
    sc = scenario(idx_sc);

    % if the Efield does not exist in the scenario results
    % We calculate it
    if not(isfield(sc.results, 'Efield'))
        sc = aloha_compute_champEmbouchure1D(sc);
    end
    Efield = sc.results.Efield;

    % get the dimensions and coordinates of the waveguides of the antenna
    if aloha_isAntennaITM(sc)
      disp(aloha_message('ITM antenna description'));
      [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
    else
      disp(aloha_message('Old-fashioned antenna description'));
      [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor,scenario]=aloha_utils_getAntennaCoordinates(architecture,scenario);
    end

    % calculate the average field in a waveguide
    for idx_wg=1:nb_g_total_ligne
        idx_z_wg = find(sc.results.abs_z(1,:)>=z(idx_wg) & sc.results.abs_z(1,:)<=(z(idx_wg)+b(idx_wg)));
        Ez_average(idx_wg) = mean(abs(Efield(1,idx_z_wg)));
    end

    % calculate the passive waveguide position of some known antenna
    switch sc.antenna.architecture
        case 'antenne_C3'
            PassWG_idx = [1:7:nb_g_total_ligne];
        case 'antenne_C4'
            PassWG_idx = [1:2:nb_g_total_ligne];
        otherwise
            PassWG_idx = [];
    end

    % save results in the scenario
    sc.results.Ez_average = Ez_average;
    sc.results.PassWG_idx = PassWG_idx; 
    scenario(idx_sc) = sc;
end % for all scenarios