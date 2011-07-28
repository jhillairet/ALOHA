function [S_plasma, rac_Zhe] = aloha_getDatasFromAsciiFile(filename, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max, ind,type_swan_aloha, architecture, freq, scenario)
%  data extraction from ASCII file S_plasma.dat (fortran90 output)
%  
%  EXAMPLE:
%  [S_plasma, rac_Zhe] = aloha_getDatasFromAsciiFile(...
%          filename, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max, ind,type_swan_aloha, architecture, freq)
%  
%  INPUTS:
%  - filename : 
%  - S_plasma : matrice S plasma
%  - rac_Zhe : 
%  - Nme, 
%  - Nmh, 
%  - nb_g_pol, 
%  - nb_g_total_ligne, 
%  - D_guide_max, 
%  - ind,
%  - type_swan_aloha
%  - architecture
%  - freq
%  OUTPUTS:
%  - S_plasma, 
%  - rac_Zhe]
%  
%  
%  AUTHOR(S) : DV,JH
%  
%  LAST UPDATE : 
%   - 02/07/2008 [creation]
%   - 18/06/2010 [add some input arguments to perform the SWAN->ALOHA impedance renormalization]
%  
% 


    tableau = zeros(1,nb_g_pol);
    tableau(ind) = 1;     

    [S_r, S_i, RZ_r, RZ_i, K_r, K_i]=textread('S_plasma2.dat', '(%f,%f) (%f,%f) (%f,%f)', 'headerlines', 2);

    S_plasma_inline = complex(S_r,S_i);
    S_plasma_ligne2 = reshape(S_plasma_inline, sqrt(length(S_plasma_inline)), sqrt(length(S_plasma_inline)));    


    rac_Zhe_inline= complex(RZ_r,RZ_i); 
    rac_Zhe_ligne2 = reshape(rac_Zhe_inline, sqrt(length(rac_Zhe_inline)), sqrt(length(rac_Zhe_inline)));


    %% conversion from SWAN->ALOHA
    % Waveguide impedance normalization definition is not the same for the two code:
    %  - in SWAN, this is parallel plate waveguide (=vacuum impedance for the TEM mode)
    %  - in ALOHA, this is rectangular waveguide
    % So, in order to compare scattering parameters, one has to renormalize them to the correct impedance.
    % In this case, we renormalize to SWAN impedances.
    if (type_swan_aloha == 0)
        % load physical constants in workspace
        aloha_constants  
        % load antenna parameters in workspace
        if aloha_isAntennaITM(scenario)
            disp(aloha_message('ITM antenna description'));
            [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
        else
            disp(aloha_message('Old-fashioned antenna description'));
            [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor,scenario]=aloha_utils_getAntennaCoordinates(architecture,scenario);
        end
    
        Zc_swan = [];
        for i_ind = 1:nb_g_total_ligne
            % waveguide impedance for parallel plate waveguide is Z0 for TEM and value below for TM
            Zc_swan = [Zc_swan,120*pi,-i*120*pi.*sqrt((1:Nme)*((celerite/(2*b(i_ind)*freq))^2-1))];
        end
        Z_ligne = rac_Zhe_ligne2*(eye((Nme+Nmh)*nb_g_total_ligne)+S_plasma_ligne2)*inv(eye((Nme+Nmh)*nb_g_total_ligne)-S_plasma_ligne2)*rac_Zhe_ligne2;
        S_plasma_ligne2 = diag(sqrt(1./Zc_swan))*(Z_ligne-diag(Zc_swan))*inv(Z_ligne+diag(Zc_swan))*diag(sqrt(Zc_swan));
    
    end

    S_plasma = S_plasma + kron(diag(tableau), S_plasma_ligne2);
    rac_Zhe = rac_Zhe + kron(diag(tableau), rac_Zhe_ligne2);