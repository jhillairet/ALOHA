function [S_plasma, rac_Zhe] = aloha_getDatasFromAsciiFile(filename, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max, ind,type_swan_aloha)
%  extraction des donnees du fichier S_plasma.dat
%  
%  EXAMPLE:
%  [S_plasma, rac_Zhe] = aloha_getDatasFromAsciiFile(...
%          filename, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max, ind,type_swan_aloha)
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
%  
%  OUTPUTS:
%  - S_plasma, 
%  - rac_Zhe]
%  
%  
%  AUTHOR(S) : DV,JH
%  
%  LAST UPDATE : 
%   - 02/07/2008 [creation]
%  
% 


    tableau = zeros(1,nb_g_pol);
    tableau(ind) = 1;     

    [S_r, S_i, RZ_r, RZ_i, K_r, K_i]=textread('S_plasma2.dat', '(%f,%f) (%f,%f) (%f,%f)', 'headerlines', 2);

    S_plasma_inline = complex(S_r,S_i);
    S_plasma_ligne2 = reshape(S_plasma_inline, sqrt(length(S_plasma_inline)), sqrt(length(S_plasma_inline)));    
    S_plasma = S_plasma + kron(diag(tableau), S_plasma_ligne2);

    rac_Zhe_inline= complex(RZ_r,RZ_i); 
    rac_Zhe_ligne2 = reshape(rac_Zhe_inline, sqrt(length(rac_Zhe_inline)), sqrt(length(rac_Zhe_inline)));
    rac_Zhe = rac_Zhe + kron(diag(tableau), rac_Zhe_ligne2);


