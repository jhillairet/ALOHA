function [S_plasma, rac_Zhe] = aloha_getDatasFromBinaryFile(filename, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max, ind,type_swan_aloha)
% extraction des donnees du fichier S_plasma.dat
% 
% 

    tableau = zeros(1,nb_g_pol);
    tableau(ind) = 1;     

    fichier_data = fopen(filename,'r');

    fread(fichier_data,1,'int32');
    
    S_ligne = [];
    for i_ind = 1:(Nme+Nmh)*nb_g_total_ligne
        for j = 1:(Nme+Nmh)*nb_g_total_ligne
            S_ligne(i_ind,j) = fread(fichier_data,1,'double') + i*fread(fichier_data,1,'double');
           if (abs(i_ind-j)>((Nme+Nmh)*(D_guide_max-5)))
               S_ligne(i_ind,j) =0;
           end
        end
    end

    rac_Zhe_ligne = [];
    fread(fichier_data,1,'double');
    for i_ind = 1:(Nme+Nmh)*nb_g_total_ligne
        for j = 1:(Nme+Nmh)*nb_g_total_ligne
            rac_Zhe_ligne(i_ind,j) = fread(fichier_data,1,'double') + i*fread(fichier_data,1,'double');
        end
    end
    

           
    if (type_swan_aloha == 0)
    
        Zc_swan = [];
        for i_ind = 1:nb_g_total_ligne
            Zc_swan = [Zc_swan,120*pi,-i*120*pi.*sqrt((1:Nme)*((celerite/(2*b(i_ind)*freq))^2-1))];
        end
    Z_ligne = rac_Zhe_ligne*(eye((Nme+Nmh)*nb_g_total_ligne)+S_ligne)*inv(eye((Nme+Nmh)*nb_g_total_ligne)-S_ligne)*rac_Zhe_ligne;
    S_ligne = diag(sqrt(1./Zc_swan))*(Z_ligne-diag(Zc_swan))*inv(Z_ligne+diag(Zc_swan))*diag(sqrt(Zc_swan));
    
    end
     
    S_plasma = S_plasma + kron(diag(tableau),S_ligne);
    rac_Zhe = rac_Zhe + kron(diag(tableau),rac_Zhe_ligne);