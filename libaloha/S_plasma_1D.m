%  Calcule la matrice S_plasma
%  
%  
% TODO : option code 1D/2D
%  
% pour ttes les versions 1D
% calcul Ys(nz) 1D periodicite (T_grill), reciprocite et decouplage (D_guide_max)
%  
% version 3
% profil lineaire, methode des residus + domaine d'integration limite e t = 1/(max_nz+1)
%  
% version 4,5 et 6 
% profil lineaire, pas de methode des residus, pertes dans le milieu, domaine d'integration limite e max_nz
%  
% version 4
% profil lineaire, une couche de plasma
%  
% version 5 
% profil lineaire, couche de vide + couche de plasma
%  
% version 6
% profil lineaire, deux couches de plasma
%  
% VErsion 7
% plasma step 

%  Determine le nom du binaire a utiliser
% TODO : use tempname for generate a temporary file - fortran bin should accept args 
binary_name = aloha_getBinaryName(aloha_getArchitecture, version);
% noms des fichiers d'échange avec le code fortran
fortranInputAsciiFileName = 'par_grill.dat';
fortranOutputBinaryFileName = 'S_plasma.dat';
fortranOutputAsciiFileName = 'S_plasma2.dat';

%  borne d'integration nz
max_nz = 100; 

%%%%%%%%%%%%%%%%%%%%%%%
k0 = 2*pi*scenario.antenna.freq/celerite;


%%%%%%%%%%%%%%%%%%%%%%%
% Changement de repertoire :
% - pour creer le fichier d'argument du code en fortran
% - lancer le code fortran
% - recuperer les resultats
% avant de revenir au repertoire actuel.
% 
% En cas d'erreur lors de l'execution du binaire, le programme affiche une alerte et continu.
chemin_retour = pwd;
chemin_binaire_fortran = '/code_1D/couplage_1D';

% change dir to the fortran binary directory
cd([aloha_utils_getRootPath, chemin_binaire_fortran]);

if (bool_lignes_identiques)
    ne0 = ones(1,nb_g_pol)*ne0;
    dne0 = ones(1,nb_g_pol)*dne0;
    % Rajout le 16/05/2007 par Izacard Olivier
    dne1 = ones(1,nb_g_pol)*dne1;
    d_couche = ones(1,nb_g_pol)*d_couche;
end


nc = ((2*pi*freq)^2)*me*Eps0/(qe^2);
X0 = ne0./nc;
D0 = k0*ne0./dne0;

S_plasma = zeros(nb_g_total_ligne*nb_g_pol*(Nmh+Nme),nb_g_total_ligne*nb_g_pol*(Nmh+Nme));
rac_Zhe = zeros(nb_g_total_ligne*nb_g_pol*(Nmh+Nme),nb_g_total_ligne*nb_g_pol*(Nmh+Nme));


disp(aloha_message(['Coupling calculation ; version ', num2str(version)]));

if bool_debug
    disp('Appuyez sur une touche pour lancer le calcul');
    pause;
end

% for all poloidal rows
for ind = 1:nb_g_pol
    textInd = ['[row #', num2str(ind),'/', num2str(nb_g_pol),'] : '];
    % execute binary only if:
    %  - this is the first row and bool_lignes_identiques == true
    %  - bool_lignes_identiques == false (then do for all rows)
    if ((bool_lignes_identiques == true) & (ind > 1))
        disp(aloha_message([textInd, '[bool_lignes_identiques=true]']));
        disp(aloha_message([textInd, '--> No need to calculate this row! taking previous result ']));
    else
        % 
        % Writing the ascii file containing binary parameters 
        % This file depends on the binary version
        % 
        if (version == 3)

            tampon1 = ne0(ind);
            tampon2 = dne0(ind);

            varCell = {Nmh, Nme, freq, tampon1, tampon2,  ...
                    nb_g_total_ligne, a, b, z, T_grill, ...
                    D_guide_max, erreur_rel, max_nz}; 

        elseif (version == 4)

            tampon1 = ne0(ind);
            tampon2 = dne0(ind);
            varCell = {Nmh, Nme, freq, tampon1, tampon2, ...
                    nb_g_total_ligne, a, b, z, T_grill, ...
                    D_guide_max, erreur_rel, pertes, max_nz};

        elseif (version == 5)

            tampon1 = ne0(ind);
            tampon2 = dne0(ind);

            varCell = {Nmh, Nme, freq, tampon1, tampon2, ...
                    nb_g_total_ligne, a, b, z, T_grill, ...
                    D_guide_max, erreur_rel, pertes, max_nz, d_vide}; 

        elseif (version == 6)

            ne1 = ne0+dne0.*d_couche;            % densite electronique couche ne2
            %if (lignes_identiques == 1)
            %    ne1 = ones(1,nb_g_pol)*ne1;% Modif du 10/05/2007
            %    dne1 = ones(1,nb_g_pol)*dne1;	
            %end
            X1 = ne1./nc;
            D1 = k0*ne1./dne1;
            tampon1 = ne0(ind);
            tampon2 = dne0(ind);
            tampon3 = d_couche(ind);
            tampon4 = dne1(ind);  

            varCell = {Nmh, Nme, freq, tampon1, tampon2, ...
                    tampon3, tampon4, nb_g_total_ligne, a, b, z, ...
                    T_grill, D_guide_max, erreur_rel, pertes max_nz};

        elseif (version == 7)

            tampon1 = ne0(ind);
            tampon2 = dne0(ind);

            varCell = {Nmh, Nme, freq, tampon1, tampon2, ...
                    B0_version7, nb_g_total_ligne, a, b, z, ...
                    T_grill, D_guide_max, erreur_rel, pertes, max_nz};
        end

        disp(aloha_message([textInd,'Writing parameters file for the binary']));
        aloha_saveFortranInputFile(fortranInputAsciiFileName, varCell);

        %
        % Fortran Binary execution
        % 
        try
            disp(aloha_message([textInd, 'Run binary ', binary_name]));
            [status,result] = system(['./',binary_name]);

            if bool_debug
                disp(result);
            end

            if (status ~= 0)
                disp(aloha_message([textInd, 'Binary execution problem : ']));
                error(aloha_message(result));
            end
        catch
            cd(chemin_retour);
            error(aloha_message('?! Binary execution problem ?!'));
        end
    end
    % determine quel fichier de resultat
    % utiliser en fonction de l'architecture
    % Les resultats sont ecrits dans des fichiers binaires lorsqu'on utilise
    % les binaires F77, soit sur la plateforme alpha (deneb)
    % Autrement, il s'agit de fichier ascii
    if strcmp(aloha_getArchitecture, 'alpha')
        disp(aloha_message([textInd, 'Reading binary result file']));
        [S_plasma, rac_Zhe] = aloha_getDatasFromBinaryFile(fortranOutputBinaryFileName, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max,ind,type_swan_aloha, architecture, freq);
    else
        disp(aloha_message([textInd, 'Reading ascii result file']));
        [S_plasma, rac_Zhe] = aloha_getDatasFromAsciiFile(fortranOutputAsciiFileName, S_plasma, rac_Zhe, Nme, Nmh, nb_g_pol, nb_g_total_ligne, D_guide_max,ind,type_swan_aloha, architecture, freq, scenario);
    end

    if isempty(S_plasma)
        error(aloha_message([textInd, 'S_plasma is empty !!']));
    else
        disp(aloha_message([textInd, 'Results reading: OK']));

        if bool_debug
            disp(aloha_message(['Sum(S_plasma)=', num2str(sum(S_plasma(:)))]));
            disp(aloha_message(['Sum(rac_Zhe)=', num2str(sum(rac_Zhe(:)))]));
        end
    end
    
end % for ind = 1:nb_g_pol

% back to working dir
cd(chemin_retour);

