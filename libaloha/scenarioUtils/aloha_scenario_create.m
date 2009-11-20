%  Creation d'un scenario a partir des variables contenues dans la session en cours
%  
%  

%%%% sous-groupe 'antenna'
antenna = aloha_setfield([], ...
            architecture, freq, ...
            a_ampl, a_phase);

%%%% sous-groupe 'plasma'
plasma = aloha_setfield([], ne0, d_couche, d_vide);
plasma.lambda_n = [lambda_n0, lambda_n1];
plasma.dne = [dne0, dne1];
plasma.B0 = B0_version7;
plasma.version = version;

%%%%% sous-groupe 'options'
options.scenario_date = date;
options.modes = [Nmh, Nme];
options = aloha_setfield(options, ...
            version_code, ...
            aloha_path, type_swan_aloha, definition_directivite, ...
            comment, scenario_filename, ...
            choc, tps_1, tps_2, TSport, ...
            nz_min, nz_max, dnz, ny_min, ny_max, dny, nbre_ny, nbre_nz, z_coord_min, z_coord_max, nbre_z_coord, x_coord_max, nbre_x_coord, ...
            pas_nz_fig_plasma, fig_Ez_ou_EzHy, lig_fig_plasma);
% ajoute toutes les variables booleennes
% grace a une boucle (pour etre sur de ne pas en oublier une!)
list_bool = whos('bool_*') ;
for ind=1:length(list_bool)
    options = setfield(options, list_bool(ind).name, eval(list_bool(ind).name));
end


% empty results sub-group
results = struct('S_acces', [], 'S_plasma', []);

%%%% scenario
scenario = aloha_setfield([], antenna, plasma, options, results);
