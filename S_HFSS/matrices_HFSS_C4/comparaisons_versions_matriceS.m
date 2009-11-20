% Comparaison des matrices S des antennes. 
% Ces matrices S sont obtenues avec HFSS, par J.Belo.
%
% Auteur : J.Hillairet 
% MAJ : 14/12/07
%

% version 1
Smatrix_mfile_v1 = strvcat('S_C4_24b', ...
			'S_C4_23b', ...
			'S_C4_22b', ...
			'S_C4_21b', ...
			'S_C4_14b', ...
			'S_C4_13b', ...
			'S_C4_12b', ...
			'S_C4_11b');  
% version 2 (J.Belo, le 13/12/07)
Smatrix_mfile_v2 = strvcat('Module_bas_plaque_long7_ext2', ...
                  'Module_bas_plaque_long6_long7', ...     
                  'Module_bas_plaque_long5_long6', ...                
                  'Module_bas_plaque_long4_long5', ...                
                  'Module_bas_plaque_long3_long4', ...                           
                  'Module_bas_plaque_long2_long3', ...                
                  'Module_bas_plaque_long1_long2', ...                
                  'Module_bas_plaque_ext1_long1');
		  			
Smatrix_mfile_v1 = strvcat('S_C4_24h', ...
			'S_C4_23h', ...
			'S_C4_22h', ...
			'S_C4_21h', ...
			'S_C4_14h', ...
			'S_C4_13h', ...
			'S_C4_12h', ...
			'S_C4_11h');    	  
Smatrix_mfile_v2 = strvcat('Module_haut_plaque_long7_ext2', ...
                  'Module_haut_plaque_long6_long7', ...     
                  'Module_haut_plaque_long5_long6', ...                
                  'Module_haut_plaque_long4_long5', ...                
                  'Module_haut_plaque_long3_long4', ...                           
                  'Module_haut_plaque_long2_long3', ...                
                  'Module_haut_plaque_long1_long2', ...                
                  'Module_haut_plaque_ext1_long1');
		  

for idx=1:size(Smatrix_mfile_v1,1)
  % ouverture des fichiers et conversion des variables
  eval(Smatrix_mfile_v1(idx,:));
  S_v1 = reshape(S, sqrt(size(S,2)), sqrt(size(S,2)));
  eval(Smatrix_mfile_v2(idx,:));
  S_v2 = squeeze(S);

  clf(figure(idx))
    set(gcf, 'Name', [Smatrix_mfile_v1(idx,:), ' VS ', Smatrix_mfile_v2(idx,:)]);
    subplot(121)	
    pcolor(abs(S_v1) - abs(S_v2)), axis equal; colorbar; axis([1 7 1 7]);
    title('|S_{ij}^{old}| - |S_{ij}^{new}|')

    subplot(122)
    pcolor(180/pi*(angle(S_v1)) - 180/pi*angle(S_v2)), axis equal; colorbar; axis([1 7 1 7]);
    title('Phase(S_{ij}^{old}) - Phase(S_{ij}^{new})')

end
