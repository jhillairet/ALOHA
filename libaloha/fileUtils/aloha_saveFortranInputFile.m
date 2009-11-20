function aloha_saveFortranInputFile(fileName, varCell)
% Ecrit les valeurs contenues dans une {Cell}
% dans un fichier texte.
% Cette fonction reproduit le comportement de la fonction 
%   save -ascii 
% de matlab, car des problèmes ont été observés avec celle-ci
% (inversion des lignes erratique...)
%
%  EXAMPLE:
% aloha_saveFortranInputFile(fileName, varCell)
%  
%  INPUTS:
%  - fileName (str) : nom du fichier texte 
%  - varCell {cell} : cell contenant tous les valeurs a enregistrer
%  
%  
%  OUTPUTS: none
%  
%  AUTHOR(S) : JH
%  
%  LAST UPDATE : 
%   - 02/07/2008 [creation]
%  
% 

    fid = fopen(fileName, 'w');

    for ind=1:length(varCell)
        fprintf(fid, '  %1.7e', varCell{ind});
        fprintf(fid, '\n');
    end

    fclose(fid);