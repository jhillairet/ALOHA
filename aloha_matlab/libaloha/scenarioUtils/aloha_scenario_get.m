function varargout = aloha_scenario_get(scenario, varargin)
% Return the specified field names of a scenario structure.
% The fields are found into the structure AND the sub-structure which
% are contained in the scenario structure
%  
% [field1, field2, ...] = aloha_scenario_get(scenario, fieldName1, fieldName2, ...)
%  
% EXAMPLE:
%  
%  [ny,nz,dP] = aloha_scenario_get(scenario(1),'ny','nz','dP') ;
%  
%  [dens, gradient] = aloha_scenario_get(scenario, 'ne0', 'dne');
%
% If there is more than one scenario which are stacked into the input scenario structure,
% then the output variable are of dimension length(scenario).
%
% INPUT:
%  - scenario [structure scenario]
%  - fieldName1, fieldName2, ... [string]
% 
% OUTPUT:
%  - field1, field2... [type of the specified fileName1, fieldName2...]
% 
% NB: the returned values are taken from from the 1st scenario structure by default.
% 
% NB2: matrix are returned with size (1,N,M), so don't forget to use the 'squeeze' function !
% 
% NB3: empty is returned each time a field is not found
% 
% AUTHOR: JH
% LAST CHANGE:
%  - 23/09/2008: the function search recursively into all the level of the structure for the searched fields.
%  - 02/09/2008: modify the function to search into the 1st sub-structure contained into the main structure
%  - 08/08/2088: creation
%   

% pour tous les scenarios imbriques dans 
% la structure sceraio passee en argument
for idx_sc=1:length(scenario)
    sc = scenario(idx_sc);
    fn = fieldnames(sc);

    % pour tous les champs que l'on recherche 
    for idx_varargin=1:length(varargin)
        varargout{idx_varargin}(idx_sc,:,:) = get_field_recursive(sc, [], varargin{idx_varargin});
    end
end


% fonction recursive pour chercher un champ dans une structure. 
% Renvoie vide si non trouvee.
% 
function F = get_field_recursive(S, F, fieldName)

    % pour tous les champs contenus dans la structure
    fn = fieldnames(S);
    id_f = 1;
    while((id_f <= length(fn)))
        actualField = getfield(S, fn{id_f});
        % plusieurs possibilites : 
        % - le champ actuel est une structure : 
        %   on cherche alors dans cette nouvelle structure (recursion)
        % - le champ actuel est une variable : 
        %   on teste s'il s'agit du champ que l'on cherche. 
        if isstruct(actualField)
            F = get_field_recursive(actualField, F, fieldName);    
        elseif isfield(S, fieldName)
            if strcmp(fn{id_f}, fieldName) %& ~isempty(actualField)
                F = actualField;
                break; 
            end
        end
        
        id_f = id_f + 1;
    end