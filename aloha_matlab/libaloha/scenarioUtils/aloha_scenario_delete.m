function scenario = aloha_scenario_delete(scenario, index)
%  Delete one (or more) scenarios contained in a scenario structure.
%  
%  scenario = aloha_scenario_delete(scenario, index)
%  
%  INPUT :
%   - scenario [struct]: scenario structure of dim=N
%   - index [array of int]: index of the scenario(s) to delete. dim=L
%  
%  OUTPUT :
%   - scenario [struct]: scenario structure of dim=N-L
%   if all scenarii have been delete, or if there was an error, 
%   the returned structure is empty. 
%   
%   AUTHOR : JH
%   LAST UPDATE : 
%    - 07/08/2008 : creation
%  

    scenario(index)=[];
    