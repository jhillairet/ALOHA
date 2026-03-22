% Load each field of the stucture named 'scenario' into the workspace
%  
%  
%  
%  aloha_scenario_loadIntoWorkspace
%  
%  NB : works only if a 'scenario' structure exist in the matlab session!
  
if ~exist('scenario', 'var')
    error('scenario structure doesn''t exist !');
end
  
fn = fieldnames(scenario);

for idx_fn=1:length(fn)
    subfn = fieldnames(scenario.(fn{idx_fn}));
    for idx_subfn=1:length(subfn)
        str = [subfn{idx_subfn},'=','scenario.' fn{idx_fn} '.' subfn{idx_subfn} ';'];
        eval(str);
    end
end

if bool_lignes_identiques
  dne0 = dne(1); 
  dne1 = dne(2);
else
  dne0 = dne(:,1); 
  dne1 = dne(:,2);
end

Nmh = modes(1);
Nme = modes(2);