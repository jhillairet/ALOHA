function S = aloha_setfield(S, varargin)
%  Ajoute ou definit des champs a une structure. 
%  Les champs ont automatiquement le nom des variables passees en argument.
%  
%  struct = aloha_setfield(struct, varargin)
%  
%  INPUT: 
%   - S : matlab structure
%   - varargin : variables to save into the structure S
%  
%  OUTPUT:
%   - S : matlab structure
%   
%  EXEMPLE:
%  var1=1; var2=false; var3='test';
%  S = aloha_setfield([], var1, var2, var3)
%  
%  S =
%      var1: 1
%      var2: 0
%      var3: 'test'
%  
% Since Matlab 2015b, there is a know incompatibility:  if the input 
% argument contains a dot indexing, the output will be an empty char array
% A workaround to this consists in wrapping the input argument with braces,
% that is : aloha_setfield([structure.field], varargin)
%    
%  AUTOR: JH
%  LAST UPDATE: 
%   - 02/09/2008 creation
%   - 20/10/2017 workaround for a know incompatibility since R2015b 
%  

    for ind=1:length(varargin)
        S = setfield(S, inputname(ind+1), varargin{ind});
    end
    %  on range les champs de la structure (pour faire plus propre)
    S = orderfields(S);