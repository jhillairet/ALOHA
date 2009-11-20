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
%  AUTOR: JH
%  LAST UPDATE: 
%   - 02/09/2008 creation
%  

    for ind=1:length(varargin)
        S = setfield(S, inputname(ind+1), varargin{ind});
    end
    %  on range les champs de la structure (pour faire plus propre)
    S = orderfields(S);