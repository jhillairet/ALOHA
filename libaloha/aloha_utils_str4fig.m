function str = aloha_utils_str4fig(str)
% Clean a string in order to be used into a figure
% 
% This functions clean some caracters in a string 
% in order to obtain a good display when this string
% is put into a title or a legend.
% 
% By example, this function replace underscore '_' by '\_', 
% so underscore are displayed as underscore in a figure, 
% not as a subscript.
% 
% EXAMPLE
%  str = aloha_utils_str4fig('antenne_C2')
%  ==> str = antenne\_C2
%  
% INPUT
%  - str [string] : original string
%  
% OUPUT
%  - str [string] : modified string
%   
% AUTHOR : JH
% LAST UPDATE:
%  - 22/10/2008: creation
  
    % replace '_' by an escape character '\_'
    str = strrep(str, '_', '\_');