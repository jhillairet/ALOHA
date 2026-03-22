function const = aloha_constants_get(constName, varargin)
%  Get some usefull plasma related constants
%  
% EXAMPLES
%  nc = aloha_constants_get('nc_LH', 3.7e9)
%  
% INPUT
%  - constName : constant's names. See below.
%  - optional inputs, which depends of the constant to be computed. See below.  
%  
% OUPUT
%  - const : the constant
% 
% Constant names            description                     required additional argument(s) 
%   'nc_LH'         LH electron cut-off density.                 generator frequency [Hz]
%  
%  
% AUTHOR(S): J.Hillairet
%     
% LAST UPDATE:
%  - 06/02/2009 : creation. nc_LH    

% if some physical constants, such as the electron mass, does not exist 
% in the matlab workspace, then we run aloha_constants, which is a script
% in which these physical constants are defined.
if ~exist('me')
    aloha_constants;
end

switch lower(constName)
    case 'nc_lh'
        if (nargin ~= 2)
            error('One more input is required : generator frequency in Hz.');
        end
        f = varargin{1}; 
        const = ((2*pi*f)^2)*me*Eps0/(qe^2);    

    otherwise
        error('unknow constant name!');
end