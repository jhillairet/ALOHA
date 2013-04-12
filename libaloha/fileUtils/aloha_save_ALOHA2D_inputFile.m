function aloha_save_ALOHA2D_inputFile(filename, varargin)
% aloha_save_ALOHA2D_inputFile(param1_name, param1, param2_name, param2, ...)
% 
% Save the input parameters list into a text file like :
%  param1_name = param1
%  param2_name = param2
% etc...
%
% INPUT
%  the number of input argument should be even. For each parameter, give
%  its name then its value
%
% OUPUT : none
%
% Author: JH
% Last update: 27/03/2013.

% check if the input argument number if even (should be!)
if (mod(length(varargin),2) ~= 0)
    error('The total number of input arguments should be even !' )
end

fid=fopen(filename, 'w');

% then process input arguments 2 by 2
for idx=1:length(varargin)/2
    arg_name{idx} = varargin{2*(idx-1)+1};
    arg_value{idx}= varargin{2*(idx-1)+2};
    
    arg_string{idx} = [arg_name{idx}, ' = ', num2str(reshape(arg_value{idx}, 1, []))];
    disp(arg_string{idx})
    fprintf(fid, '%s\n', arg_string{idx});
end

fclose(fid);