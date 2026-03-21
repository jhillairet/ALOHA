function [data, time] = aloha_experimental_get(machine, pulse, signal, varargin)
% Get experimental datas from a machine database
% 
% 
% EXAMPLE
%  [RC, t_RC] = aloha_experimental_get('TS', 
% 
% INPUT
%  - machine <string> : tokamak name : 'TS', ...
%  - pulse <integer>: pulse number
%  - signal <string>: signal name
%  - t1 <real> [optionnal] : start time
%  - t2 <real> [optionnal] : stop time
% 
% OUPUT
%  - data : data vector
%  - time : time vector
%               
% AUTHOR:  J.Hillairet
% LAST UPDATE:
%  - 06/02/2009 : work for Tore Supra
%  

switch lower(machine)

    case {'ts', 'toresupra'}        
        [data, time] = tsbase(pulse, signal);
        % keep only one time vector
        time=time(:,1);
        
    otherwise
        error('undefinited machine!');
end

% if the time interval exists, we crop the signal to [t1,t2]
if (nargin > 3)
    t1 = varargin{1};
    t2 = varargin{2};

    % cut signal to [t1,t2]
    id_sup = find(time>t1);
    id_inf = find(time<t2);
    id = intersect(id_sup, id_inf);

    time = time(id);
    data = data(id,:);
end
