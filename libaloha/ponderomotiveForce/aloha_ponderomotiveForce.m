function [sc,varargout]=aloha_ponderomotiveForce(scenario, iter_nb, varargin)
% ALOHA ponderomotive force effect evaluation.
% 
% INPUT :
%  - scenario <struct (1x1)> : ALOHA scenario (without results). 
%                       This scenario should be unique (not a vector).
%                       If a vector is passed, only the first element is considered
%                       as the reference scenario.
%  - iter_nb <integer> : number of iteration to make
% 
% OPTIONAL INPUT ARGUMENTS :
%  - T : temperature = Te+Ti [eV] Default : 10+20=30 eV
% 
% OUTPUT :
%  - sc <struct (iter_nbx1)>: resulting ALOHA scenarios. 
%                           There is iter_nb scenarios returned, one for each iteration.
%  
%  AUTHOR : J.Hillairet
%  LAST UPDATE:
%  - 23/02/2010 : creation


eps0 = 8.85e-12; % permittivity of free space [F/m]
e = 1.602176487e-19; % elementary charge in [C]
T = 10+20; % defaut temperature in [eV] (T=Te+Ti)

% cut-off density 
nc = 0.0124*aloha_scenario_get(scenario, 'freq')^2;


% parsing optionnal input arguments
if nargin == 3
    T = varargin{1};
end

% temperature must be expressed in [J]
% so we convert [eV]->[J]
T = T*e;


% loop procedure : depending on the intensity of the Efield
% at the grill mouth, the electron density is changed

% initial values 
ne0 = aloha_scenario_get(scenario, 'ne0');
delta = zeros(1,iter_nb);
Ez_mean = zeros(1,iter_nb);
Ez_max = zeros(1,iter_nb);

% we duplicate the reference scenario iter_nb times
sc = repmat(scenario, 1, iter_nb);


h = waitbar(0,'Please wait...');
for iter=1:iter_nb
    % new density loop
    if iter > 1
        % get average electric field from (n-1) th scenario
        [Ez_mean(iter), Ez_max(iter)] = aloha_compute_averageEz(sc(iter-1)); 

        % calculate the new electron density as given by the
        % following ponderomotive expression
        delta(iter) = eps0 * abs(Ez_mean(iter)).^2 / (4*nc*T);
    end
    ne_iter = ne0 * exp(-delta(iter));

    % updating nth scenario electron density
    sc(iter).plasma.ne0 = ne_iter;

   % processing nth ALOHA scenario
    sc(iter) = aloha_scenario(sc(iter));

    waitbar(iter/(iter_nb),h);
end
close(h)

if nargout > 1
    varargout{1} = delta;
end
