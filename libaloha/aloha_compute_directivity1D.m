function [scenario, varargout] = aloha_compute_directivity1D(scenario, varargin)
%  scenario = aloha_compute_directivity1D(scenario)
%  or 
%  [scenario, [D]] = aloha_compute_directivity1D(scenario [, definition_directivite] [, n_parallel_0])
%  
%  compute the directivity of the spectrum nz.
%  
%  WARNING: the spectrum must have been computed before in the scenario!
% 
%  The directivity calculation depends on the definition we use. 
%  The definition is set using the field (integer) 'definition_directivite' in scenario.option structure
%  or by setting an optional input argument.
%  
%  Definition are
%   
%  * definition_directivite=1 : 
%        D = (1/P)*int{nz=[+1,+inf]}(dP)
%   This is the most usual definition : ratio of the positive part of the spectrum 
%   over total power. This results may be considered in percents.
%  
%  * definition_directivite=2 : 
%       D = delta_cd = (1-R)*nz0^2/P*( int{nz=+1,+inf}(dP/nz^2) - int{nz=-inf,-1}(dP/nz^2))
%  This is the 'weighted directivity' as defined in [Litaudon & Moreau, Nucl Fusion 30 (1990), 471]
%  This directivity definition is based on the Fish current drive efficiency of the LH waves.
%  It is a quantitative measure, based on the CD efficiency which decrease as 1/n//^2, 
%  and of the quality of the launched spectrum with regard to an ideal Dirac spectrum, 
%  centered on n//0.
%  
%   * definition_directivite=3 : 
%       D = (1/P)*int{+1,+inf}(dP/nz^2);
%       
%   where P is the toral launched power : P = (int{+inf,+inf}(dP)
%   and R is the averaged reflection coefficient.
%  
%  INPUT:
%   scenario: aloha scenario
%   [optionnal] definition_directivite : integer
%   [optionnal] n_parallel_0 : assumed central peak of the spectrum. 
%               May be used for the computation of the weighted directivity.
%               If not set, assumes it is the max peak value of the spectrum.
%
%  OUTPUT
%   scenario: corresponds to the input scenario plus the directivity field inside.
%  [optionnal] : directivity
%  
% AUTHOR: J.Hillairet
% LAST UPDATE:
%  - 06/2009: put into a function.  
%  - 09/2009: cleanning the different directivity definition 
%             and adding some optionnal input/output arguments

% if many scenarios
for idx_sc = 1:length(scenario)
    sc = scenario(idx_sc);

    % input arguments section
    if nargin == 1
        definition_directivite = aloha_scenario_get(sc, 'definition_directivite');
    elseif nargin == 2
        definition_directivite = varargin{1};
    elseif nargin == 3
        definition_directivite = varargin{1};
        n_parallel_0 = varargin{2};
    end

    % input tests     
    if not(isfield(sc.results, 'dP_nz'))
        error('The spectrum must have been computed before ! see help aloha_compute_spectrum');
    end
    
    % retrieve data from scenario 
    % parallel index
    nz = aloha_scenario_get(sc, 'nz');
    dnz = aloha_scenario_get(sc, 'dnz');
    % power density in the nz space
    dP_nz = aloha_scenario_get(sc, 'dP_nz');
    % total power launched
    P = sum(dP_nz)*dnz; 
    % reflected power coefficient
    RC = aloha_scenario_get(sc, 'CoeffRefPuiss');
    
    D = 0;
    D_cumulated = zeros(1,length(nz));

    switch definition_directivite
        case 1
            D = sum(dP_nz(find(nz>1)))*dnz/P ;    

            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                D_cumulated(idx) = 1-dnz*sum(dP_nz([idx:end]))/P;
            end


        case 2   
            % In this calculation, it is assumed that the spectrum is centred on n_parallel_0.
            % However, if this variable had not been set to the present function,
            % we assume that this value corresponds to the nz for which the spectrum is max.
            if not(exist('n_parallel_0'))
                [dummy, idx_max] = max(abs(dP_nz));
                n_parallel_0 = nz(idx_max);
            end 
        
            if n_parallel_0 > 0
                D = (1-mean(RC,2)/100).*n_parallel_0.^2 * (1/P) .*  ...
                    (dnz*sum(dP_nz(find(nz>1))./nz(find(nz>1)).^2) ...
                    - dnz*sum(dP_nz(find(nz<1))./nz(find(nz<1)).^2));
            else
                D = (1-mean(RC,2)/100).*n_parallel_0.^2 * (1/P) .*  ...
                    (dnz*sum(dP_nz(find(nz<1))./nz(find(nz<1)).^2) ...
                    - dnz*sum(dP_nz(find(nz>1))./nz(find(nz>1)).^2));
            end    

            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                 D_cumulated(idx) = 1 - (1/P)*dnz*(sum(dP_nz([idx:end])./nz([idx:end]).^2));
            end

        case 3
            D = (1/P)*dnz*(sum(dP_nz(find(nz>0))./nz(find(nz>0)).^2));

            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                D_cumulated(idx) = 1 - (1/P)*dnz*(sum(dP_nz([idx:end])./nz([idx:end]).^2));
            end

    
        otherwise
            error('Bad or non existent directivity definition !') ;
    end
    
    disp(aloha_message(...
    ['Directivity (definition #', num2str(definition_directivite),'): D=', num2str(D)]));
    
    % save results into the output scenario
    sc.results.directivite = D; 
    sc.results.directivite_cumulee = D_cumulated;
    
    directivite(idx_sc) = D;
    scenario(idx_sc) = sc;
end

% output variables
if nargout == 2
    varargout{1} = directivite;
end 