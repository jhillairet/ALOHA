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
    P = sum(real(dP_nz))*dnz; 
    % reflected power coefficient
    RC = aloha_scenario_get(sc, 'CoeffRefPuiss');
    
    D = 0;
    D_cumulated = zeros(1,length(nz));

    % In the following calculations, it is assumed that the main peak of the spectrum is 'n_parallel_0'
    % However, if this variable had not been provided by the user in the input arguments to the present function,
    % we assume that this value corresponds to the nz for which the spectrum is maximum.
    if not(exist('n_parallel_0'))
        [dummy, idx_max] = max(real(dP_nz));
        n_parallel_0 = nz(idx_max);
    end 
    disp(aloha_message(['The maximum peak in the spectrum is located at n//= ',num2str(n_parallel_0)]))

    switch definition_directivite
        case 1
            % Positive n//0
            if n_parallel_0 > 0
                D = sum(real(dP_nz(find(nz>1))))*dnz/P ;    
            else
                D = sum(real(dP_nz(find(nz<1))))*dnz/P;
            end
            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                D_cumulated(idx) = 1-dnz*sum(real(dP_nz([idx:end])))/P;
            end


        case 2   
            % Positive n//0
            if n_parallel_0 > 0
                D = (1-mean(RC,2)/100).*n_parallel_0.^2 * (1/P) .*  ...
                    (dnz*sum(real(dP_nz(find(nz>1)))./nz(find(nz>1)).^2) ...
                    - dnz*sum(real(dP_nz(find(nz<1)))./nz(find(nz<1)).^2));
            else % negative n//0
                D = (1-mean(RC,2)/100).*n_parallel_0.^2 * (1/P) .*  ...
                    (dnz*sum(real(dP_nz(find(nz<1)))./nz(find(nz<1)).^2) ...
                    - dnz*sum(real(dP_nz(find(nz>1)))./nz(find(nz>1)).^2));
            end

            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                 D_cumulated(idx) = 1 - (1/P)*dnz*(sum(real(dP_nz([idx:end]))./nz([idx:end]).^2));
            end

        case 3
            D = (1/P)*dnz*(sum(real(dP_nz(find(nz>0)))./nz(find(nz>0)).^2));

            % compute the cumulated directivity, function of nz
            for idx=1:length(nz)
                D_cumulated(idx) = 1 - (1/P)*dnz*(sum(real(dP_nz([idx:end]))./nz([idx:end]).^2));
            end

    
        otherwise
            error('Bad or non existent directivity definition !') ;
    end
    
    disp(aloha_message(...
    ['Directivity (definition #', num2str(definition_directivite),'): D= ', num2str(D)]));
    
    %% calculate the relative height in % of the secondary lobe
    % 
    [dummy, idx_0] = min(abs(nz - n_parallel_0));
    if n_parallel_0 > 0
        % Find the location of the maximum peak in the negative region : this is the opposite lobe
        [dummy, idx_opp] = max(real(dP_nz.*(nz<1)));
        n_parallel_opposite = nz(idx_opp);
    else
        % Find the location of the maximum peak in the positive region : this is the opposite lobe
        [dummy, idx_opp] = max(real(dP_nz.*(nz>1)));
        n_parallel_opposite = nz(idx_opp);
    end
    % Opposite peak relative height  to the main peak 
    opposite_peak_relative_height = real(dP_nz(idx_opp))/real(dP_nz(idx_0))*100;

    disp(aloha_message(['Opposite peak located at n//= ',num2str(n_parallel_opposite)]));
    disp(aloha_message(['Opposite peak relative height to the main peak: ', num2str(opposite_peak_relative_height), '%']));

    %% save results into the output scenario
    sc.results.opposite_peak_location = n_parallel_opposite;
    sc.results.opposite_peak_relative_height = opposite_peak_relative_height;
    sc.results.main_peak_location = nz(idx_0);
    sc.results.directivite = D; 
    sc.results.directivite_cumulee = D_cumulated;
    
    directivite(idx_sc) = D;
    scenario(idx_sc) = sc;
end

% output variables
if nargout == 2
    varargout{1} = directivite;
end 