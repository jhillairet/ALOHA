function scenario=aloha_compute_spectrum1D(scenario)
%  Compute the power spectrum of the antenna
%  
%  INPUT
%   - ALOHA scenario
%  OUPUT
%   - ALOHA scenario with suppl. fields into the sub-field 'results"
%     corresponding to the spectrum parameters (dP, dP_nz, nz, ny)
%  
%  NB : All the options such as : 
%    dny, dnz, ny_min, ny_max, nz_min, nz_max, nbre_ny, nbre_nz
%  should be defined into the 'option' sub-field of the input scenario
%  
%  LAST UPDATE:
%   01/10/2009 : integration of aloha 1D and aloha 2D

switch scenario.options.version_code
    case '1D'
        scenario = aloha_compute_spectrum1D(scenario);
    case '2D'
        scenario = aloha_compute_spectrum2D(scenario);
    otherwise
        error('Wrong or no code_version parameter into the scenario !')
end
