function [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(antenna_lh)
% Extract the pertinent information for ALOHA from the ITM CPO 'antenna_lh'
% 
% 



%% Shortcuts
mod= antenna_lh.setup.modules;
wg = antenna_lh.setup.modules.waveguides;

% (total) number of waveguides per row
% = (nb wg in a module) + 2*(nb ext wg) + (nb of wg between modules)
nwr = mod.nma_phi * wg.nwm_phi  + 2*wg.npwe_phi + (mod.nma_phi-1)*wg.npwbm_phi;

% (total) number of waveguides per column
nwc = mod.nma_theta*wg.nwm_theta;

% total number of waveguides
% JH 04/12/2013
% BUG discovered when mutli poloidal rows... to be cleaned... one day...
% maybe
nwa = nwr*nwc;% was nwa = nwr*mod.nma_theta;

% waveguide height - supposed constant for all the waveguides of the antenna
a = wg.hw_theta;

 

%% b
% Make the array b which contains all the waveguide width 
% of a row of waveguides
if not(isfield(wg, 'b'))
    b_module = wg.mask.*wg.bwa + not(wg.mask).*wg.biwp; % waveguide width inside a module
    b_edge = repmat(wg.bewp, 1, wg.npwe_phi);  % passive wg width on each side
    b_inter= repmat(wg.biwp, 1, wg.npwbm_phi); % passive wg width between modules

    b = [b_edge, kron(ones(1,mod.nma_phi-1),[b_module, b_inter]),b_module, b_edge];
else
    b = wg.b;
end
%% e
% Make the array e which contains all the waveguide septum width
% of a row of waveguides 
e = wg.e_phi;

%% z
% Make the array z which contains all the waveguide positions 
% in the toroidal direction
z = zeros(1,nwr);
for ind = 2:nwr
    z(ind) = z(ind-1) + b(ind-1) + e(ind-1);
end

%% y
% Make the array y which contains all the waveguide positions
% in the poloidal direction
h = [kron(ones(1,mod.nma_theta-1), [repmat(wg.sw_theta,1,wg.nwm_theta),mod.sm_theta]), repmat(wg.sw_theta,1,wg.nwm_theta)]; 
y = zeros(1, nwc);
for ind = 2:nwc
    y(ind) = y(ind-1) +  h(ind-1) + a;
end

% output variables
nb_g_total_ligne = nwr;
nbre_guides = nwa;
% index of active waveguides in a module
act_module_tor = find(wg.mask == 1);
