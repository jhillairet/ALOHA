function scenario=aloha_compute_champEmbouchure1D_v2(scenario)
% Compute the parallel component of electric field (Ez) in the  mouth of the antenna.
% Calcule le champ electrique parallele (Ez) dans l'embouchure de l'antenne
% 
% scenario=aloha_compute_champEmbouchure1D_v2(scenario)
% 
% 1D case : the electric field correspond to the modes of a parallel-plate waveguide, ie 1 TEM and n TM modes.
% However, the E & H fields amplitudes are normalized to the waveguide largest side (a), 
% ie. the amplitudes would be similar in a full 3D modeling. 
% 
% If one wants to compare with a 2D modeling, he has to activate the 'SWAN' mode in ALOHA.
% 
% INPUT
%   - scenario : ALOHA scenario (coupling allready calculated)
% OUPUT
%   - scenario : ALOHA scenario with the additionnal fields in the 'results' part:
%           - abs_z(Nb_poloidal_lines, Nb_z)  : parallel abcisse z
%           - Efield(Nb_poloidal_lines, Nb_z) : parallel electric field (Ez)
%           - E_mouth(3, Nb_z, Nb_poloidal_lines) : Electric field [Ex,Ey,Ez]
%           - H_mouth(3, Nb_z, Nb_poloidal_lines) : Magnetic field [Hx,Hy,Hz]
%           - Poynting_mouth(3, Nb_z, Nb_poloidal_lines) : Poynting vector [Px,Py,Pz]
%           - Poynting_mouth_total : total power flux integrated along the toroidal direction
%  
% AUTHOR: JH
% LAST UPDATE:
%  - 01/12/2012: corrected a bug in the z sampling, which was  

% for easier matlab manipulation, load all the fields of the input scenario 'scenario'
% into matlab workspace
aloha_scenario_loadIntoWorkspace;

% load physical constants
aloha_constants;

k0 = 2*pi*freq/celerite;

% test if the main results, such as the plasma scattering matrix, 
% had been calculated. If not, warn user.
if not(exist('S_plasma'))
  error('The plasma scattering matrix has not been calculated ! Please run aloha_scenario before.');
end

% Load into workspace the geometrical parameter of the scenario architecture 
% and get the main geometrical parameters.
if aloha_isAntennaITM(scenario)
    disp(aloha_message('assuming ITM antenna description')); 
    aloha_utils_ITM2oldAntenna; % convert ITM data to old fashioned ALOHA global parameters
    [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinatesFromCPO(aloha_getAntenna(scenario));
else
    eval(architecture);
    [b,h,z,y,nb_g_total_ligne,nbre_guides,act_module_tor]=aloha_utils_getAntennaCoordinates(architecture);
end

% Modification JH 04/2009
% Precedemment, l'echantillonnage n'était pas régulier : les champs dans les guides etaient 
% decrit sur 100 points, tandis que le champs au niveau des septums etaient decrit avec 1 seul point.
% Lorsqu'on faisait la transformee de fourier du champ total, on obtenait alors un resultat "suprenant".
% 
% La modification du code consiste donc a exprimer les champs avec un pas constant et controle.
% 
% ATTENTION : doit etre absolument inferieur a la plus petite decimale des dimensions b ou e 
% pour resoudre correctement spatialement en z.
% Ex for the C4 launcher: dz=1e-6 since e=4.325mm
dz = 1e-5; 

% The previous comment means that the total number of point should be an integer. 
% ie: dz should be sufficiently small !
%  number_z_points = (b(end)+z(end))/dz+1; 


% Clear all output fields variables
% PS: This is recommanded here since these variables can be overloaded from the 
% input scenario.
Efield  = [];
E_mouth = [];
H_mouth = [];
Poynting_mouth = [];
abs_z = [];

% The E and H field amplitudes are normalized to the waveguide largest side. 
% In order to have the same unit than Efield calculated with the 2D ALOHA mode (or any 3D FEM code).
% If one would compare this Efield with a 2D FEM code, then he has to remove the /a factor below. 
% This is automatically done when one activates the "SWAN" option in ALOHA.
%  
% Bug correction JH 13/11/2011
% The magnetic field amplitude was lower than expected. This has an impact on
% the amplitude of the spectrum when calculated from Electric & Magnetic field by FFT.
% It seems that the ratio 1/Z0=Y0 has been taken into account two times, so one time too much !
% Below I've added the inv(rac_Zhe) term to the 'poids_H' term (H field coupling weight), 
% and removed the YO term in the H-field calculation.
% Moreover, I've added the if bloc to take into account the 2D/3D Efield normalization 
% as explained before.
if type_swan_aloha == 1 % "a" normalization (rectangular waveguide, 3D)
    poids_E = sqrt(2/a)*sqrt(2)*rac_Zhe*(a_plasma + b_plasma);
    poids_H = sqrt(2/a)*sqrt(2)*inv(rac_Zhe)*(a_plasma - b_plasma);
elseif type_swan_aloha == 0 % no "a" dependance (parallel-plate waveguide, 2D)
    poids_E = sqrt(2)*rac_Zhe*(a_plasma + b_plasma);
    poids_H = sqrt(2)*inv(rac_Zhe)*(a_plasma - b_plasma);
end



disp(aloha_message('Calculating Electric and Magnetic fields at the mouth of the launcher')); 
% for all poloidal lines
for idx_pol = 1:nb_g_pol 
    abs_z_idx_pol = [];
    Etot=[];
    Htot=[];
    % for all waveguides in a poloidal line
    for idx_tor = 1:nb_g_total_ligne

        % local waveguide abscisse 
        z_g = [0:dz:b(idx_tor)];

        % ALOHA mode index
        mode_idx=1+(idx_tor-1)*(Nme+Nmh)+(idx_pol-1)*nb_g_total_ligne*(Nme+Nmh):Nme+Nmh+(idx_tor-1)*(Nme+Nmh)+(idx_pol-1)*nb_g_total_ligne*(Nme+Nmh);

        E_wg_m = [];
        H_wg_m = [];
        coeffx_wg = [];
        coeffy_wg = [];
        coeffz_wg = [];
        S = [];
        % compute fields for all modes
        for m = [0:Nmh-1+Nme]; % mode index

            %% electric field
            coeffx_wg = j/Y0*m*pi./k0./b(idx_tor).*Phi_m_prime(m, b(idx_tor), z_g);
            coeffy_wg = zeros(size(z_g));
            coeffz_wg = Phi_m(m, b(idx_tor), z_g);

            E_wg_m(:, :, m+1) = [coeffx_wg;coeffy_wg;coeffz_wg]*poids_E(mode_idx(m+1));

            %% magnetic field
            coeffx_wg = zeros(size(z_g));
            coeffy_wg = -Phi_m(m, b(idx_tor), z_g);
            coeffz_wg = zeros(size(z_g));
            % guided wavenumber 
            km = (k0^2>(m*pi/b(idx_tor)).^2).*sqrt(k0^2-(m*pi/b(idx_tor)).^2) ...
                -j*(k0^2<(m*pi/b(idx_tor)).^2).*sqrt((m*pi/b(idx_tor)).^2 - k0^2);
            % modal normalization
            Ym = k0/km; % = 1 if m=0 (TEM)

            H_wg_m(:, :, m+1) = Ym*[coeffx_wg;coeffy_wg;coeffz_wg]*poids_H(mode_idx(m+1));
        end
        % the field in a waveguide is the sum of all modes
        E_wg = sum(E_wg_m, 3);
        H_wg = sum(H_wg_m, 3);

        abs_z_idx_pol = [abs_z_idx_pol, z(idx_tor) + z_g];

        % Set fields to zero on septa
        if (idx_tor < nb_g_total_ligne) % sauf pour le dernier
          z_e = [z(idx_tor)+b(idx_tor)+dz:dz:z(idx_tor+1)-dz];
          abs_z_idx_pol = [abs_z_idx_pol, z_e];
          E_wg = [E_wg, zeros(size(E_wg,1),length(z_e))];
          H_wg = [H_wg, zeros(size(E_wg,1),length(z_e))];  
        end

        % AJOUT JH 10/2008 in order to save the results
        if (idx_tor == 1)
            E_idx_pol = E_wg;
            H_idx_pol = H_wg;
        else
            E_idx_pol = [E_idx_pol , E_wg];
            H_idx_pol = [H_idx_pol, H_wg];
        end

    end % idx_tor

   % variables pour sortie & sauvegarde
   abs_z(idx_pol,:) = abs_z_idx_pol;
   Efield(idx_pol,:)= E_idx_pol(3,:); % z

   E_mouth(:,:,idx_pol) = E_idx_pol;
   H_mouth(:,:,idx_pol) = H_idx_pol; 

   if type_swan_aloha == 1 % "a" normalization (rectangular waveguide, 3D)
    Poynting_mouth(:,:,idx_pol) = (a/2)*1/2*real(cross(E_idx_pol, conj(H_idx_pol)));
    % JH 12/09/2011
    % I have to admit that the 1/2 in the (a/2) term before is mainly ad-hoc, found by comparizon with
    % the power balance. The a term is logical in order to transform the power density (W/m) into
    % power (W), but I have no clue for the additional /2 factor (which does not comes from 
    % the RMS definition, since there is already a 1/2 factor before the Poynting vector...).  
   elseif type_swan_aloha == 0 % no "a" dependance (parallel-plate waveguide, 2D)
    Poynting_mouth(:,:,idx_pol) = 1/2*real(cross(E_idx_pol, conj(H_idx_pol)));
   end
   Poynting_mouth_total = trapz(abs_z(idx_pol,:), Poynting_mouth(1,:,idx_pol));   

   disp(aloha_message(['Summation of the poynting vector in x direction : ', num2str(Poynting_mouth_total), ' W']))

end % idx_pol



% add the calculated results to the scenario
scenario.results.abs_z=abs_z; 
scenario.results.dz=dz;
scenario.results.Efield=Efield; % compatibility output paramerters with previous ALOHA version
scenario.results.E_mouth=E_mouth;
scenario.results.H_mouth=H_mouth;
scenario.results.Poynting_mouth=Poynting_mouth;

% modal function 
function Phi_m = Phi_m(m,b,z)
    % eps = 1 for fundamental TEM mode m=0
    % eps = 2 for sup. TM modes m>=1
    eps = (m==0) + 2*(m>=1);

    Phi_m = sqrt(eps./b).*cos(m.*pi.*z./b);

% z derivative of the modal function
function Phi_m = Phi_m_prime(m,b,z)
    if (m == 0) % fundamental TEM mode
        Phi_m = zeros(size(z));
    elseif (m >= 1) % sup. TM modes
        Phi_m = -m.*pi./b.*sqrt(2./b).*sin(m.*pi.*z./b);
    else
        error('m index must be >= 0');
    end
