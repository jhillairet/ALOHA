function scenario=aloha_compute_champEmbouchure2D(scenario, varargin)
% Compute the parallel component of electric field (Ez) in the  mouth of the antenna.
% Calcule le champ electrique parallele (Ez) dans l'embouchure de l'antenne
% 
% 2D case : the electric field correspond to the modes of a parallel-plate waveguide, ie 1 TEM and n TM modes.
%  

if nargin == 1
    modes_file = [pwd, '/ALOHA2D.out.modes.dat'];
elseif nargin == 2
    modes_file = varargin{1};
end

% for easier matlab manipulation, load all the fields of the input scenario 'scenario'
% into matlab workspace
aloha_scenario_loadIntoWorkspace;

% test if the main results, such as the plasma scattering matrix, 
% had been calculated. If not, warn user.
if not(exist('S_plasma'))
  error('The plasma scattering matrix has not been calculated !');
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

% open the modes ascii file given by the fortran binary for aloha-2D
fid = fopen(modes_file,'r');
    % read nb of mode and nb of guides
    C = textscan(fid, '%f %f', 1);
    nbre_guides= C{1};
    nbre_modes = C{2};
    C = textscan(fid, '%f', nbre_guides);
    a_grill = C{:};
    C = textscan(fid, '%f', nbre_guides);
    b_grill = C{:};
    C = textscan(fid, '%f', nbre_guides);
    y_grill = C{:};
    C = textscan(fid, '%f', nbre_guides);
    z_grill = C{:};
    C = textscan(fid, '%f', nbre_guides*nbre_modes);
    tab_ind_m = C{:};
    C = textscan(fid, '%f', nbre_guides*nbre_modes);
    tab_ind_n = C{:};
    C = textscan(fid, '%c', nbre_guides*nbre_modes);
    tab_TE_TM = C{:};
fclose(fid);



poids_E = sqrt(2)*(scenario.results.rac_Zhe)*(scenario.results.a_plasma + scenario.results.b_plasma);

z_guide=(0:0.04:1);
y_guide=(0:0.04:1);

% variable allocation
y = zeros(nbre_guides, length(y_guide));
z = zeros(nbre_guides, length(z_guide));
Ey = zeros(nbre_guides, length(y_guide), length(z_guide));
Ez = zeros(nbre_guides, length(y_guide), length(z_guide));


y_guide_tot=y_guide'*ones(1,length(z_guide));
z_guide_tot=ones(length(y_guide),1)*z_guide;

% for all waveguides
for ind_guide=1:nbre_guides
    
    a1=a_grill(ind_guide);
    b1=b_grill(ind_guide);
    y1=y_grill(ind_guide);
    z1=z_grill(ind_guide);
    
    
    Ey_guide=zeros(length(y_guide),length(z_guide));
    Ez_guide=zeros(length(y_guide),length(z_guide));
    
    for ind_modes=1:nbre_modes
        
        poids_mode=poids_E(ind_modes+(ind_guide-1)*nbre_modes);
        m1=tab_ind_m(ind_modes+(ind_guide-1)*nbre_modes);
        n1=tab_ind_n(ind_modes+(ind_guide-1)*nbre_modes);
        TE_TM_1=tab_TE_TM(ind_modes+(ind_guide-1)*nbre_modes);
        

    
        % coefficients de normalisation
        if (TE_TM_1=='H') 
      
            coeff_mode_1y=(n1/b1)/sqrt((m1^2)*b1/a1+(n1^2)*a1/b1);
            coeff_mode_1z=-(m1/a1)/sqrt((m1^2)*b1/a1+(n1^2)*a1/b1);
            if (m1>0) 
                coeff_mode_1y=coeff_mode_1y*sqrt(2);
                coeff_mode_1z=coeff_mode_1z*sqrt(2);
            end
            if (n1>0) 
                coeff_mode_1y=coeff_mode_1y*sqrt(2);
                coeff_mode_1z=coeff_mode_1z*sqrt(2);
            end
      
        else 
            
            coeff_mode_1y=-2.*(m1/a1)/sqrt((m1^2)*b1/a1+(n1^2)*a1/b1);
            coeff_mode_1z=-2.*(n1/b1)/sqrt((m1^2)*b1/a1+(n1^2)*a1/b1);
      
        end
      
        % expression des modes
        
        cos_m1=cos(m1*pi*y_guide_tot);
        cos_n1=cos(n1*pi*z_guide_tot);
        sin_m1=sin(m1*pi*y_guide_tot);
        sin_n1=sin(n1*pi*z_guide_tot);
     
        Ey_guide=Ey_guide+poids_mode*coeff_mode_1y*cos_m1.*sin_n1;
        Ez_guide=Ez_guide+poids_mode*coeff_mode_1z*sin_m1.*cos_n1;
      
    end
    
    z(ind_guide,:) = b1*z_guide+z1;
    y(ind_guide,:) = a1*y_guide+y1;
    Ey(ind_guide,:,:) = Ey_guide;
    Ez(ind_guide,:,:) = Ez_guide;

end

zz = [];
EE = [];
% calculate the norm of the E field in the middle of the waveguide
for idx_wg = 1:nbre_guides
    Ey_g = squeeze(Ey(idx_wg,:,:));
    Ez_g = squeeze(Ez(idx_wg,:,:));
    z_g = z(idx_wg,:);
    y_g = y(idx_wg,:);

    % look for the middle of the waveguide (in y direction)
    % i.e. where the E field should be maximum
    y_middle = (max(y_g)-min(y_g))/2;
    % take the closest vector index from geometric middle
    [dummy, idx_y_middle] = min(abs(y_g-y_middle)); 
    % take fields for this value
    Ey_g = Ey_g(idx_y_middle,:);
    Ez_g = Ez_g(idx_y_middle,:);
    norm_E = sqrt(abs(Ey_g).^2+abs(Ez_g).^2);
    
    zz = [zz, z_g];
    EE = [EE, norm_E];
end

abs_z = zz;
scenario.results = aloha_setfield(scenario.results, y, z, abs_z, EE, Ez_g, Ey, Ez, nbre_guides); 

