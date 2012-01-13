% Calcule de la matrice S antenne
%  
%  
%  NB JH : A priori, il n'est plus utile de se deplacer 
%  dans le repertoire contenant les matrices S, 
%  puisque tous les fichiers contenus ce repertoire et ses sous repertoires
%  ont ete inclus dans le PATH MATLAB lors de l'initialisation d'ALOHA

S_ant_11 = [];
S_ant_12 = zeros((nb_g_pol/nb_g_module_pol)*nb_modules_tor, length(S_plasma));
S_ant_21 = zeros(length(S_plasma), (nb_g_pol/nb_g_module_pol)*nb_modules_tor);

% on introduit les voies passives

S_ant_22 = zeros(1,length(S_plasma));

% MODIF JH 29/10/2008 : ajout de la condition ~isempty(pass_tot) qui permet 
% de modeliser des antennes sans voies passives.
if ~isempty(pass_tot)
    if (antenne_standard == 1) 
        S_ant_22(1,pass_tot) = -exp(+i*4*pi*lcc);
    else
        lcc=kron(ones(1,nb_g_pol),lcc); % replicates for the number of poloidal rows
        S_ant_22(1,pass_tot) = -exp(+i*4*pi*lcc);
%          % passive waveguide scattering parameter for a non ideal depth of lambdag/4.
%          S_ant_22(1,pass_tot) = 9.910618E-001 - i*1.334034E-001;
    end
end
S_ant_22 = diag(S_ant_22);

% extractions des matrices S_modules des fichiers HFSS
% 
% Pour tous les modules situes sur une ligne poloidale de modules
for ind = 1:(nb_g_pol/nb_g_module_pol)*nb_modules_tor
    % introduit les variables S [, Z et f] lues dans fichier issu de HFSS
    nom = nom_fichiers(ind,:);
    
    [dummy,dummy,ext]=fileparts(nom);
    if strcmp(ext, '') || strcmp(ext, '.m')
      % if the name has the .m extension, load as a script
      eval(nom);
    
    elseif strcmp(ext, '.mat')
      % otherwise if it has a .mat extension, load it a matlab matrix 
      load(nom);
      
    elseif regexp(ext, '.s.p')
      % touchstone file
      % Problem is : how many header lines does the file have ? Depends on the RF software which creates the file....
      %[S,f]=aloha_touchstone_read(nom, 1); % sometime doesn't work properly
      [freq, S, freq_noise, data_noise, Zo] = SXPParse(nom);
      S=reshape(S, 1, length(S)*length(S));
    else
        error('Bad or unknow file name/extension');
    end
    
    % J.Hillairet 14/12/07
    % Le format des nouvelles matrices calculees par J.Belo le 13/12/07 est different des anciennes : 
    % (1,7,7) au lieu de (1,49).
    % On procede alors a une petite modification pour assurer la compatibilite avec le reste du code :
    if ndims(S) ==  3
        S_module = squeeze(S);
    else
    % Matrice S format precedent
        S_module = reshape(S,sqrt(length(S)),sqrt(length(S)));
    end
    
    % Temporary piece of code for the conceptual study of the alcator CMod LH3 
    %
    % If the bijLength field exists in the scenario structure, 
    % then we proceed to a phase desembedding of the scattering matrix of the launcher
    % in order to change the electric length of the secondary waveguide after the bijunction.
    if isfield(scenario.antenna_lh.setup.modules.Sparameters, 'bijLength')
      beta = sqrt(k0^2 - (pi/60e-3)^2);
      bijLength = scenario.antenna_lh.setup.modules.Sparameters.bijLength;
      E = diag([1,repmat(exp(-j*beta*bijLength),1,8)]);
      S_module = E*S_module*E;
    end
    
    S_module_11 = S_module(1,1); 
    S_module_12 = S_module(1,:); 
    S_module_12(:,1) = [];
    S_module_21 = S_module(:,1);
    S_module_21(1,:) = [];
    S_module_22 = S_module;
    S_module_22(1,:) = [];
    S_module_22(:,1) = [];

    %  D.Voyer 05/06/2008 :
    %  Lorsqu'on utilise les données issues de la mesure
    %  Les references des phases doivent etre corrigees
    %  pour prendre en compte le dephasage du a la longueur
    %  des guides entre les fenetres HF et 
    %  l'entree du diviseur de puissance. 
    if (bool_mesure)   
      disp('[ALOHA] (INFO) : correction de la phase pour chaque module entre entree fenetres et disp. div. puissance');

        % if ITM antenna description 
        if aloha_isAntennaITM(scenario)
            phase_rallonge = scenario.antenna_lh.setup.modules.Sparameters.phase_deembedded;
        else % old antenna description
            % phase_rallonge has allready been defined in the antenna description
        end

        S_module_11 = S_module_11*exp(+i*2*phase_rallonge(ind)); 
        S_module_12 = S_module_12*exp(+i*phase_rallonge(ind)); 
        S_module_21 = S_module_21*exp(+i*phase_rallonge(ind));
    end
   
    S_ant_11 = [S_ant_11,S_module_11];
    S_ant_12(ind,modules_act(ind,:)) = S_module_12;
    S_ant_21(modules_act(ind,:),ind) = S_module_21;
    S_ant_22(modules_act(ind,:),modules_act(ind,:)) = S_module_22;
   
    
end

S_ant_11 = diag(S_ant_11);

%  cd(chemin_retour)