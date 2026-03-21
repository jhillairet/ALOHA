function ntpar = aloha_ntnpar(scenario)
%  ALOHA_NTNPAR(scenario)
% Adaptation de la fonction ntnpar(pulsenb,occur,com)
% de R. MASSET  (version du 28 octobre 2003)
%
% Cette fonction enregistre les resultats calcules par ALOHA 
% dans la base de donne TS, dans le cadre du traitement a la demande. 
% (traitement TNPAR)
% 
% 
% Cette fonction necessite les fonctions : 
%  - comment
% qui sont dans la toolbox IRFM. Ces fonctions doivent donc etre
% accessibles a matlab.

% ==================================
% creation des signaux a enregistrer
% ==================================
%  occur=0;com='';
%  tserv=clock;



% recuperation des informations concernant le choc
pulsenb = scenario.options.choc
t_start = scenario.options.tps_1;
t_end   = scenario.options.tps_2;
TSport  = scenario.options.TSport;

disp(aloha_message(['TNPAR treatment on TS pulse #' num2str(pulsenb)]));
disp(aloha_message(['between t=',num2str(t_start),' s and t=', num2str(t_end), 's.']));

% moyenne des temps d'acquisition : 
% sert de temps pour lequel le spectre 
% ainsi que les consignes de phases sont donnes.
tmean = (t_start+t_end)/2;


%Lancement des programmes de calcul gtnpar
%[tdata,phcons,phmes,nparmax,npar,spectre1,spectre2,cert,scom]=gtnpar(pulsenb); 






% Preparation du signal commentaire
% signal commentaire : SCHYBNPAR
commentaire=['Spectres N// des coupleurs Hybride' 10];
commentaire=[commentaire 'Indice 1 : coupleur C2 en Q6A' 10];
commentaire=[commentaire 'Indice 2 : coupleur C3 en Q6B'];


%  sctnpar=comment(scom);
%  ncar=length(sctnpar);

% Dephasages demandes
[phcons1,tcons1] = tsbase(pulsenb,'gtesthyb%3'); % Q6A
[phcons2,tcons2] = tsbase(pulsenb,'gtesthyb%4'); % Q6B
ph1cons = mean(phcons1(tcons1>t_start & tcons1<t_end));
ph2cons = mean(phcons2(tcons2>t_start & tcons2<t_end));

% Dephasages mesures
[phmes1,tmes1] = tsbase(pulsenb,'gphashyb%1'); % Q6A
[phmes2,tmes2] = tsbase(pulsenb,'gphashyb%2'); % Q6B 
ph1mes = mean(phmes1(tmes1>t_start & tmes1<t_end));
ph2mes = mean(phmes2(tmes2>t_start & tmes2<t_end));


% spectre et n//
dP    = squeeze(aloha_scenario_get(scenario, 'dP_nz'));
n_par = squeeze(aloha_scenario_get(scenario, 'nz'));
Nb_dP   = size(dP,1);
Nb_n_par= size(n_par,1);

% maximum du spectre : GHYBNPAR
[max_dP, idx_max_dP] = max(dP);
n_par0 = n_par(idx_max_dP);

% directivite ponderee en n//^2 
[dummy, direc] = aloha_compute_directivity1D(scenario, 2, n_par0);

% directivite nette, ponderee en n//^2
directivite_cumulee = scenario.results.directivite_cumulee;


% ==============================
% creation des donnees produites
% ==============================

GHYBPHCONS(1) = ph1cons;
GHYBPHCONS(2) = ph2cons;

GHYBPHMES(1) = ph1mes;
GHYBPHMES(2) = ph2mes;

GHYBNPAR(1) = n_par0;
GHYBNPAR(2) = n_par0; % TODO

GHYBEFFN(1) = direc;
GHYBEFFN(2) = direc; % TODO

GHYBSPEC1 = dP(:);
GHYBSPEC2 = dP(:); % TODO

GHYBCDIR1 = directivite_cumulee(:);
GHYBCDIR2 = directivite_cumulee(:); % TODO

SCHYBNPAR = comment(commentaire);

% ==============================
%  Taille des vecteurs
% ==============================
% Deuxieme coordonnee des GHYBSPEC1 et GHYBSPEC2 
nbt1=size(GHYBSPEC1,1);
t1=tmean; 

%  TODO : GHYBSPEC2
nbt2=size(GHYBSPEC2,1);
t2=tmean; 

% Deuxieme coordonnee des GHYBDIR1 et GHYBDIR2 
nbt3=size(GHYBCDIR1,1);
t3=tmean; 

% TODO : GHYBDIR2
nbt4=size(GHYBCDIR2,1);
t4=tmean;

% Lecture du temps ignitron
tigni=tsbase(pulsenb,'RIGNITRON');

% Conversion des temps en microsecondes
temps=(tmean+tigni)*1e6;


% =========================================================
%  Certifications
% =========================================================

cert1=-2;     %phcons1
cert2=-2;     %phcons2
cert3=-2;     %phmes1
cert4=-2;     %phmes2
cert5=-2;     %nparmax1
cert6=-2;     %nparmax2
cert7=-2;     %eff1
cert8=-2;     %eff2
cert9=-2;     %spec1
cert10=-2;    %spec2
cert11=-2;    %dir1
cert12=-2;    %dir2
cert13=-2;    %scom

% Certifications eclatees pour GHYBPHCONS, GHYBPHMES, GHYBNPAR, GHYBDIRN
certa=[101;cert1;cert2;101;cert3;cert4;101;cert5;cert6;101;cert7;cert8];
certb=[cert9;cert10;cert11;cert12;cert13];
cert=[certa;certb];


% ecriture de la structure de sortie
ntpar.cert = cert;   
ntpar.tdata = temps;
ntpar.npar = n_par;     

ntpar.phcons =  GHYBPHCONS;
ntpar.phmes = GHYBPHMES;   
ntpar.scom = SCHYBNPAR;
ntpar.spectre1 = GHYBSPEC1;
ntpar.spectre2 = GHYBSPEC2;
ntpar.direc1 = GHYBCDIR1;   
ntpar.direc2 = GHYBCDIR2; 
ntpar.effic = GHYBEFFN;  
ntpar.nparmax = GHYBNPAR;
  
fileName = ['ntpar_TS', num2str(pulsenb), '_', ...
            num2str(round(t_start)), '_', num2str(round(t_end)), '_', ...
            scenario.options.TSport, '_', scenario.antenna.architecture];

save(fileName, 'ntpar');

%  % Transposition des groupes pour l'ecriture
%  phc=GHYBPHCONS';
%  phm=GHYBPHMES';
%  par=GHYBNPAR';
%  eff=GHYBEFFN';
%  sp1=GHYBSPEC1';
%  sp2=GHYBSPEC2';
%  dir1=GHYBCDIR1';
%  dir2=GHYBCDIR2';
%  
%  nbt=size(temps,1);
%  npr=size(n_par,1);
%  ncar = size(SCHYBNPAR,1);
%  % Deuxieme coordonnee des groupes
%  nphc=size(phc,1);phc2=[1:nphc]';
%  nphm=size(phm,1);phm2=[1:nphm]';
%  npmx=size(par,1);npar2=[1:npmx]';
%  neff=size(eff,1);neff2=[1:neff]';
%  
%  
%  
%  %Structure du tableau a ecrire
%  % 1 groupe heterogene HTNPAR1
%  %  Signaux temporels
%  %  1 -    nbt  I4  : temps tdata
%  %  Groupes homogenes
%  %  2 - (nbt,ncoup)  R4  :  GHYBPHCONS
%  %  3 - (nbt,ncoup)  R4  :  GHYBPHMES
%  %  4 - (nbt,ncoup)  R4  :  GHYBNPAR
%  %  5 - (nbt,ncoup)  R4  :  GHYBEFFN
%  % 1 groupe heterogene HTNPAR2
%  %  5 -  npr       R4  :  N paralleles
%  % 4 Groupes de 2 dimensions
%  %  5 - (npr,nbt1)  R4  :  GHYBSPEC1
%  %  6 - (npr,nbt2)  R4  :  GHYBSPEC2
%  %  7 - (npr,nbt3)  R4  :  GHYBCDIR1
%  %  7 - (npr,nbt4)  R4  :  GHYBCDIR2
%  % Signal Commentaire
%  %  50 -  C4 : SCHYBNPAR
%  tab=[temps(:); ...
%       phc(:); phc2; ...
%       phm(:); phm2; ...
%       par(:); npar2; ...
%       eff(:); neff2; ...
%       n_par(:); ... % 
%       sp1(:); t1(:); ...
%       sp2(:); t2(:); ...
%       dir1(:); t3(:); ...
%       dir2(:); t4(:); ...
%       SCHYBNPAR];
%  
%  % Types donnees (temps I4,dontrait R4,commentaire C4)
%  typ=[1;2;3];
%  
%  % Longueurs des types donnees traitees
%  ltyp=[nbt;(nphc+nphm+npmx+neff)*(nbt+1)+npr+(npr+1)*(nbt1+nbt2+nbt3+nbt4);ncar];
%  
%  % Longueurs des coordonnees
%  lcoo=[nbt;nphc;nphm;npmx;neff;npr;nbt1;nbt2;nbt3;nbt4;ncar/4];
%  
%  % Controle de la structure a ecrire
%  test = isinf(tab) | isnan(tab);
%  if any(any(test))
%     disp('tab contient des infinis ou des NaN''s');
%     disp('Mise a zero de ces valeurs');
%     tab(test)=zero(tab(test));
%  end
%  
%  
%  
%  % =============================================
%  %   Ecriture dans la base
%  % =============================================
%  nom = 'TNPAR';
%  tserv = clock;
%  
%  %   [cr]=nbwrite(nom,choc,lcoo,cert,tab,typ,ltyp,[occur,comment])
%  %    Routine d'ecriture de toutes les donnee traitees
%  %     d'un traitement pour un choc donne.
%  %    R. MASSET le 9 Fevrier 1993.
%  %    Ajout du format caractere*4 le 7 Juillet 1993.
%  %    Entrees :
%  %     nom   : nom du traitement
%  %     choc  : numero de choc
%  %     lcoo  : longueurs des coordonnees   [ncoo,1]
%  %     cert  : certifications des donnees  [ncert,1]
%  %     tab   : tableau a ecrire            [ntab,1]
%  %     typ   : types des sequences a ecrire (1:I4, 2:R4, 3:C4)  [ntyp,1]
%  %     ltyp  : longueurs des sequences a ecrire (nbre de reel8) [ntyp,1]
%  %     occur : numero occurrence a ecrire  [1,1]
%  %     comment : commentaire a ecrire (32 caracteres maxi) [1,ncar]
%  %    Sorties :
%  %     cr    : compte rendu d'ecriture
%  occur = 0; com = '';
%  %  keyboard
%  %cr = nbwrite(nom,pulsenb,lcoo,cert,tab,typ,ltyp,occur,com);
%  
%  disp(aloha_message(['temps serveur (ecriture) : ' num2str(etime(clock,tserv))]))

