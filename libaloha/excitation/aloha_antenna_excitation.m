function [a_ampl, a_phase] = aloha_antenna_excitation(varargin) 
%  Retrieve amplitude and phase exciatation at the input of each module of an LH antenna
%  
%  EXAMPLE
%  2 calls are possible :
%   retrieve predefinite configurations :
%   [a_ampl, a_phase] = aloha_antenna_excitation(architecture) 
%  Or,
%   retrieve measured data from the TS database:
%   [a_ampl, a_phase] = aloha_antenna_excitation(choc, tps_1, tps_2, TSport) 
% 
% PREDEFINTE EXCITATION
% Predefinite excitation are, for 25 MW/m^2 : 
%   Pin = { 2.0672 , 4.0320 , 2.6725 } MW for { C2 , C3 , C4 }
%   (fields amplitudes are deduced sqrt(Pin/number_of_modules/number_of_power_divider), eq:/8/2 on TS)
%   Phasing = { 0 , -90 , 180 } MW for { C2 , C3 , C4 } for nominal design values, which are {1.83, 2.03, xxx }
%   (C3 : -150 pour n//=1.87 as C2)
%  
% INPUTS
%  - architecture [string] : nom de l'architecture pre-definies
%  Or,
%  - choc [int]  : numero du choc TS
%  - tps_1 [real]: temps lecture debut
%  - tps_2 [real]: temps lecture fin
%  - TSport [string] : % 'Q6A' [C2 before 2009, C3 after] ou 'Q6B' [C3 before 2009, C4 after]
%  
% OUPUT 
%  - a_ampl
%  - a_phase
% 
% AUTHOR: OI/DV/JH
% LAST CHANGE: 
%  - 09/09/2008: creation
%   


%  
% 1 argument : 
% Excitations predefinies
% 
if (nargin == 1)
    architecture = varargin{1};

    % Pour 25 MW/me mettre Pi = { 2.0672 , 4.0320 , 2.6725 } MW pour { C2 , C3 , C4 }
    % Pour 20MW pour l'antenne d'ITER mettre 5 MW pour 1/4 d'antenne
    % Dephasage = { 0 , -90 , 180 } MW pour { C2 , C3 , C4 } pour respecter la periodicite 
    % (C3 : -150 pour n//=1.87 comme C2)
    switch architecture
        case 'antenne_C2'
            a_ampl = sqrt(2.0672e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(0*pi/180)*(0:7)';

        case 'antenne_COMPASS_old2'
            a_ampl = sqrt(1)*ones(8,1);   
            a_phase = 1*(-60*pi/180)*(0:7)';
            
         case 'antenne_COMPASS_old2_TEST'
            a_ampl = sqrt(2.0672e6/16)*ones(20,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(60*pi/180)*(0:19)';   
            
        case 'antenne_COMPASS_old3'
            a_ampl = sqrt(2.0672e6/16)*ones(6,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(60*pi/180)*(0:5)';
        
        case 'antenne_COMPASS_6wg_HAMAC'
            a_ampl = sqrt(2.0672e6/16)*ones(2,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(0*pi/180)*(0:1)';
        
        case 'antenne_COMPASS_8wg_HAMAC'
            a_ampl = sqrt(2.0672e6/16)*ones(2,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(0*pi/180)*(0:1)';   
         
        case 'antenne_COMSOL_16WG_HAMAC'
            a_ampl = sqrt(2.0672e6/16)*ones(2,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(0*pi/180)*(0:1)';   
          
        case {'antenne_C3', 'antenne_C3_referenceModule', 'antenna_C3_ITM'}
            a_ampl = sqrt(4.0320e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = -90*(pi/180)*(0:7)';

        case {'antenna_C3_WEST_ITM'}
            a_ampl = sqrt(4.0320e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = +90*(pi/180)*(0:7)';

        case 'antenne_C3_COMSOL'
            a_ampl = 1*ones(16,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(-90*pi/180)*(0:15)';
  
        case 'antenne_C3_Coupled_COMSOL'
            a_ampl = sqrt(1e5)*ones(8,1);   % En champ elect 
            a_phase = 1*(-90*pi/180)*(0:7)';
            
        case 'antenne_PAM_Coupled_COMSOL'
            a_ampl = sqrt(1e5)*ones(8,1);   % En champ elect 
            a_phase = 1*(270*pi/180)*(0:7)';
 
        case 'antenne_C3_16WG_Coupled_COMSOL'
            a_ampl = sqrt(4e5)*ones(2,1);   % En champ elect 
            a_phase = 1*(-90*pi/180)*(0:1)';
 
        case 'antenne_C3_COMSOL_SupportCase'
            a_ampl = 1*ones(4,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(-90*pi/180)*(0:3)';

        case {'antenne_C4', 'antenna_C4_ITM'}
            a_ampl = sqrt(2.6725e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:7)';
        
        case {'antenne_C4_dphi_D'}
            a_ampl = sqrt(2.6725e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(-140*pi/180)*(0:7)';
            
        case {'antenne_C4_mod'}
            a_ampl = sqrt(2.6725e6/16)*ones(8,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:7)';    

        case {'antenne_C4_1module'}
            a_ampl = sqrt(2.6725e6/16)*ones(1,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:0)';      
            
        case {'antenne_C4_perfect_1module'}
            a_ampl = sqrt(2.6725e6/16)*ones(1,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:0)';      
        
        case {'antenne_C4_perfect2_1module'}
            a_ampl = sqrt(2.6725e6/16)*ones(1,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:0)';      
         
        case {'antenne_C4_3guides'}
            a_ampl = sqrt(2.6725e6/16)*ones(1,1);   % En champ elect pour 25 MW/m^2 !!!
            a_phase = 1*(180*pi/180)*(0:0)';    

        case 'antenne_ITER' 
            a_ampl = sqrt(5e6/12)*ones(3,1);   % Antenne ITER
            a_phase = 1*(0*pi/180)*(0:2)';   % Antenne ITER 

        case 'antenne_reduite' % test 
            a_ampl = [];
            a_phase = [];

        case 'antenne_orso1' % test
            a_ampl = sqrt(1)*ones(8,1); % 1 W by waveguide, 8 wgs
            a_phase = (90*pi/180)*(0:7)';      

        case 'antenne_orso2' % test
            a_ampl = ones(2,1);
            a_phase = (pi/180)*(0:1)';    

        case {'antenne_PAM_TOPLHA', 'antenne_C4_1module'} % test
            a_ampl = sqrt(1); % 1 W
            a_phase = 0;
            
        case {'antenne_PAM_TOPLHA_4modules'} % test
            a_ampl = sqrt(0.113e6)*ones(4,1);
            a_phase = zeros(4,1);

        case {'antenne_JET', 'antenne_JET_1row'} % 8 modules of 4(t)x2(p) waveguides each
           a_ampl = sqrt(1e3*130)*ones(8,1);  % 25 MW/m^2 
           a_phase = 0*(pi/180)*(0:7)';
           
  
            
        otherwise
          error(aloha_message(['Antenna excitation unknow for architecture : ', architecture, '. Use an home made excitation']));
    end

% 
% 4 arguments : il s'agit de la recuperation des donnees isssues 
% d'un choc
% 
elseif (nargin == 4)
    choc = varargin{1};
    tps_1 = varargin{2};
    tps_2 = varargin{3}; 
    TSport = varargin{4};
    
    switch (TSport)
        case 'Q6A'
            [a_ampl, a_phase] = retrieve_tsbase_Q6A(choc, tps_1, tps_2);

        case 'Q6B'
            [a_ampl, a_phase] = retrieve_tsbase_Q6B(choc, tps_1, tps_2);

        otherwise
            error('Bad argument : TSport')
    end
else
    error('Bad number of arguments!');
end



function [a_ampl, a_phase] = retrieve_tsbase_Q6A(choc, tps_1, tps_2)
    % Lecture des puissances incidentes :
    disp(aloha_message('Reading power measurements on Q6A [gpinjc1]'));
    [p_inc_mes,tps_pow]=tsbase(choc,'gpinjc1');
    % MODIF 12/11/2008 JH
    % Le dephasage du aux rallonges (entre fenetres et jonction hybride) a ete
    % prise en compte dans l'architecture de C2. Par consequent, la lecture de la phase
    % incident doit se faire au niveau des fenetres. La phase en bout etant calculee grace
    % aux dephasages des rallonges et a la modelisation HFSS de C2
    % 
    % Lecture des phases incidentes (entree du module):
    disp(aloha_message('Reading phase measurements on Q6A (before RF windows) [gphic1]'));
    [phase_inc_mes,tps_pha]=tsbase(choc,'gphic1');
 
    % keep only one and positive time vector
    tps_pow=tps_pow(:,1);
    tps_pow_pos=find(tps_pow>0);
    tps_pha=tps_pha(:,1);
    tps_pha_pos=find(tps_pha>0);

    % keep measured data for positive time
    p_inc_mes=p_inc_mes(tps_pow_pos,:);
    tps_pow = tps_pow(tps_pow_pos);
    phase_inc_mes=phase_inc_mes(tps_pha_pos,:);
    tps_pha = tps_pha(tps_pha_pos);

    % Lecture des puissances reflechies :
    [p_refl_mes,tps_pow2]=tsbase(choc,'gcrefc1');
    % Lecture des phases reflechies (entr??? du module):
    [phase_refl_mes,tps_pha2]=tsbase(choc,'gphrc1');

    % keep only one and positive time vector
    tps_pow2=tps_pow2(:,1);  
    tps_pow2_pos=find(tps_pow2>0);
    tps_pha2=tps_pha2(:,1);  
    tps_pha2_pos=find(tps_pha2>0);

    % keep measured data for positive time
    p_refl_mes=p_refl_mes(tps_pow2_pos,:);
    tps_pow2 = tps_pow2(tps_pow2_pos);
    phase_refl_mes=phase_refl_mes(tps_pha2_pos,:);
    tps_pha2 = tps_pha2(tps_pha2_pos);

    % calculate the average of the data between t1 and t2
    p_inc_mes = mean(p_inc_mes(tps_pow>tps_1 & tps_pow<tps_2,:));
    phase_inc_mes = mean(phase_inc_mes(tps_pha>tps_1 & tps_pha<tps_2,:));
    p_refl_mes= mean(p_refl_mes(tps_pow2>tps_1 & tps_pow2<tps_2,:));
    phase_refl_mes= mean(phase_refl_mes(tps_pha2>tps_1 & tps_pha2<tps_2,:));

%  % OLD CODE
%      tps=tps(:,1);  
%      tps_pos=find(tps>0);
%      p_refl_mes=p_refl_mes(tps_pos,:);
%      phase_refl_mes=phase_refl_mes(tps_pos,:);
%      
%      % calcul des moyennes sur l'intervalle [t1 t2]:
%      
%      p_inc_complex=p_inc_mes.*exp(i*phase_inc_mes*pi./180);
%      p_inc_real=mean(real(p_inc_complex(tps>tps_1&tps<tps_2,:)));
%      p_inc_imag=mean(imag(p_inc_complex(tps>tps_1&tps<tps_2,:)));
%      p_inc_mes=abs(p_inc_real+i*p_inc_imag);
%      phase_inc_mes=angle(p_inc_real+i*p_inc_imag);
%      
%      p_refl_complex=p_refl_mes.*exp(i*phase_refl_mes*pi./180);
%      p_refl_real=mean(real(p_refl_complex(tps>tps_1&tps<tps_2,:)));
%      p_refl_imag=mean(imag(p_refl_complex(tps>tps_1&tps<tps_2,:)));
%      p_refl_mes=abs(p_refl_real+i*p_refl_imag);
%      phase_refl_mes=angle(p_refl_real+i*p_refl_imag);

    % rangement des klystrons ds le code ALOHA, 
    % de Gauche a Droite face a l'antenne vu du plasma :
    % indices 1 -> 8  coupleur haut 
    % indices 9 -> 16 coupleur bas 
    p_inc_mes=p_inc_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    p_refl_mes=p_refl_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    phase_inc_mes=phase_inc_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    phase_refl_mes=phase_refl_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);

%      % coupleur haut 
%      a_ampl = sqrt(p_inc_mes(1:8)*1e3)';
%      a_phase = pi/180*phase_inc_mes(1:8)';
    % coupleur bas 
    a_ampl = sqrt(p_inc_mes(9:16)*1e3)';
    a_phase = pi/180*phase_inc_mes(9:16)';


function [a_ampl, a_phase] = retrieve_tsbase_Q6B(choc, tps_1, tps_2)
    % Lecture des puissances incidentes :
    disp(aloha_message('Reading power measurements on Q6B [gpinjc2]'));
    [p_inc_mes,tps_pow]=tsbase(choc,'gpinjc2');
    
    % Lecture des phases incidentes (entree du module):
    disp(aloha_message('Reading phase measurements on Q6B (before RF windows) [gphic2]'));
    % 'gphic2' : phase before the antenna
    % 'gphbc2' : phase calulated at the mouth of the antenna
    [phase_inc_mes,tps_pha]=tsbase(choc,'gphic2');  

    % keep only one and positive time vector
    tps_pow=tps_pow(:,1);
    tps_pow_pos=find(tps_pow>0);
    tps_pha=tps_pha(:,1);
    tps_pha_pos=find(tps_pha>0);

    % keep measured data for positive time
    p_inc_mes=p_inc_mes(tps_pow_pos,:);
    tps_pow = tps_pow(tps_pow_pos);
    phase_inc_mes=phase_inc_mes(tps_pha_pos,:);
    tps_pha = tps_pha(tps_pha_pos);

    % read measured reflected power
    [p_refl_mes,tps_pow2]=tsbase(choc,'gcrefc2');
    % read measured reflected phase
    [phase_refl_mes,tps_pha2]=tsbase(choc,'gphrc2');

    % keep only one and positive time vector
    tps_pow2=tps_pow2(:,1);  
    tps_pow2_pos=find(tps_pow2>0);
    tps_pha2=tps_pha2(:,1);  
    tps_pha2_pos=find(tps_pha2>0);

    % keep measured data for positive time
    p_refl_mes=p_refl_mes(tps_pow2_pos,:);
    tps_pow2 = tps_pow2(tps_pow2_pos);
    phase_refl_mes=phase_refl_mes(tps_pha2_pos,:);
    tps_pha2 = tps_pha2(tps_pha2_pos);

    % calculate the average of the data between t1 and t2
    p_inc_mes = mean(p_inc_mes(tps_pow>tps_1 & tps_pow<tps_2,:));
    phase_inc_mes = mean(phase_inc_mes(tps_pha>tps_1 & tps_pha<tps_2,:));
    p_refl_mes= mean(p_refl_mes(tps_pow2>tps_1 & tps_pow2<tps_2,:));
    phase_refl_mes= mean(phase_refl_mes(tps_pha2>tps_1 & tps_pha2<tps_2,:));

%  % OLD CODE
%      p_inc_complex=p_inc_mes.*exp(i*phase_inc_mes*pi./180);
%      p_inc_real=mean(real(p_inc_complex(tps>tps_1&tps<tps_2,:)));
%      p_inc_imag=mean(imag(p_inc_complex(tps>tps_1&tps<tps_2,:)));
%      p_inc_mes=abs(p_inc_real+i*p_inc_imag);
%      phase_inc_mes=angle(p_inc_real+i*p_inc_imag);
%      
%      p_refl_complex=p_refl_mes.*exp(i*phase_refl_mes*pi./180);
%      p_refl_real=mean(real(p_refl_complex(tps>tps_1&tps<tps_2,:)));
%      p_refl_imag=mean(imag(p_refl_complex(tps>tps_1&tps<tps_2,:)));
%      p_refl_mes=abs(p_refl_real+i*p_refl_imag);
%      phase_refl_mes=angle(p_refl_real+i*p_refl_imag);
    
    
    % rangement des klystrons ds le code ALOHA, 
    % de Gauche a Droite face a l'antenne vu du plasma :
    % indices 1 -> 8  coupleur haut 
    % indices 9 -> 16 coupleur bas 
    p_inc_mes=p_inc_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    p_refl_mes=p_refl_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    phase_inc_mes=phase_inc_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);
    phase_refl_mes=phase_refl_mes([15 13 11 9 7 5 3 1 16 14 12 10 8 6 4 2]);

    % coupleur haut
    a_ampl = sqrt(p_inc_mes(1:8)*1e3)';
    a_phase = pi/180*phase_inc_mes(1:8)';
%      % coupleur bas
%      a_ampl = sqrt(p_inc_mes(9:16)*1e3)';
%      a_phase = pi/180*phase_inc_mes(9:16)';