function aloha_ondemand_plotPhase(TSpulse, TSport, TStime1, TStime2)
% ALOHA
% 
% aloha_ondemand_plotPhase(TSpulse, TSport)
%  or
% aloha_ondemand_plotPhase(TSpulse, TSport, tstart, tstop)
% 
% Plot the phase difference between modules for both upper and lower modules for a TS pulse
% 
% INPUT
%  - TSpulse : TS pulse number
%  - TSport : 'Q6A' or 'Q6B'
%  - [optionnal] tstart : time to start the display
%  - [optionnal] tstop : time to stop the display
% 
% OUTPUT : none
% 
% AUTHOR: JH
%  
%  TODO : get the requested phase shift automatically


%% Select the appropriate signal to plot depending on the TS port
switch upper(TSport)
    case 'Q6A' % C2 or C3
        sig_phase = 'GPHIC1';
        sig_ahyb = 'AHYB;Ph_Cp1;VALEURS';

    case 'Q6B' % C3 or C4
        sig_phase = 'GPHIC2';
        sig_ahyb = 'AHYB;Ph_Cp2;VALEURS';

    otherwise
        error('bad port definition. See help');
end

%% retrieve the phase signal from the TS database
[phase, t_phase] = tsbase(TSpulse, sig_phase);
t_phase = t_phase(:,1);
phase = phase*pi/180; % -> put in radians
% retrieve the requested values from the TS database
phase_request=tsmat(TSpulse,sig_ahyb)
% TODO : find a better way to get this value ! in correspondance with the TStime1 & 2 !
defPhaseShift = phase_request(1,2);%interp1(phase_request(:,1), phase_request(:,2), TStime1)


% if tstart and tstop provided -> crop the signal
cropIdx = find(t_phase >= TStime1 & t_phase <= TStime2);
t_phase = t_phase(cropIdx);
phase = phase(cropIdx, :);

%% average the phase
phase_avg = mean(phase);

upperModules = [15 13 11  9 7 5 3 1];
lowerModules = [16 14 12 10 8 6 4 2];

%% generate data for upper and lower modules
for idx=1:2
        if idx==1 % upper modules
            PP{idx} = phase_avg(:,upperModules);
        elseif idx == 2 % lower modules
            PP{idx} = phase_avg(:,lowerModules);
        end
        
        % wrap the phase
        PP2{idx} = (PP{idx}<pi).*PP{idx} + (PP{idx}>=pi).*(PP{idx}-2*pi);
        % calculates the ideal phase progression between modules
        idealPhase{idx} = angle(exp(i*(PP2{idx}(1)+[0:7]*defPhaseShift*pi/180)));

        % calculates the difference to this ideal phase progression
        phaseDiff{idx} = idealPhase{idx}-PP2{idx};
        % correct bot the +/-360 ambiguity
        PP2{idx}(phaseDiff{idx}<-pi)=PP2{idx}(phaseDiff{idx}<-pi)-2*pi;
        PP2{idx}(phaseDiff{idx}>+pi)=PP2{idx}(phaseDiff{idx}>+pi)+2*pi;

        phaseDiff{idx} = idealPhase{idx}-PP2{idx};
        % unwrapped phase
        PP2_unwrap{idx}=unwrap(PP2{idx});       
        idealPhase_unwrap{idx} = unwrap(idealPhase{idx});
        phaseDiff_unwrap{idx}=unwrap(phaseDiff{idx});
end

%% plot phase at the output of each module
aloha_plot_figure(figure)
for idx=1:2
    subplot(2,1,idx);
        plot(1:8, PP2{idx}*180/pi, '.', ...
            1:8,  idealPhase{idx}*180/pi, 'k.--', 'MarkerSize',30);   
        ylabel('Phase [deg]');
        grid on;
        if idx == 1 
            title(['TS pulse #',num2str(TSpulse), ' - TSport ', TSport, ' - \delta\phi=',num2str(defPhaseShift),' deg - Module output phase [deg]']);
            set(gca, 'XTickLabel', {upperModules});
            xlabel('Upper Modules #');
        elseif idx == 2
            set(gca, 'XTickLabel', {lowerModules});
            xlabel('Lower Modules #');
        end
        % write the phase difference on the plot
        for idx2=1:8
            text(idx2, idealPhase{idx}(idx2)*180/pi+10, num2str(phaseDiff{idx}(idx2)*180/pi))
        end
end

%% Same thing unwrapped phase
aloha_plot_figure(figure)
for idx=1:2
    subplot(2,1,idx);
        plot(1:8, PP2_unwrap{idx}*180/pi, '.', ...
            1:8,  idealPhase_unwrap{idx}*180/pi, 'k.--', 'MarkerSize',30);   
        ylabel('Phase [deg]');
        grid on;
        if idx == 1 
            title(['TS pulse #',num2str(TSpulse), ' - TSport ', TSport, ' - \delta\phi=',num2str(defPhaseShift),' deg - Module output phase (unwrapped) [deg]']);
            set(gca, 'XTickLabel', {upperModules});
            xlabel('Upper Modules #');
        elseif idx == 2
            set(gca, 'XTickLabel', {lowerModules});
            xlabel('Lower Modules #');
        end
        % write the phase difference on the plot
        for idx2=1:8
            text(idx2, idealPhase_unwrap{idx}(idx2)*180/pi+10, num2str(phaseDiff{idx}(idx2)*180/pi))
        end
end


%% Display ALOHA input phase values
disp('Top modules  --- Bot modules');
[PP2{1}'*180/pi, PP2{2}'*180/pi]
disp('Phase Difference between modules - Top modules  --- Bot modules')
difference_H = ( PP2{1}(2:end)- PP2{1}(1:end-1) ).' *180/pi;
difference_B = ( PP2{2}(2:end)- PP2{2}(1:end-1) ).' *180/pi;
[difference_H, difference_B]
unwrap([difference_H, difference_B]*pi/180)*180/pi
disp('average phase difference between modules')
mean([unwrap(difference_H*pi/180)*180/pi, unwrap(difference_B*pi/180)*180/pi])


