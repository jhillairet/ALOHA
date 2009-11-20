function aloha_plot_spectrum(scenarios)  
% Affiche le spectre en 1D et 2D
%  
% EXAMPLE
%  aloha_plot_spectrum(scenario)
%  
% INPUT
%  - scenario [struct(1)] : structure scenario
%  
% AUTHOR: JH
% LAST CHANGES
%  - 23/10/2009 works with ALOHA-2D
%  - 10/09/2008 (démarrage du LHC!) 
%    

for id_scen=1:length(scenarios) % for many scenarios
    scenario=scenarios(id_scen);

    % load scenario var into workspace 
    % (for compatibility reasons, it's simpler than to change evvery variables in the code...)
    aloha_scenario_loadIntoWorkspace

    if ~exist('dP')
        error(aloha_message('Le spectre n''a pas ete calcule !'))
    end

    letitre1 = ['Spectre ', aloha_utils_str4fig(architecture)];


    if (bool_mesure)
        letitre1 = strcat(letitre1,' du choc ', num2str(choc),' de ', num2str(tps_1), ' a ', num2str(tps_2), ' s');
    else
        letitre1 = strcat(letitre1,' sans mesure (cas ideal)');
    end
    if aloha_scenario_get(scenario, 'bool_lignes_identiques')
        letitre2 = ['ne0 = ',num2str(ne0(1),'%1.2e'),', d\_couche = ',num2str(d_couche(1),'%1.2e'),','];
        switch aloha_scenario_get(scenario, 'version')
            case 1 || 2 || 3
                letitre3 = ['dne0 = ', num2str(dne0(1),'%1.2e')];
            case 6
                letitre3 = ['dne0 = ', num2str(dne0(1),'%1.2e'),' et dne1 = ',num2str(dne1(1),'%1.2e') ];
        end
    else
        letitre2 = strcat('ne0 = [',num2str(ne0(1),'%1.2e'),';',num2str(ne0(2),'%1.2e'),'], d-couche = [',num2str(d_couche(1),'%1.2e'),';',num2str(d_couche(2),'%1.2e'),']');
        if (version == 3)
            letitre3 = strcat('dne0 = [',num2str(dne0(1),'%1.2e'),';',num2str(dne0(2),'%1.2e'),']');
        elseif (version == 6)
            letitre3 = strcat('dne0 = [',num2str(dne0(1),'%1.2e'),';',num2str(dne0(2),'%1.2e'),'] et dne1 = [',num2str(dne1(1),'%1.2e'),';',num2str(dne1(2),'%1.2e'),']' );
        end
    end

        h=aloha_plot_figure(figure, 'ALOHA : nz spectrum');
        plot(nz, real(dP_nz));
        xlabel('n_{//}');
        ylabel('Power density');
        title( {letitre1;letitre2} );  
%          text( nz(1) , max(real(dP(find(ny==0),:)))*.95 , strcat('Nmh(TE) = ',num2str(Nmh)) );
%          text( nz(1) , max(real(dP(find(ny==0),:)))*.88 , strcat('Nme(TM) = ',num2str(Nme)) );


    % figure spectre 2D  
    % Rajout le 11/06/2007 par Izacard Olivier :
        h=aloha_plot_figure(figure, 'ALOHA : ny,nz spectrum');
         
        mesh(nz , ny , real(dP));
        xlabel('n_z');
        ylabel('n_y');
        zlabel('dP');


%      % calcul et affichage du spectre de l'antenne entiere
%      % (spectre des deux demi-antennes combine)
%      % 
%      % TODO
%      dP_all = aloha_spectrumAllAntenna(freq, ny, dP, 3*a+2*espacement_g_pol+80e-3);
%      [NYY,NZZ]=ndgrid(ny,nz);
%  
%      figure
%      subplot(211)
%          pcolor(NZZ,NYY, real(dP)); 
%          shading interp;
%          xlabel('n_z'); ylabel('n_y');
%          axis equal; colorbar; cax=caxis;
%      subplot(212)
%          pcolor(NZZ,NYY, real(dP_all)); 
%          shading interp;
%          xlabel('n_z'); ylabel('n_y');
%          axis equal; colorbar; caxis(cax);
end % id_scen
