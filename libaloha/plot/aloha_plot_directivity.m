
function aloha_plot_directivity(scenario)
% Plot the antennea direcivity
% 
% EXAMPLE
%  aloha_plot_directivity(scenario)
%  
% INPUT
%  - scenario [struct(1)] : structure scenario
%   
% NB : in order to plot the antenna directivity, 
% the antenna spectrum must have been computed.
%   
% AUTHOR: JH
% LAST CHANGES
%  - 10/09/2008 (démarrage du LHC!) 
%    


    % load scenario var into workspace 
    % (for compatibility reasons, it's simpler than to change evvery variables in the code...)
    aloha_scenario_loadIntoWorkspace


    if ~exist('directivite_cumulee')
        error(aloha_message('La directivite n''a pas ete calculee !'));
    end

    % Rajout le 16/05/2007 par Izacard Olivier pour les titres :
        letitre1 = ['Directivity of ', aloha_utils_str4fig(architecture)];

        if(bool_mesure)
           letitre1 = [letitre1,' du choc ',num2str(choc),' de ',num2str(tps_1),' a ',num2str(tps_2),'s'];
        else
           letitre1 = [letitre1,' sans mesure (cas ideal)'];
        end

        if(bool_lignes_identiques)
           letitre2 = strcat('n_{e0} = ',num2str(ne0(1),'%1.2e'),', d-couche = ',num2str(d_couche(1),'%1.2e'),',');
           letitre3 = strcat('dne0 = ',num2str(dne0(1),'%1.2e'),' et dne1 = ',num2str(dne1(1),'%1.2e') );
        else
           letitre2 = strcat('ne0 = [',num2str(ne0(1),'%1.2e'),';',num2str(ne0(2),'%1.2e'),'], d-couche = [',num2str(d_couche(1),'%1.2e'),';',num2str(d_couche(2),'%1.2e'),']');
           letitre3 = strcat('dne0 = [',num2str(dne0(1),'%1.2e'),';',num2str(dne0(2),'%1.2e'),'] et dne1 = [',num2str(dne1(1),'%1.2e'),';',num2str(dne1(2),'%1.2e'),']' );
        end


        h=aloha_plot_figure(figure, '-= ALOHA : directivity =-');  
  
        plot(nz, directivite_cumulee)
        xlabel('indice n_z');
        ylabel('directivite normalisee');

        title( {letitre1;letitre2;letitre3} );  %title('directivite de l''antenne')
        text( dnz , max(directivite_cumulee)*.95 , strcat('Nmh(TE) = ',num2str(Nmh)) );
        text( dnz , max(directivite_cumulee)*.78 , strcat('Nme(TM) = ',num2str(Nme)) );

