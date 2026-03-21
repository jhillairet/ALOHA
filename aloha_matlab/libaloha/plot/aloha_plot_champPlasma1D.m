function aloha_plot_champPlasma1D(scenario)
%  Plot the spectrum-propagated Ez field into the plasma
%  
%  INPUT
%   - scenario : ALOHA scenario
%  OUPUT
%   none
%   
% AUTHORS: D.Voyer / O.Izacard / J.Hillairet
% LAST UPDATE:
%  - 04/2009: use scenario

h=aloha_plot_figure(figure, '-= ALOHA : Electric field in plasma =-');

% retrieve computed fields parameters in the input scenario
[x_coord, z_coord, Ez_x_z, Hy_x_z] = aloha_scenario_get(scenario, 'x_coord', 'z_coord', 'Ez_x_z', 'Hy_x_z');

z_coord=squeeze(z_coord);
Ez_x_z=squeeze(Ez_x_z);
Hy_x_z=squeeze(Hy_x_z);

if (aloha_scenario_get(scenario, 'fig_Ez_ou_EzHy') == 1)
       % Rajout le 16/05/2007 par Izacard Olivier pour les titres :
       letitre0 = '';
       %pcolor(z_coord, x_coord, log10(1+abs(Ez_x_z)));
       pcolor(z_coord, x_coord, abs(Ez_x_z));
else
       % Rajout le 16/05/2007 par Izacard Olivier pour les titres :
       letitre0 = 'Vecteur de poynthing dans le plasma (V/m) a l'' embouchure';
       pcolor(z_coord, x_coord, 20*log10(abs(Ez_x_z.*(Hy_x_z.')')));
end
     
   shading flat
   shading interp
   colorbar
   ylabel('x [m]')
   xlabel('z [m]')
   axis equal;
    
    letitre1 = ['Electric field in plasma for ', aloha_utils_str4fig(aloha_scenario_get(scenario, 'architecture'))];

    if(aloha_scenario_get(scenario, 'bool_mesure'))
       letitre1 = [letitre1, ...
        ' from pulse ',num2str(aloha_scenario_get(scenario, 'choc')), ...
        ' from ',num2str(aloha_scenario_get(scenario, 'tps_1')), ...
        ' to ',num2str(aloha_scenario_get(scenario, 'tps_2')),'s'];
    else
       letitre1 = strcat(letitre1,' sans mesure (cas ideal)');;
    end

%  [ne0, dne, d_couche] = aloha_scenario_get(scenario, 'ne0', 'dne0', 'dne1', 'd_couche');
%      
%      if(aloha_scenario_get(scenario, 'bool_lignes_identiques') )
%         letitre2 = strcat('ne0 = ',num2str(ne0(1),'%1.2e'),', d-couche = ',num2str(d_couche(1),'%1.2e'),',');
%         letitre3 = strcat('dne0 = ',num2str(dne0(1),'%1.2e'),' et dne1 = ',num2str(dne1(1),'%1.2e') );
%      else
%         letitre2 = strcat('ne0 = [',num2str(ne0(1),'%1.2e'),';',num2str(ne0(2),'%1.2e'),'], d-couche = [',num2str(d_couche(1),'%1.2e'),';',num2str(d_couche(2),'%1.2e'),']');
%         letitre3 = strcat('dne0 = [',num2str(dne0(1),'%1.2e'),';',num2str(dne0(2),'%1.2e'),'] et dne1 = [',num2str(dne1(1),'%1.2e'),';',num2str(dne1(2),'%1.2e'),']' );
%      end
    title( {strcat(letitre0,letitre1)});%;letitre2;letitre3} );

