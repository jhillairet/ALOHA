function aloha_plot_champEmbouchure1D(scenario)
  % load scenario var into workspace 
  % (for compatibility reasons, it's simpler than to change evvery variables in the code...)
  aloha_scenario_loadIntoWorkspace
  % Load into workspace the geometrical parameter of the scenario architecture 
  % and get the main geometrical parameters.
  eval(architecture);

  %  Trace le champ dans l'embouchure (si il existe)
  %  
  if not(exist('abs_z')) || not(exist('Efield'))
    error('Fields haven''t been calculated !')
  end

  h=aloha_plot_figure(figure, 'ALOHA : Parallel Electric field amplitude (Ez) in the antenna mouth');

% for all poloidal lines
  for idx_pol = 1:nb_g_pol 
    subplot(nb_g_pol,1,idx_pol)
    plot(abs_z(idx_pol,:), abs(Efield(idx_pol,:)));
    grid on;
    xlabel('z [m]');
    ylabel('|E_z| (V/m)');
  end % idx_pol


  letitre1 = ['Electric field amplitude in ', aloha_utils_str4fig(architecture), ' mouth.'];

  if bool_mesure
    letitre1 = strcat(letitre1,' du choc ',num2str(choc),' de ',num2str(tps_1),' a ',num2str(tps_2),'s');
  else
    letitre1 = strcat(letitre1,' sans mesure (cas ideal)');;
  end

  if bool_lignes_identiques
    letitre2 = strcat('ne0 = ',num2str(ne0(1),'%1.2e'),', d-couche = ',num2str(d_couche(1),'%1.2e'),',');
    letitre3 = strcat('dne0 = ',num2str(dne0(1),'%1.2e'),' et dne1 = ',num2str(dne1(1),'%1.2e') );
  else
    letitre2 = strcat('ne0 = [',num2str(ne0(1),'%1.2e'),';',num2str(ne0(2),'%1.2e'),'], d-couche = [',num2str(d_couche(1),'%1.2e'),';',num2str(d_couche(2),'%1.2e'),']');
    letitre3 = strcat('dne0 = [',num2str(dne0(1),'%1.2e'),';',num2str(dne0(2),'%1.2e'),'] et dne1 = [',num2str(dne1(1),'%1.2e'),';',num2str(dne1(2),'%1.2e'),']' );
  end


  title({letitre1;letitre2;letitre3}); 


% phase
  h=aloha_plot_figure(figure, 'ALOHA parallel electric field phase in the antenna mouth');
  % for all poloidal lines
  nb_g_pol=size(Efield,1);
  for idx_pol = 1:nb_g_pol 
    subplot(nb_g_pol,1,idx_pol)
    plot(abs_z(idx_pol,:), 180/pi*angle(Efield(idx_pol,:)), ...
         abs_z(idx_pol,:), 180/pi*Efield_average_phase(idx_pol,:));
    grid on;
    xlabel('z [m]')
    ylabel('\angle E_z (degree)')
  end % idx_pol