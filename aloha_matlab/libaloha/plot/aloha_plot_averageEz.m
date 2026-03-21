function aloha_plot_averageEz(scenarios)
%  
% Plot the average electric field at the mouth

for idx_sc=1:length(scenarios)
    sc=scenarios(idx_sc);
    
    if isfield(sc.results, 'Ez_average')
        aloha_plot_figure(figure())
            plot(1:length(sc.results.Ez_average), sc.results.Ez_average);
            xlabel('z [m]');
            ylabel('Average Parallel Efield [V/m]');
       
    else
        disp('The average electric field has not been calculated !');
        
    end
end

