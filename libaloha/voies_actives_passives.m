

pass_ligne = 1:nb_g_passifs_bord;
act_ligne = [];

for ind = 1:nb_modules_tor-1

    pass_ligne = [pass_ligne,pass_module_tor + nb_g_passifs_bord + (ind-1)*(nb_g_module_tor+nb_g_passifs_inter_modules)];   
  
    pass_ligne = [pass_ligne,(1:nb_g_passifs_inter_modules) + nb_g_module_tor + ...
                                    		nb_g_passifs_bord + (ind-1)*(nb_g_module_tor+nb_g_passifs_inter_modules)];

			    
    act_ligne = [act_ligne,act_module_tor + nb_g_passifs_bord + (ind-1)*(nb_g_module_tor+nb_g_passifs_inter_modules)];
    
end

pass_ligne = [pass_ligne,pass_module_tor + nb_g_passifs_bord + (nb_modules_tor-1)*(nb_g_module_tor+nb_g_passifs_inter_modules)];   
act_ligne = [act_ligne,act_module_tor + nb_g_passifs_bord + (nb_modules_tor-1)*(nb_g_module_tor+nb_g_passifs_inter_modules)];

pass_ligne = [pass_ligne,(1:nb_g_passifs_bord) + (nb_g_total_ligne - nb_g_passifs_bord)];


pass_tot = [];
act_tot = [];

for ind = 1:nb_g_pol

    pass_tot = [pass_tot,pass_ligne + (ind-1)*nb_g_total_ligne];
    act_tot = [act_tot,act_ligne + (ind-1)*nb_g_total_ligne];
    
end


act_module_tor = (nb_g_module_tor-length(pass_module_tor));

act_tot = reshape(act_tot,nb_modules_tor*act_module_tor,nb_g_pol)';

pass_tot_tamp = pass_tot;
pass_tot = pass_tot*(Nme+Nmh) - (Nme+Nmh-1);

incr = 1;
modules_act = [];

for ind = 1: (nb_g_pol/nb_g_module_pol)

    for k = 1:nb_modules_tor

    tampon = act_tot(1+(ind-1)*nb_g_module_pol:ind*nb_g_module_pol,1+(k-1)*act_module_tor:k*act_module_tor);
	tampon = reshape(tampon',1,nb_g_module_pol*act_module_tor);
	modules_act(incr,:) = tampon;
	incr = incr + 1;
	
    end

end

act_tot_tamp = act_tot;
act_tot = act_tot*(Nme+Nmh) - (Nme+Nmh-1);

modules_act_tamp = modules_act;
modules_act = modules_act*(Nme+Nmh) - (Nme+Nmh-1);


%disp(act_tot);
%disp(modules_act)






