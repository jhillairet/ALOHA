%  Calcule les vecteurs qui definissent les dimensions de l'antenne


% nombre total de guides par ligne
nb_g_total_ligne = nb_g_module_tor*nb_modules_tor + 2*nb_g_passifs_bord + nb_g_passifs_inter_modules*(nb_modules_tor - 1);
%  nombre total de guides 
nbre_guides=nb_g_total_ligne*nb_g_pol;

b_module = zeros(1,nb_g_module_tor);	% largeur des guides dans le sens toroidal

act_module_tor = 1:nb_g_module_tor;
act_module_tor(:,pass_module_tor) = [];

if (antenne_standard == 1) 

   b_module(act_module_tor) = b_g_actif;
   b_module(pass_module_tor) = b_g_pass;

   b_bord = b_g_pass*ones(1,nb_g_passifs_bord);			% largeur des guides passifs du bord dans le sens toroidal
   b_inter_modules = b_g_pass*ones(1,nb_g_passifs_inter_modules);	% largeur des guides passifs du bord dans le sens toroidal

   b = [b_bord,kron(ones(1,nb_modules_tor-1),[b_module,b_inter_modules]),b_module,b_bord];

   z = zeros(1,nb_g_total_ligne);

   for ind = 2:nb_g_total_ligne

       z(1,ind) = z(1,ind-1) + b(ind-1) + e;
	
   end
   
else
   
   z = zeros(1,nb_g_total_ligne);
   for ind = 2:nb_g_total_ligne

       z(1,ind) = z(1,ind-1) + b(ind-1) + e(ind-1);
	
   end
    
end

h = espacement_g_pol*ones(1,nb_g_pol-1);

y = zeros(1,nb_g_pol);

for ind = 2:nb_g_pol

    y(ind) = y(ind-1) +  h(ind-1) + a;
    
end

