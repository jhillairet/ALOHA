% matrice S theorique du module "02B"
% provenant du document TS9789 du 11/4/89 de Ph.Bibet et J.Achard

% module
S_mod11 = [0.03 0.35 0.35 0.35 0.35; ...
           0.35 0.50 0.62 0.25 0.25; ...
					 0.35 0.62 0.50 0.25 0.25; ...
					 0.35 0.25 0.25 0.50 0.62; ...
					 0.35 0.25 0.25 0.62 0.50];
					 
S_mod12 = [0.35*ones(1,4);zeros(4,4)];

S_mod = [S_mod11,  S_mod12; ...
         S_mod12', S_mod11(2:5,2:5)];
				 					
% phase en degres						
S_phase11 = [+034 -098 -008 +088 +178; ...
             -098 -013 -160 +066 +160;  ...
						 -008 -160 +170 +160 +110; ...
						 +088 +066 +160 -014 -160; ...
						 +178 +160 +110 -160 +170];
						
S_phase22 = [-100 +110 -025 +065; ...
             +110 +076 +065 +160; ...
						 -025 +065 -100 +110; ...
						 +065 +160 +110 +076];
						
S_phase12 = [+082 +172 -088 +002; zeros(4,4)];
						
S_phase = [S_phase11, S_phase12; ...
					 S_phase12', S_phase22];

% matrice S complexe
S = S_mod .* exp(i*pi./180.*S_phase);					 

S = reshape(S,1,length(S)*length(S)); 
