% matrice S de la multijonction
S_mod_multi_jonction = [0.028  0.5    0.5    0.5    0.5
                0.5    0.448  0.652  0.249  0.249
                0.5    0.652  0.448  0.249  0.249
                0.5    0.249  0.249  0.448  0.652
                0.5    0.249  0.249  0.652  0.448 ];



S_phase_multi_jonction = [ 2.61    -112.84  -22.84    67.16    157.16
                -112.84  -54.80    164.94   34.79    124.79
                -22.84    164.94   125.20   124.79  -145.21
                 67.16    34.79    124.79  -54.80    164.94
                 157.16   124.79  -145.21   164.94   125.20 ];


S_M = S_mod_multi_jonction.*exp(i*S_phase_multi_jonction*pi./180);

% matrice S de la jonction hybride avec une charge adaptée
S_H = (1/sqrt(2))*[0  i  1
                   i  0  0
                   1  0  0];


S_tot = [[S_H(1,1),   zeros(1,4),   zeros(1,4),   S_H(1,2),   S_H(1,3),   0,          0         ];
         [zeros(4,1), S_M(2:5,2:5), zeros(4,4),   zeros(4,1), zeros(4,1), S_M(2:5,1), zeros(4,1)];
         [zeros(4,1), zeros(4,4),   S_M(2:5,2:5), zeros(4,1), zeros(4,1), zeros(4,1), S_M(2:5,1)];        
         [S_H(2,1),   zeros(1,4),   zeros(1,4),   S_H(2,2),   S_H(2,3),   0,          0         ];
         [S_H(3,1),   zeros(1,4),   zeros(1,4),   S_H(3,2),   S_H(3,3),   0,          0         ];
         [0,          S_M(1,2:5),   zeros(1,4),   0,          0,          S_M(1,1),   0         ];	 	 
         [0 ,         zeros(1,4),   S_M(1,2:5),   0,          0,          0,          S_M(1,1)  ] ];


	 
S_connex=[0 0 1 0
          0 0 0 1
	  1 0 0 0
	  0 1 0 0 ];




S = S_tot(1:9,1:9) + S_tot(1:9,10:13)*inv(S_connex - S_tot(10:13,10:13))*S_tot(10:13,1:9);

% prise en cpte de la forme incurvée en bout de grill
S(1,6:9) = S(1,6:9)*exp(i*pi/2);
S(2:5,6:9) = S(2:5,6:9)*exp(i*pi/2);
S(6:9,1) = S(6:9,1)*exp(i*pi/2);
S(6:9,2:5) = S(6:9,2:5)*exp(i*pi/2);
S(6:9,6:9) = S(6:9,6:9)*exp(i*pi);


S = reshape(S,1,length(S)*length(S)); 

