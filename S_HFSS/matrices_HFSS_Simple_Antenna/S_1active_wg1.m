% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,2,2);

f = [4600000000 ];
%  Version HFSS
S(1,:,:) = [   -2.782861E-004 -    8.239460E-004i,   -2.342642E-001 -    9.721726E-001i;  -2.342642E-001 -    9.721726E-001i,    1.275582E-004 +    8.602668E-004i];
%  % Version theorique
%  S(1,:,:) = [   -0,   1;  1,    0];