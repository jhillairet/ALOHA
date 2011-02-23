% Matlab m-File exported from HFSS12.1.1 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,9,9);

f = [3700000000 ];
S(1,:,:) = [   -8.863874E-003 +    2.789834E-002i,    2.019595E-001 +    2.920204E-001i,   -2.943553E-001 +    1.989317E-001i,   -1.967686E-001 -    2.912201E-001i,    2.931165E-001 -    1.943377E-001i,    2.031092E-001 +    2.912629E-001i,   -2.931616E-001 +    2.005585E-001i,   -1.970731E-001 -    2.911801E-001i,    2.930566E-001 -    1.946498E-001i;   2.019595E-001 +    2.920204E-001i,   -2.628794E-001 +    1.359621E-001i,   -4.880372E-001 +    6.711026E-001i,    7.653699E-002 +    1.017475E-001i,   -1.024783E-001 +    7.569233E-002i,   -7.372740E-002 -    1.051437E-001i,    1.058328E-001 -    7.280673E-002i,    7.088665E-002 +    1.055537E-001i,   -1.062291E-001 +    7.000786E-002i;  -2.943553E-001 +    1.989317E-001i,   -4.880372E-001 +    6.711026E-001i,    2.649212E-001 -    1.308900E-001i,   -1.026274E-001 +    7.548547E-002i,   -7.463248E-002 -    1.033495E-001i,    1.059953E-001 -    7.263779E-002i,    7.170919E-002 +    1.066748E-001i,   -1.063750E-001 +    6.979103E-002i,   -6.890450E-002 -    1.070413E-001i;  -1.967686E-001 -    2.912201E-001i,    7.653699E-002 +    1.017475E-001i,   -1.026274E-001 +    7.548547E-002i,   -2.652902E-001 +    1.374544E-001i,   -4.965110E-001 +    6.664822E-001i,    7.243762E-002 +    1.044711E-001i,   -1.051484E-001 +    7.152261E-002i,   -6.968646E-002 -    1.048177E-001i,    1.054823E-001 -    6.881341E-002i;   2.931165E-001 -    1.943377E-001i,   -1.024783E-001 +    7.569233E-002i,   -7.463248E-002 -    1.033495E-001i,   -4.965110E-001 +    6.664822E-001i,    2.668108E-001 -    1.334326E-001i,   -1.051676E-001 +    7.156667E-002i,   -7.064518E-002 -    1.058373E-001i,    1.054903E-001 -    6.881083E-002i,    6.793144E-002 +    1.061476E-001i;   2.031092E-001 +    2.912629E-001i,   -7.372740E-002 -    1.051437E-001i,    1.059953E-001 -    7.263779E-002i,    7.243762E-002 +    1.044711E-001i,   -1.051676E-001 +    7.156667E-002i,   -2.586140E-001 +    1.307117E-001i,   -4.885848E-001 +    6.733781E-001i,    7.614990E-002 +    1.020542E-001i,   -1.027735E-001 +    7.530448E-002i;  -2.931616E-001 +    2.005585E-001i,    1.058328E-001 -    7.280673E-002i,    7.170919E-002 +    1.066748E-001i,   -1.051484E-001 +    7.152261E-002i,   -7.064518E-002 -    1.058373E-001i,   -4.885848E-001 +    6.733781E-001i,    2.604875E-001 -    1.263765E-001i,   -1.027641E-001 +    7.525754E-002i,   -7.440546E-002 -    1.034760E-001i;  -1.970731E-001 -    2.911801E-001i,    7.088665E-002 +    1.055537E-001i,   -1.063750E-001 +    6.979103E-002i,   -6.968646E-002 -    1.048177E-001i,    1.054903E-001 -    6.881083E-002i,    7.614990E-002 +    1.020542E-001i,   -1.027641E-001 +    7.525754E-002i,   -2.665501E-001 +    1.373044E-001i,   -4.987018E-001 +    6.642849E-001i;   2.930566E-001 -    1.946498E-001i,   -1.062291E-001 +    7.000786E-002i,   -6.890450E-002 -    1.070413E-001i,    1.054823E-001 -    6.881341E-002i,    6.793144E-002 +    1.061476E-001i,   -1.027735E-001 +    7.530448E-002i,   -7.440546E-002 -    1.034760E-001i,   -4.987018E-001 +    6.642849E-001i,    2.681128E-001 -    1.332486E-001i];