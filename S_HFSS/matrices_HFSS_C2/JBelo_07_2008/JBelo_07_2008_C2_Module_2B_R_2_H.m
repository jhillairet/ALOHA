% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,9,9);

f = [3700000000 ];
S(1,:,:) = [   -5.773154E-002 +    2.398109E-003i,   -2.725723E-001 +    2.154446E-001i,   -1.912586E-001 -    2.891170E-001i,    2.882071E-001 -    1.952828E-001i,    1.800569E-001 +    2.969969E-001i,   -3.148073E-001 +    1.711869E-001i,   -1.441327E-001 -    3.271630E-001i,    3.273711E-001 -    1.484017E-001i,    1.315266E-001 +    3.335294E-001i;  -2.725723E-001 +    2.154446E-001i,   -2.651626E-001 -    4.004733E-001i,    1.620031E-001 +    7.337636E-001i,    9.762581E-002 -    6.810376E-002i,    6.293911E-002 +    1.007017E-001i,   -7.695856E-002 +    8.936943E-002i,   -8.248121E-002 -    8.391368E-002i,    8.348040E-002 -    8.381953E-002i,    7.934257E-002 +    8.734198E-002i;  -1.912586E-001 -    2.891170E-001i,    1.620031E-001 +    7.337636E-001i,    1.974481E-001 +    4.388916E-001i,    5.949383E-002 +    1.027891E-001i,   -1.054126E-001 +    5.410037E-002i,   -8.237495E-002 -    8.403123E-002i,    9.036637E-002 -    7.494152E-002i,    7.630853E-002 +    9.004821E-002i,   -9.351067E-002 +    7.153266E-002i;   2.882071E-001 -    1.952828E-001i,    9.762581E-002 -    6.810376E-002i,    5.949383E-002 +    1.027891E-001i,   -2.154490E-001 -    3.957706E-001i,    1.165839E-001 +    7.602537E-001i,    8.346863E-002 -    8.365750E-002i,    7.626325E-002 +    8.991321E-002i,   -8.957850E-002 +    7.763248E-002i,   -7.287501E-002 -    9.310871E-002i;   1.800569E-001 +    2.969969E-001i,    6.293911E-002 +    1.007017E-001i,   -1.054126E-001 +    5.410037E-002i,    1.165839E-001 +    7.602537E-001i,    1.780513E-001 +    4.151120E-001i,    7.914947E-002 +    8.737936E-002i,   -9.342810E-002 +    7.145746E-002i,   -7.283858E-002 -    9.316357E-002i,    9.644133E-002 -    6.792056E-002i;  -3.148073E-001 +    1.711869E-001i,   -7.695856E-002 +    8.936943E-002i,   -8.237495E-002 -    8.403123E-002i,    8.346863E-002 -    8.365750E-002i,    7.914947E-002 +    8.737936E-002i,    4.568928E-001 -    3.290695E-001i,   -5.716812E-001 -    2.246554E-002i,    2.668357E-002 +    2.920780E-001i,   -2.923039E-001 +    1.202725E-002i;  -1.441327E-001 -    3.271630E-001i,   -8.248121E-002 -    8.391368E-002i,    9.036637E-002 -    7.494152E-002i,    7.626325E-002 +    8.991321E-002i,   -9.342810E-002 +    7.145746E-002i,   -5.716812E-001 -    2.246554E-002i,   -5.053571E-001 +    2.516153E-001i,   -2.926003E-001 +    2.345922E-003i,    1.224445E-002 -    2.916114E-001i;   3.273711E-001 -    1.484017E-001i,    8.348040E-002 -    8.381953E-002i,    7.630853E-002 +    9.004821E-002i,   -8.957850E-002 +    7.763248E-002i,   -7.283858E-002 -    9.316357E-002i,    2.668357E-002 +    2.920780E-001i,   -2.926003E-001 +    2.345922E-003i,    4.669608E-001 -    2.813941E-001i,   -5.856673E-001 -    5.440046E-002i;   1.315266E-001 +    3.335294E-001i,    7.934257E-002 +    8.734198E-002i,   -9.351067E-002 +    7.153266E-002i,   -7.287501E-002 -    9.310871E-002i,    9.644133E-002 -    6.792056E-002i,   -2.923039E-001 +    1.202725E-002i,    1.224445E-002 -    2.916114E-001i,   -5.856673E-001 -    5.440046E-002i,   -4.931361E-001 +    2.363421E-001i];


% modif JH 25/09/2008 : compatibilite avec ALOHA
S = reshape(S,1,length(S)*length(S)); 

