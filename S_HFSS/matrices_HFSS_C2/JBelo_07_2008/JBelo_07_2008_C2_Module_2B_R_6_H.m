% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,9,9);

f = [3700000000 ];
S(1,:,:) = [   -5.721485E-002 +    1.448541E-004i,   -2.869517E-001 +    1.976682E-001i,   -1.923836E-001 -    2.894913E-001i,    2.892789E-001 -    1.959319E-001i,    1.806690E-001 +    2.979295E-001i,   -3.287569E-001 +    1.396790E-001i,   -1.224428E-001 -    3.347259E-001i,    3.313450E-001 -    1.370148E-001i,    1.306538E-001 +    3.328106E-001i;  -2.869517E-001 +    1.976682E-001i,   -2.089912E-001 -    4.330451E-001i,    1.108820E-001 +    7.420292E-001i,    1.026629E-001 -    6.188238E-002i,    5.649627E-002 +    1.053500E-001i,   -9.102288E-002 +    7.517759E-002i,   -7.027805E-002 -    9.453529E-002i,    9.209347E-002 -    7.457692E-002i,    7.269940E-002 +    9.315615E-002i;  -1.923836E-001 -    2.894913E-001i,    1.108820E-001 +    7.420292E-001i,    1.981119E-001 +    4.392836E-001i,    6.001468E-002 +    1.034246E-001i,   -1.060151E-001 +    5.459792E-002i,   -7.346881E-002 -    9.203584E-002i,    9.545767E-002 -    6.852363E-002i,    7.285191E-002 +    9.309367E-002i,   -9.412236E-002 +    7.096161E-002i;   2.892789E-001 -    1.959319E-001i,    1.026629E-001 -    6.188238E-002i,    6.001468E-002 +    1.034246E-001i,   -2.143668E-001 -    3.959097E-001i,    1.134536E-001 +    7.599805E-001i,    9.185952E-002 -    7.466016E-002i,    6.971981E-002 +    9.534265E-002i,   -9.292824E-002 +    7.404944E-002i,   -7.215851E-002 -    9.397893E-002i;   1.806690E-001 +    2.979295E-001i,    5.649627E-002 +    1.053500E-001i,   -1.060151E-001 +    5.459792E-002i,    1.134536E-001 +    7.599805E-001i,    1.777922E-001 +    4.149940E-001i,    6.976372E-002 +    9.522930E-002i,   -9.845123E-002 +    6.466878E-002i,   -6.910195E-002 -    9.626318E-002i,    9.721503E-002 -    6.716599E-002i;  -3.287569E-001 +    1.396790E-001i,   -9.102288E-002 +    7.517759E-002i,   -7.346881E-002 -    9.203584E-002i,    9.185952E-002 -    7.466016E-002i,    6.976372E-002 +    9.522930E-002i,    5.132498E-001 -    2.350327E-001i,   -5.602303E-001 -    1.132400E-001i,   -1.161975E-002 +    2.928532E-001i,   -2.917654E-001 -    1.684340E-002i;  -1.224428E-001 -    3.347259E-001i,   -7.027805E-002 -    9.453529E-002i,    9.545767E-002 -    6.852363E-002i,    6.971981E-002 +    9.534265E-002i,   -9.845123E-002 +    6.466878E-002i,   -5.602303E-001 -    1.132400E-001i,   -5.351635E-001 +    1.838903E-001i,   -2.912394E-001 -    2.649888E-002i,    3.164889E-002 -    2.898894E-001i;   3.313450E-001 -    1.370148E-001i,    9.209347E-002 -    7.457692E-002i,    7.285191E-002 +    9.309367E-002i,   -9.292824E-002 +    7.404944E-002i,   -6.910195E-002 -    9.626318E-002i,   -1.161975E-002 +    2.928532E-001i,   -2.912394E-001 -    2.649888E-002i,    4.837596E-001 -    2.509050E-001i,   -5.845622E-001 -    7.176261E-002i;   1.306538E-001 +    3.328106E-001i,    7.269940E-002 +    9.315615E-002i,   -9.412236E-002 +    7.096161E-002i,   -7.215851E-002 -    9.397893E-002i,    9.721503E-002 -    6.716599E-002i,   -2.917654E-001 -    1.684340E-002i,    3.164889E-002 -    2.898894E-001i,   -5.845622E-001 -    7.176261E-002i,   -4.929487E-001 +    2.366399E-001i];


% modif JH 25/09/2008 : compatibilite avec ALOHA
S = reshape(S,1,length(S)*length(S)); 
