% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,9,9);

f = [3700000000 ];
S(1,:,:) = [   -5.719674E-002 +    3.222848E-002i,   -6.838087E-002 +    3.537330E-001i,   -3.460874E-001 -    9.686129E-002i,    9.367597E-002 -    3.489715E-001i,    3.430940E-001 +    1.107692E-001i,   -1.358389E-001 +    3.167629E-001i,   -3.087318E-001 -    1.516879E-001i,    1.605876E-001 -    3.065947E-001i,    2.913196E-001 +    1.855193E-001i;  -6.838087E-002 +    3.537330E-001i,   -3.836204E-001 -    4.125511E-001i,    5.422447E-002 +    5.715808E-001i,    2.837353E-001 -    6.328743E-002i,    4.890904E-002 +    2.859160E-001i,    6.929654E-002 +    9.382516E-002i,   -9.708131E-002 +    6.424151E-002i,   -6.209591E-002 -    9.931534E-002i,    1.039860E-001 -    5.337243E-002i;  -3.460874E-001 -    9.686129E-002i,    5.422447E-002 +    5.715808E-001i,    3.140346E-001 +    4.695718E-001i,    3.975384E-002 +    2.872474E-001i,   -2.882416E-001 +    2.528127E-002i,   -9.893544E-002 +    6.123241E-002i,   -5.594103E-002 -    1.017599E-001i,    1.038057E-001 -    5.362555E-002i,    4.457162E-002 +    1.077369E-001i;   9.367597E-002 -    3.489715E-001i,    2.837353E-001 -    6.328743E-002i,    3.975384E-002 +    2.872474E-001i,   -3.389086E-001 -    4.270450E-001i,    2.538497E-002 +    5.900134E-001i,   -6.261778E-002 -    9.880862E-002i,    1.017048E-001 -    5.732836E-002i,    5.502240E-002 +    1.037863E-001i,   -1.078353E-001 +    4.596219E-002i;   3.430940E-001 +    1.107692E-001i,    4.890904E-002 +    2.859160E-001i,   -2.882416E-001 +    2.528127E-002i,    2.538497E-002 +    5.900134E-001i,    2.971452E-001 +    4.587864E-001i,    1.015932E-001 -    5.747079E-002i,    5.205506E-002 +    1.042154E-001i,   -1.061747E-001 +    4.965317E-002i,   -4.042210E-002 -    1.097576E-001i;  -1.358389E-001 +    3.167629E-001i,    6.929654E-002 +    9.382516E-002i,   -9.893544E-002 +    6.123241E-002i,   -6.261778E-002 -    9.880862E-002i,    1.015932E-001 -    5.747079E-002i,    3.222503E-001 -    3.481952E-001i,   -5.138016E-001 -    1.931361E-001i,    1.587932E-001 +    3.539647E-001i,   -3.653413E-001 +    1.280569E-001i;  -3.087318E-001 -    1.516879E-001i,   -9.708131E-002 +    6.424151E-002i,   -5.594103E-002 -    1.017599E-001i,    1.017048E-001 -    5.732836E-002i,    5.205506E-002 +    1.042154E-001i,   -5.138016E-001 -    1.931361E-001i,   -3.572500E-001 +    3.151463E-001i,   -3.609688E-001 +    1.400604E-001i,   -1.088401E-001 -    3.707267E-001i;   1.605876E-001 -    3.065947E-001i,   -6.209591E-002 -    9.931534E-002i,    1.038057E-001 -    5.362555E-002i,    5.502240E-002 +    1.037863E-001i,   -1.061747E-001 +    4.965317E-002i,    1.587932E-001 +    3.539647E-001i,   -3.609688E-001 +    1.400604E-001i,    3.344438E-001 -    3.234016E-001i,   -4.948033E-001 -    2.522098E-001i;   2.913196E-001 +    1.855193E-001i,    1.039860E-001 -    5.337243E-002i,    4.457162E-002 +    1.077369E-001i,   -1.078353E-001 +    4.596219E-002i,   -4.042210E-002 -    1.097576E-001i,   -3.653413E-001 +    1.280569E-001i,   -1.088401E-001 -    3.707267E-001i,   -4.948033E-001 -    2.522098E-001i,   -3.853462E-001 +    2.644792E-001i];


% modif JH 25/09/2008 : compatibilite avec ALOHA
S = reshape(S,1,length(S)*length(S)); 
