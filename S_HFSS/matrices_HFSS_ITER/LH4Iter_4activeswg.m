% Matlab m-File exported from HFSS11.1.3 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,5,5);

f = [5000000000 ];
S(1,:,:) = [    1.430606E-002 +    1.246157E-002i,    5.033351E-001 +    3.565191E-002i,    2.980960E-002 -    5.011048E-001i,   -4.969749E-001 -    3.413495E-002i,   -3.552709E-002 +    4.935755E-001i;   5.033351E-001 +    3.565191E-002i,   -1.641503E-001 -    7.019299E-001i,    2.167344E-001 -    2.236054E-001i,    3.605305E-002 -    2.520050E-001i,   -2.502203E-001 -    3.663495E-002i;   2.980960E-002 -    5.011048E-001i,    2.167344E-001 -    2.236054E-001i,    1.823047E-001 +    7.002605E-001i,   -2.510930E-001 -    3.303283E-002i,   -3.363174E-002 +    2.493240E-001i;  -4.969749E-001 -    3.413495E-002i,    3.605305E-002 -    2.520050E-001i,   -2.510930E-001 -    3.303283E-002i,    2.281223E-001 +    2.041120E-001i,   -7.070201E-001 +    1.712503E-001i;  -3.552709E-002 +    4.935755E-001i,   -2.502203E-001 -    3.663495E-002i,   -3.363174E-002 +    2.493240E-001i,   -7.070201E-001 +    1.712503E-001i,   -2.290207E-001 -    2.149899E-001i];
