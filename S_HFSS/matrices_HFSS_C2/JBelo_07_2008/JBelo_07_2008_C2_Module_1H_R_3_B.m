% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,9,9);

f = [3700000000 ];
S(1,:,:) = [   -5.830483E-002 +    4.133540E-003i,   -3.196168E-001 +    1.607799E-001i,   -1.440363E-001 -    3.267173E-001i,    3.315170E-001 -    1.369475E-001i,    1.085001E-001 +    3.409784E-001i,   -2.861382E-001 +    1.980975E-001i,   -1.733144E-001 -    3.008456E-001i,    3.004854E-001 -    1.770507E-001i,    1.613497E-001 +    3.083102E-001i;  -3.196168E-001 +    1.607799E-001i,    4.791219E-001 -    2.984463E-001i,   -5.704930E-001 -    4.045854E-002i,    6.395504E-003 +    2.924575E-001i,   -2.912650E-001 -    1.802506E-002i,   -8.556540E-002 +    8.102280E-002i,   -7.348567E-002 -    9.176296E-002i,    9.145012E-002 -    7.471753E-002i,    6.986566E-002 +    9.487688E-002i;  -1.440363E-001 -    3.267173E-001i,   -5.704930E-001 -    4.045854E-002i,   -5.072253E-001 +    2.505477E-001i,   -2.918164E-001 -    8.460704E-003i,    3.274015E-002 -    2.893891E-001i,   -7.641496E-002 -    8.939285E-002i,    9.518762E-002 -    6.858834E-002i,    6.983200E-002 +    9.493831E-002i,   -9.810764E-002 +    6.482232E-002i;   3.315170E-001 -    1.369475E-001i,    6.395504E-003 +    2.924575E-001i,   -2.918164E-001 -    8.460704E-003i,    4.834785E-001 -    2.499683E-001i,   -5.799729E-001 -    1.100827E-001i,    9.157901E-002 -    7.463420E-002i,    6.663730E-002 +    9.721403E-002i,   -9.699298E-002 +    6.789218E-002i,   -6.278624E-002 -    1.000577E-001i;   1.085001E-001 +    3.409784E-001i,   -2.912650E-001 -    1.802506E-002i,    3.274015E-002 -    2.893891E-001i,   -5.799729E-001 -    1.100827E-001i,   -5.191237E-001 +    1.686927E-001i,    6.655865E-002 +    9.726103E-002i,   -1.021961E-001 +    5.813914E-002i,   -5.940504E-002 -    1.020810E-001i,    1.047019E-001 -    5.407371E-002i;  -2.861382E-001 +    1.980975E-001i,   -8.556540E-002 +    8.102280E-002i,   -7.641496E-002 -    8.939285E-002i,    9.157901E-002 -    7.463420E-002i,    6.655865E-002 +    9.726103E-002i,   -2.087938E-001 -    4.333161E-001i,    6.206495E-002 +    7.481497E-001i,    1.059361E-001 -    5.477565E-002i,    4.926721E-002 +    1.083112E-001i;  -1.733144E-001 -    3.008456E-001i,   -7.348567E-002 -    9.176296E-002i,    9.518762E-002 -    6.858834E-002i,    6.663730E-002 +    9.721403E-002i,   -1.021961E-001 +    5.813914E-002i,    6.206495E-002 +    7.481497E-001i,    1.379718E-001 +    4.618327E-001i,    4.570775E-002 +    1.098476E-001i,   -1.117539E-001 +    4.003499E-002i;   3.004854E-001 -    1.770507E-001i,    9.145012E-002 -    7.471753E-002i,    6.983200E-002 +    9.493831E-002i,   -9.699298E-002 +    6.789218E-002i,   -5.940504E-002 -    1.020810E-001i,    1.059361E-001 -    5.477565E-002i,    4.570775E-002 +    1.098476E-001i,   -1.608017E-001 -    4.212002E-001i,    1.356237E-002 +    7.685898E-001i;   1.613497E-001 +    3.083102E-001i,    6.986566E-002 +    9.487688E-002i,   -9.810764E-002 +    6.482232E-002i,   -6.278624E-002 -    1.000577E-001i,    1.047019E-001 -    5.407371E-002i,    4.926721E-002 +    1.083112E-001i,   -1.117539E-001 +    4.003499E-002i,    1.356237E-002 +    7.685898E-001i,    1.211568E-001 +    4.353362E-001i];


% modif JH 25/09/2008 : compatibilite avec ALOHA
S = reshape(S,1,length(S)*length(S)); 
