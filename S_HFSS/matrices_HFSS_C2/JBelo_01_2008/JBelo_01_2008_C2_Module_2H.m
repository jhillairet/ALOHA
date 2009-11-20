% Matlab m-File exported from HFSS11.0 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,10,10);

f = [3700000000 ];
S(1,:,:) = [   -5.550888E-002 +    3.574086E-002i,    3.535767E-001 +    8.389274E-003i,   -1.009158E-002 -    3.542961E-001i,   -3.525658E-001 +    7.042659E-003i,    2.528608E-002 +    3.524328E-001i,    2.635805E-001 -    2.326183E-001i,   -2.459691E-001 -    2.522252E-001i,   -2.542952E-001 +    2.415101E-001i,    2.547651E-001 +    2.421957E-001i,    3.525854E-002 +    1.488981E-002i;   3.535767E-001 +    8.389274E-003i,    3.936561E-001 +    1.346217E-001i,    3.434891E-001 +    5.745560E-001i,   -9.264502E-002 +    2.358337E-001i,    2.407811E-001 +    8.050573E-002i,    1.001481E-002 -    4.063149E-003i,   -4.566064E-003 -    9.821399E-003i,   -9.836029E-003 +    4.415639E-003i,    4.923356E-003 +    9.620053E-003i,   -6.389134E-002 +    3.471586E-001i;  -1.009158E-002 -    3.542961E-001i,    3.434891E-001 +    5.745560E-001i,   -4.031436E-001 -    9.447614E-002i,    2.408652E-001 +    8.038753E-002i,    6.797991E-002 -    2.451818E-001i,   -4.590022E-003 -    9.810314E-003i,   -9.590453E-003 +    5.083224E-003i,    4.933444E-003 +    9.612961E-003i,    9.370265E-003 -    5.430272E-003i,    3.507750E-001 +    4.579015E-002i;  -3.525658E-001 +    7.042659E-003i,   -9.264502E-002 +    2.358337E-001i,    2.408652E-001 +    8.038753E-002i,    4.130736E-001 +    1.445540E-001i,    3.508621E-001 +    5.549821E-001i,   -9.798868E-003 +    4.483508E-003i,    4.976042E-003 +    9.584306E-003i,    9.605432E-003 -    4.826840E-003i,   -5.323175E-003 -    9.368183E-003i,    4.852291E-002 -    3.485900E-001i;   2.528608E-002 +    3.524328E-001i,    2.407811E-001 +    8.050573E-002i,    6.797991E-002 -    2.451818E-001i,    3.508621E-001 +    5.549821E-001i,   -4.236438E-001 -    1.026961E-001i,    4.993360E-003 +    9.573356E-003i,    9.333175E-003 -    5.475118E-003i,   -5.326911E-003 -    9.362032E-003i,   -9.098952E-003 +    5.811298E-003i,   -3.513288E-001 -    3.052161E-002i;   2.635805E-001 -    2.326183E-001i,    1.001481E-002 -    4.063149E-003i,   -4.590022E-003 -    9.810314E-003i,   -9.798868E-003 +    4.483508E-003i,    4.993360E-003 +    9.573356E-003i,    4.098638E-001 +    1.359318E-001i,    3.416383E-001 +    5.593199E-001i,   -1.060075E-001 +    2.368292E-001i,    2.424925E-001 +    9.399752E-002i,   -1.875180E-001 -    2.996651E-001i;  -2.459691E-001 -    2.522252E-001i,   -4.566064E-003 -    9.821399E-003i,   -9.590453E-003 +    5.083224E-003i,    4.976042E-003 +    9.584306E-003i,    9.333175E-003 -    5.475118E-003i,    3.416383E-001 +    5.593199E-001i,   -4.190131E-001 -    9.612255E-002i,    2.423292E-001 +    9.429878E-002i,    8.199555E-002 -    2.473990E-001i,   -2.905909E-001 +    2.026264E-001i;  -2.542952E-001 +    2.415101E-001i,   -9.836029E-003 +    4.415639E-003i,    4.933444E-003 +    9.612961E-003i,    9.605432E-003 -    4.826840E-003i,   -5.326911E-003 -    9.362032E-003i,   -1.060075E-001 +    2.368292E-001i,    2.423292E-001 +    9.429878E-002i,    4.260327E-001 +    1.419999E-001i,    3.486215E-001 +    5.422484E-001i,    1.978630E-001 +    2.919127E-001i;   2.547651E-001 +    2.421957E-001i,    4.923356E-003 +    9.620053E-003i,    9.370265E-003 -    5.430272E-003i,   -5.323175E-003 -    9.368183E-003i,   -9.098952E-003 +    5.811298E-003i,    2.424925E-001 +    9.399752E-002i,    8.199555E-002 -    2.473990E-001i,    3.486215E-001 +    5.422484E-001i,   -4.358654E-001 -    9.956249E-002i,    2.820844E-001 -    2.129986E-001i;   3.525854E-002 +    1.488981E-002i,   -6.389134E-002 +    3.471586E-001i,    3.507750E-001 +    4.579015E-002i,    4.852291E-002 -    3.485900E-001i,   -3.513288E-001 -    3.052161E-002i,   -1.875180E-001 -    2.996651E-001i,   -2.905909E-001 +    2.026264E-001i,    1.978630E-001 +    2.919127E-001i,    2.820844E-001 -    2.129986E-001i,   -2.229203E-002 -    1.728259E-002i];

% MODIF JH 25/02/2008
% on ne recupere que les 9 coefficients (on ne s'interesse pas au 10 : voie retour vers la charge
S = S(1,1:end-1, 1:end-1);
S = reshape(S,1,length(S)*length(S)); 

S=conj(S);
