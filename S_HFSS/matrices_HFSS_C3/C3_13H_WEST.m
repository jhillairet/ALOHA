% Matlab m-File exported from HFSS13.0.2 
% Note: In three-dimensional arrays, like S(i,j,k), the first index corresponds to the frequency.
%       So, S(N,j,k) is an S(j,k) matrix for frequency N.

f = zeros(1,1);
S = zeros(1,19,19);
Zo = zeros(1,19);

f = [3700000000 ];
S(1,:,:) = [    7.223876E-002 +    5.700207E-002i,    2.269091E-001 -    4.340610E-002i,   -3.897936E-002 -    2.245644E-001i,   -2.259779E-001 +    4.192887E-002i,    3.812561E-002 +    2.247725E-001i,    2.187450E-001 -    7.059013E-002i,   -6.762968E-002 -    2.187267E-001i,   -2.377221E-001 +    4.980048E-002i,    5.204069E-002 +    2.350507E-001i,    2.380870E-001 -    4.585633E-002i,   -4.952600E-002 -    2.355813E-001i,   -2.299741E-001 +    7.723718E-002i,    8.121103E-002 +    2.272277E-001i,    2.159302E-001 -    6.816552E-002i,   -8.566223E-002 -    2.106773E-001i,   -2.171748E-001 +    6.324650E-002i,    8.260536E-002 +    2.117365E-001i,    2.068770E-001 -    9.226520E-002i,   -1.118970E-001 -    1.980152E-001i;   2.269091E-001 -    4.340610E-002i,    3.808416E-001 +    5.395265E-001i,    7.287032E-002 +    4.936229E-001i,    2.541542E-001 -    2.857343E-002i,   -2.444788E-002 -    2.525139E-001i,    1.584039E-002 -    1.804751E-001i,   -1.795538E-001 -    1.798427E-002i,   -1.199355E-002 +    5.023207E-002i,    4.992163E-002 +    1.132037E-002i,    1.278511E-002 -    4.994436E-002i,   -4.979549E-002 -    1.185208E-002i,   -6.014520E-003 +    5.123165E-002i,    5.106410E-002 +    5.001164E-003i,    6.962262E-002 +    2.736960E-003i,   -2.876373E-003 -    6.992271E-002i,   -6.947102E-002 -    4.290923E-003i,    1.880930E-003 +    6.991127E-002i,    6.950958E-002 -    5.183960E-003i,   -1.179119E-002 -    6.898687E-002i;  -3.897936E-002 -    2.245644E-001i,    7.287032E-002 +    4.936229E-001i,   -3.758400E-001 -    5.515139E-001i,   -2.388715E-002 -    2.511890E-001i,   -2.495011E-001 +    1.984532E-002i,   -1.777583E-001 -    1.867788E-002i,   -2.077708E-002 +    1.768133E-001i,    4.934763E-002 +    1.268041E-002i,    1.201111E-002 -    4.905280E-002i,   -4.905044E-002 -    1.345635E-002i,   -1.253347E-002 +    4.891937E-002i,    5.043478E-002 +    6.799405E-003i,    5.796965E-003 -    5.028664E-002i,    3.877368E-003 -    6.863165E-002i,   -6.902260E-002 +    1.654723E-003i,   -5.407682E-003 +    6.845582E-002i,    6.899448E-002 -    6.729804E-004i,   -3.937986E-003 -    6.865410E-002i,   -6.825024E-002 +    1.046441E-002i;  -2.259779E-001 +    4.192887E-002i,    2.541542E-001 -    2.857343E-002i,   -2.388715E-002 -    2.511890E-001i,    6.974714E-001 +    2.845821E-001i,   -7.475576E-002 +    2.247925E-001i,    2.475588E-001 -    5.936463E-002i,   -5.609696E-002 -    2.472870E-001i,    1.220910E-002 -    4.990660E-002i,   -4.960149E-002 -    1.153767E-002i,   -1.299497E-002 +    4.961601E-002i,    4.947306E-002 +    1.206594E-002i,    6.266436E-003 -    5.093405E-002i,   -5.077297E-002 -    5.257384E-003i,   -6.924823E-002 -    3.107362E-003i,    2.475347E-003 +    6.957778E-002i,    6.908883E-002 +    4.652469E-003i,   -1.485104E-003 -    6.956089E-002i,   -6.917952E-002 +    4.773307E-003i,    1.134933E-002 +    6.869600E-002i;   3.812561E-002 +    2.247725E-001i,   -2.444788E-002 -    2.525139E-001i,   -2.495011E-001 +    1.984532E-002i,   -7.475576E-002 +    2.247925E-001i,   -6.923571E-001 -    3.081714E-001i,   -5.508834E-002 -    2.464439E-001i,   -2.461243E-001 +    5.185155E-002i,   -4.931162E-002 -    1.287351E-002i,   -1.220291E-002 +    4.901928E-002i,    4.901137E-002 +    1.364851E-002i,    1.272489E-002 -    4.888382E-002i,   -5.042167E-002 -    6.995171E-003i,   -5.991902E-003 +    5.027736E-002i,   -4.142375E-003 +    6.863450E-002i,    6.904684E-002 -    1.389645E-003i,    5.672409E-003 -    6.845274E-002i,   -6.901493E-002 +    4.077559E-004i,    3.674917E-003 +    6.868703E-002i,    6.830816E-002 -    1.020458E-002i;   2.187450E-001 -    7.059013E-002i,    1.584039E-002 -    1.804751E-001i,   -1.777583E-001 -    1.867788E-002i,    2.475588E-001 -    5.936463E-002i,   -5.508834E-002 -    2.464439E-001i,    5.046605E-001 +    4.316324E-001i,    1.901770E-001 +    4.591763E-001i,   -5.703458E-003 +    5.106482E-002i,    5.067603E-002 +    5.076699E-003i,    6.520196E-003 -    5.087747E-002i,   -5.061645E-002 -    5.617120E-003i,    3.223925E-004 +    5.132110E-002i,    5.103182E-002 -    1.302502E-003i,    6.907976E-002 -    5.806251E-003i,   -1.138552E-002 -    6.869004E-002i,   -6.911998E-002 +    4.253343E-003i,    1.040122E-002 +    6.880040E-002i,    6.800010E-002 -    1.361354E-002i,   -2.007363E-002 -    6.667650E-002i;  -6.762968E-002 -    2.187267E-001i,   -1.795538E-001 -    1.798427E-002i,   -2.077708E-002 +    1.768133E-001i,   -5.609696E-002 -    2.472870E-001i,   -2.461243E-001 +    5.185155E-002i,    1.901770E-001 +    4.591763E-001i,   -4.979917E-001 -    4.432877E-001i,    5.078920E-002 +    6.305136E-003i,    5.676147E-003 -    5.040964E-002i,   -5.059261E-002 -    7.116289E-003i,   -6.213661E-003 +    5.034369E-002i,    5.111815E-002 +    3.067116E-004i,   -6.729860E-004 -    5.084203E-002i,   -4.937801E-003 -    6.887231E-002i,   -6.855242E-002 +    1.049935E-002i,    3.390661E-003 +    6.889337E-002i,    6.865029E-002 -    9.517667E-003i,   -1.272682E-002 -    6.789251E-002i,   -6.665327E-002 +    1.917706E-002i;  -2.377221E-001 +    4.980048E-002i,   -1.199355E-002 +    5.023207E-002i,    4.934763E-002 +    1.268041E-002i,    1.220910E-002 -    4.990660E-002i,   -4.931162E-002 -    1.287351E-002i,   -5.703458E-003 +    5.106482E-002i,    5.078920E-002 +    6.305136E-003i,    4.982310E-001 +    5.272445E-001i,    1.874124E-001 +    4.269088E-001i,    1.958616E-001 -    1.029108E-001i,   -1.054547E-001 -    1.927041E-001i,    5.609209E-002 -    1.475816E-001i,   -1.478406E-001 -    5.295909E-002i,   -6.032162E-003 +    5.019487E-002i,    5.073865E-002 +    1.989817E-003i,    7.147231E-003 -    4.999476E-002i,   -5.067242E-002 -    2.709100E-003i,   -2.968164E-004 +    5.057417E-002i,    5.058072E-002 -    4.512253E-003i;   5.204069E-002 +    2.350507E-001i,    4.992163E-002 +    1.132037E-002i,    1.201111E-002 -    4.905280E-002i,   -4.960149E-002 -    1.153767E-002i,   -1.220291E-002 +    4.901928E-002i,    5.067603E-002 +    5.076699E-003i,    5.676147E-003 -    5.040964E-002i,    1.874124E-001 +    4.269088E-001i,   -5.181612E-001 -    5.127396E-001i,   -1.042074E-001 -    1.929625E-001i,   -1.898043E-001 +    1.066932E-001i,   -1.469048E-001 -    5.392922E-002i,   -5.082110E-002 +    1.471262E-001i,    4.981750E-002 +    5.412301E-003i,    1.399693E-003 -    5.031085E-002i,   -4.963175E-002 -    6.519733E-003i,   -2.113339E-003 +    5.025332E-002i,    5.012873E-002 -    2.764290E-004i,   -5.042897E-003 -    5.008096E-002i;   2.380870E-001 -    4.585633E-002i,    1.278511E-002 -    4.994436E-002i,   -4.905044E-002 -    1.345635E-002i,   -1.299497E-002 +    4.961601E-002i,    4.901137E-002 +    1.364851E-002i,    6.520196E-003 -    5.087747E-002i,   -5.059261E-002 -    7.116289E-003i,    1.958616E-001 -    1.029108E-001i,   -1.042074E-001 -    1.929625E-001i,    7.777389E-001 +    2.418881E-001i,   -6.542378E-004 +    1.823493E-001i,    1.822618E-001 -    1.250487E-001i,   -1.278460E-001 -    1.788639E-001i,    6.834207E-003 -    5.000381E-002i,   -5.061208E-002 -    2.808188E-003i,   -7.943961E-003 +    4.978601E-002i,    5.053432E-002 +    3.525060E-003i,    1.115670E-003 -    5.047534E-002i,   -5.055979E-002 +    3.684353E-003i;  -4.952600E-002 -    2.355813E-001i,   -4.979549E-002 -    1.185208E-002i,   -1.253347E-002 +    4.891937E-002i,    4.947306E-002 +    1.206594E-002i,    1.272489E-002 -    4.888382E-002i,   -5.061645E-002 -    5.617120E-003i,   -6.213661E-003 +    5.034369E-002i,   -1.054547E-001 -    1.927041E-001i,   -1.898043E-001 +    1.066932E-001i,   -6.542378E-004 +    1.823493E-001i,   -7.886527E-001 -    2.156817E-001i,   -1.272027E-001 -    1.788311E-001i,   -1.754109E-001 +    1.299224E-001i,   -4.975444E-002 -    5.943535E-003i,   -1.936619E-003 +    5.029057E-002i,    4.955689E-002 +    7.048867E-003i,    2.649576E-003 -    5.022544E-002i,   -5.012636E-002 -    2.587308E-004i,    4.507742E-003 +    5.012949E-002i;  -2.299741E-001 +    7.723718E-002i,   -6.014520E-003 +    5.123165E-002i,    5.043478E-002 +    6.799405E-003i,    6.266436E-003 -    5.093405E-002i,   -5.042167E-002 -    6.995171E-003i,    3.223925E-004 +    5.132110E-002i,    5.111815E-002 +    3.067116E-004i,    5.609209E-002 -    1.475816E-001i,   -1.469048E-001 -    5.392922E-002i,    1.822618E-001 -    1.250487E-001i,   -1.272027E-001 -    1.788311E-001i,    6.088562E-001 +    3.942644E-001i,    2.856522E-001 +    3.693683E-001i,   -1.055340E-004 +    5.049665E-002i,    5.056267E-002 -    3.967905E-003i,    1.235048E-003 -    5.042873E-002i,   -5.058121E-002 +    3.246666E-003i,    5.627998E-003 +    5.020126E-002i,    4.964461E-002 -    1.039907E-002i;   8.121103E-002 +    2.272277E-001i,    5.106410E-002 +    5.001164E-003i,    5.796965E-003 -    5.028664E-002i,   -5.077297E-002 -    5.257384E-003i,   -5.991902E-003 +    5.027736E-002i,    5.103182E-002 -    1.302502E-003i,   -6.729860E-004 -    5.084203E-002i,   -1.478406E-001 -    5.295909E-002i,   -5.082110E-002 +    1.471262E-001i,   -1.278460E-001 -    1.788639E-001i,   -1.754109E-001 +    1.299224E-001i,    2.856522E-001 +    3.693683E-001i,   -6.280116E-001 -    3.676527E-001i,    5.022010E-002 -    8.611629E-004i,   -4.913395E-003 -    5.020783E-002i,   -5.017417E-002 -    2.634196E-004i,    4.196488E-003 +    5.024007E-002i,    4.981665E-002 -    6.557415E-003i,   -1.129152E-002 -    4.917178E-002i;   2.159302E-001 -    6.816552E-002i,    6.962262E-002 +    2.736960E-003i,    3.877368E-003 -    6.863165E-002i,   -6.924823E-002 -    3.107362E-003i,   -4.142375E-003 +    6.863450E-002i,    6.907976E-002 -    5.806251E-003i,   -4.937801E-003 -    6.887231E-002i,   -6.032162E-003 +    5.019487E-002i,    4.981750E-002 +    5.412301E-003i,    6.834207E-003 -    5.000381E-002i,   -4.975444E-002 -    5.943535E-003i,   -1.055340E-004 +    5.049665E-002i,    5.022010E-002 -    8.611629E-004i,    4.968207E-001 +    4.445266E-001i,    2.203137E-001 +    4.466822E-001i,    2.403647E-001 -    8.129863E-002i,   -1.027342E-001 -    2.333413E-001i,   -2.252129E-002 -    1.777128E-001i,   -1.754952E-001 +    3.940980E-002i;  -8.566223E-002 -    2.106773E-001i,   -2.876373E-003 -    6.992271E-002i,   -6.902260E-002 +    1.654723E-003i,    2.475347E-003 +    6.957778E-002i,    6.904684E-002 -    1.389645E-003i,   -1.138552E-002 -    6.869004E-002i,   -6.855242E-002 +    1.049935E-002i,    5.073865E-002 +    1.989817E-003i,    1.399693E-003 -    5.031085E-002i,   -5.061208E-002 -    2.808188E-003i,   -1.936619E-003 +    5.029057E-002i,    5.056267E-002 -    3.967905E-003i,   -4.913395E-003 -    5.020783E-002i,    2.203137E-001 +    4.466822E-001i,   -5.584553E-001 -    3.615834E-001i,   -1.007814E-001 -    2.340801E-001i,   -2.253196E-001 +    1.216748E-001i,   -1.760984E-001 +    3.688301E-002i,    5.361187E-002 +    1.725159E-001i;  -2.171748E-001 +    6.324650E-002i,   -6.947102E-002 -    4.290923E-003i,   -5.407682E-003 +    6.845582E-002i,    6.908883E-002 +    4.652469E-003i,    5.672409E-003 -    6.845274E-002i,   -6.911998E-002 +    4.253343E-003i,    3.390661E-003 +    6.889337E-002i,    7.147231E-003 -    4.999476E-002i,   -4.963175E-002 -    6.519733E-003i,   -7.943961E-003 +    4.978601E-002i,    4.955689E-002 +    7.048867E-003i,    1.235048E-003 -    5.042873E-002i,   -5.017417E-002 -    2.634196E-004i,    2.403647E-001 -    8.129863E-002i,   -1.007814E-001 -    2.340801E-001i,    7.421450E-001 +    1.404655E-001i,   -2.195464E-003 +    2.393262E-001i,    2.297761E-001 -    1.082337E-001i,   -1.300360E-001 -    2.193860E-001i;   8.260536E-002 +    2.117365E-001i,    1.880930E-003 +    6.991127E-002i,    6.899448E-002 -    6.729804E-004i,   -1.485104E-003 -    6.956089E-002i,   -6.901493E-002 +    4.077559E-004i,    1.040122E-002 +    6.880040E-002i,    6.865029E-002 -    9.517667E-003i,   -5.067242E-002 -    2.709100E-003i,   -2.113339E-003 +    5.025332E-002i,    5.053432E-002 +    3.525060E-003i,    2.649576E-003 -    5.022544E-002i,   -5.058121E-002 +    3.246666E-003i,    4.196488E-003 +    5.024007E-002i,   -1.027342E-001 -    2.333413E-001i,   -2.253196E-001 +    1.216748E-001i,   -2.195464E-003 +    2.393262E-001i,   -7.536271E-001 -    7.679052E-003i,   -1.287508E-001 -    2.203494E-001i,   -2.080124E-001 +    1.496479E-001i;   2.068770E-001 -    9.226520E-002i,    6.950958E-002 -    5.183960E-003i,   -3.937986E-003 -    6.865410E-002i,   -6.917952E-002 +    4.773307E-003i,    3.674917E-003 +    6.868703E-002i,    6.800010E-002 -    1.361354E-002i,   -1.272682E-002 -    6.789251E-002i,   -2.968164E-004 +    5.057417E-002i,    5.012873E-002 -    2.764290E-004i,    1.115670E-003 -    5.047534E-002i,   -5.012636E-002 -    2.587308E-004i,    5.627998E-003 +    5.020126E-002i,    4.981665E-002 -    6.557415E-003i,   -2.252129E-002 -    1.777128E-001i,   -1.760984E-001 +    3.688301E-002i,    2.297761E-001 -    1.082337E-001i,   -1.287508E-001 -    2.203494E-001i,    5.843656E-001 +    3.158951E-001i,    3.256063E-001 +    3.805471E-001i;  -1.118970E-001 -    1.980152E-001i,   -1.179119E-002 -    6.898687E-002i,   -6.825024E-002 +    1.046441E-002i,    1.134933E-002 +    6.869600E-002i,    6.830816E-002 -    1.020458E-002i,   -2.007363E-002 -    6.667650E-002i,   -6.665327E-002 +    1.917706E-002i,    5.058072E-002 -    4.512253E-003i,   -5.042897E-003 -    5.008096E-002i,   -5.055979E-002 +    3.684353E-003i,    4.507742E-003 +    5.012949E-002i,    4.964461E-002 -    1.039907E-002i,   -1.129152E-002 -    4.917178E-002i,   -1.754952E-001 +    3.940980E-002i,    5.361187E-002 +    1.725159E-001i,   -1.300360E-001 -    2.193860E-001i,   -2.080124E-001 +    1.496479E-001i,    3.256063E-001 +    3.805471E-001i,   -6.312164E-001 -    2.029265E-001i];
Zo(1,:) = [ 
   2.854885E+002 +    1.561245E-002i,    6.104198E+001 -    4.727280E-004i,    6.104056E+001 -    4.729457E-004i,    6.103904E+001 -    4.729899E-004i,    6.103746E+001 -    4.731246E-004i,    6.103593E+001 -    4.730769E-004i,    6.103386E+001 -    4.735228E-004i,    6.341451E+001 -    2.701512E-004i,    6.341365E+001 -    2.705448E-004i,    6.341295E+001 -    2.705410E-004i,    6.341215E+001 -    2.705571E-004i,    6.341156E+001 -    2.706630E-004i,    6.341073E+001 -    2.708408E-004i,    6.477919E+001 -    1.431469E-004i,    6.477911E+001 -    1.432676E-004i,    6.477894E+001 -    1.431259E-004i,    6.477902E+001 -    1.429240E-004i,    6.477868E+001 -    1.432389E-004i,    6.477846E+001 -    1.432041E-004i ];