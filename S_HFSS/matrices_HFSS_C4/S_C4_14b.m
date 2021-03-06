% Matlab m-File from Maxwell project
f = zeros(1,1);
S = zeros(1,49);

f(1) = 3700000000;

S = [
0.1085257984222925 + -0.07658868701540085i, 0.1219272727723845 + 0.3852421286088176i, 0.3910874450929067 + -0.0818314184920374i, -0.07771003872730925 + -0.4041929813752048i, -0.4000621893758501 + 0.08777592075993992i, 0.07017924876192062 + 0.3876059765243775i, 0.4015048008079927 + -0.07679835170166496i, 0.121927272772375 + 0.3852421286088206i, -0.3226779451179118 + -0.05299119739304809i, -0.2565278220322628 + -0.6614379091344479i, -0.1593729809385723 + -0.1022036205802037i, -0.09753877672129963 + 0.1611863182331602i, -0.2649569676212863 + 0.0670577853760852i, 0.07230218655823666 + 0.2742622310878199i, 0.391087445092935 + -0.0818314184920583i, -0.2565278220322655 + -0.6614379091344899i, 0.321041246321038 + 0.1119448079052039i, -0.1163259181580413 + 0.1466842096066962i, 0.1489299333992269 + 0.1119158738247047i, 0.03975247727063165 + 0.2673145394585049i, 0.2769885534297012 + -0.04399125617160434i, -0.07771003872729961 + -0.4041929813751923i, -0.1593729809385766 + -0.1022036205802026i, -0.1163259181580257 + 0.1466842096066869i, -0.7173171443170512 + -0.1962212642180561i, -0.2469171757901005 + -0.2739682135262255i, -0.1412588884732997 + -0.1185278439695894i, -0.1215502940527183 + 0.1478010110023652i, -0.4000621893758609 + 0.08777592075994228i, -0.09753877672130284 + 0.1611863182331781i, 0.1489299333992307 + 0.1119158738246897i, -0.2469171757901044 + -0.2739682135262394i, 0.7291874527964528 + 0.1561168710481969i, -0.1142468396574183 + 0.1435900847583034i, 0.1501762519532742 + 0.1170839460068607i, 0.0701792487619344 + 0.387605976524376i, -0.2649569676212814 + 0.06705778537608868i, 0.0397524772706276 + 0.2673145394584944i, -0.1412588884733109 + -0.1185278439695946i, -0.1142468396574183 + 0.1435900847583085i, -0.3316942509826617 + -0.1274413371953419i, -0.1563762939341991 + -0.6918893872817427i, 0.401504800807935 + -0.07679835170164738i, 0.07230218655822071 + 0.2742622310877869i, 0.2769885534296615 + -0.04399125617158643i, -0.1215502940527114 + 0.1478010110023515i, 0.1501762519532528 + 0.1170839460068475i, -0.1563762939341917 + -0.6918893872816542i, 0.2829569285029747 + 0.1329675838116993i;
];
Z = [
-4.352307816033595e-013 + -1091.043744958272i, -1.838678297013057e-013 + 18.10277382954064i, 5.528548444340003e-013 + -472.5645613709812i, 9.247622326931466e-014 + -25.80690364256305i, -5.941491363843968e-015 + 530.8329329370428i, -1.947841776693381e-013 + 7.737870718482471i, 7.16208832660224e-013 + -441.1782957995206i, 2.363389300053404e-013 + 18.10277382955744i, 1.289291378965546e-013 + 18.64945464871284i, -3.94928152892272e-013 + -105.5911672014101i, 3.952011978747955e-015 + 0.8137272706315891i, -2.495246369655522e-013 + -16.73789454248384i, 3.012111646200212e-014 + -0.06561763407650832i, 1.536232664136143e-013 + 3.741212132138972i, 5.865336905636554e-014 + -472.5645613710361i, -4.157630238343214e-014 + -105.5911672014238i, -5.938791498739486e-015 + -48.51292889191091i, 3.266596786008051e-013 + -21.24197841752578i, -1.068665143304331e-012 + 436.9350856406185i, 1.083424972961315e-014 + 1.712917645773765i, 5.754523225731104e-015 + -97.66259411486446i, 3.207098473617504e-013 + -25.80690364256426i, -6.151892739903853e-014 + 0.8137272706316169i, 2.378514421481748e-013 + -21.24197841752651i, -1.712853394196685e-013 + 7.720979349700905i, 7.74465847732354e-014 + -111.4994379118382i, 1.244516558577985e-014 + 0.3477073674535581i, 3.161786313199427e-013 + -19.82471588032013i, 2.574901347561629e-013 + 530.8329329370224i, 1.211948127720128e-013 + -16.73789454247651i, 1.646619881016179e-013 + 436.9350856405836i, 1.831779968309413e-013 + -111.4994379118433i, 5.372537711788717e-013 + 167.3187469286699i, -6.855509946940158e-014 + -7.152137414340522i, 1.196646533802716e-013 + 407.7828239213718i, -7.922415072624557e-014 + 7.737870718472289i, 2.298715155566745e-014 + -0.06561763407605521i, -4.500447802829035e-014 + 1.712917645769445i, 5.096290073702075e-014 + 0.3477073674532293i, 1.29202346495546e-013 + -7.152137414332427i, -2.305464329334255e-014 + 2.476778338023018i, 6.654842202802823e-014 + -94.278164942852i, -1.041717990253853e-012 + -441.1782957994764i, -1.556842672888237e-013 + 3.741212132129943i, -3.173888523642791e-013 + -97.66259411483419i, 2.773545137022876e-013 + -19.82471588031671i, -2.469630457546589e-013 + 407.7828239213319i, -7.250087653756908e-014 + -94.27816494283688i, -3.649361892561394e-013 + -55.72278133611158i;
];
