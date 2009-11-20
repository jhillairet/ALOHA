% Matlab m-File from Maxwell project
f = zeros(1,1);
S = zeros(1,361);

f(1) = 3700000000;

S = [
0.09177861171469806 + 0.0438404281421927i, -0.1632632490607834 + -0.1616214455283667i, 0.1770945396596759 + -0.1441581771762217i, 0.1353593693040233 + 0.1854813786154635i, -0.1959123579154299 + 0.115270247326658i, -0.1348728216699423 + -0.1865348602265849i, 0.1952332425323216 + -0.1178537100199396i, 0.1525968069538569 + 0.1918033025474161i, -0.1961650118983873 + 0.1457842897593855i, -0.1202548022673236 + -0.2128927938006668i, 0.2154972257377076 + -0.114736730614554i, 0.1189288951605413 + 0.2137123831495121i, -0.216137430854771 + 0.1161298649238381i, -0.1448292478234278 + -0.1821851739600214i, 0.1779630932908142 + -0.1427140458171142i, 0.1130243007702981 + 0.2002863942089489i, -0.1963544656890609 + 0.1139279967111884i, -0.112065581024801 + -0.2015280958573793i, 0.197520602193427 + -0.1154039852576388i, -0.1632632490607983 + -0.1616214455283645i, -0.7080491308035075 + -0.05231528047680007i, 0.3956778284304666 + 0.2628079948706081i, -0.1378160654310065 + 0.2069087540927264i, -0.1904000292314345 + -0.1559287253902001i, 0.06560461647811046 + 0.1491213238028084i, -0.1526109204735072 + 0.05253515498147469i, -0.0362101285818425 + -0.03735140544292893i, 0.03841398593767856 + -0.03486153122936208i, 0.02981566718969075 + 0.04247674792409812i, -0.04314129552047825 + 0.02870406590744952i, -0.02955257758696737 + -0.04267736000877357i, 0.04324764199659999 + -0.02901160838497849i, -0.03935537304529548 + 0.05450865582785714i, -0.05359061798842866 + -0.03834638056057324i, 0.04687371181058612 + -0.04707719629200345i, 0.047021409549382 + 0.04570961469856551i, -0.04729391562069667 + 0.04690755733334411i, -0.04752298655269371 + -0.04591863332607602i, 0.1770945396596534 + -0.1441581771762209i, 0.395677828430416 + 0.2628079948705884i, 0.6843585527051753 + 0.1973888428977537i, -0.1898245423923793 + -0.1582107825620122i, 0.1743551383804714 + -0.1715825512913321i, -0.1543526472273038 + 0.04897346920096429i, -0.03568594191997022 + -0.1564111202539382i, 0.04076554172092805 + -0.03181268282776435i, 0.03036685699557779 + 0.04167221791250211i, -0.04515063395283044 + 0.02494794672353174i, -0.02377868003346392 + -0.04568915089366146i, 0.04532090989257586 + -0.02466660000838227i, 0.02407130670199363 + 0.04582696551263938i, -0.04968364134528652 + -0.04469249669008455i, 0.04359767205298627 + -0.04888368946935429i, 0.04153950746586189 + 0.05133223835087507i, -0.05017584657756787 + 0.04160820276079882i, -0.04132715783867426 + -0.05172947371212822i, 0.05043577037151029 + -0.0420816689489407i, 0.1353593693040245 + 0.1854813786154713i, -0.1378160654309903 + 0.2069087540927335i, -0.1898245423923678 + -0.1582107825620307i, -0.7728939320567163 + 0.1349284578552361i, 0.04518069085268198 + 0.2100360480209975i, -0.169067413920115 + 0.1790415657178582i, -0.1628948738933707 + -0.1816119404971363i, 0.02978157550996891 + 0.04262334705550256i, -0.04345706462582026 + 0.02828162038081413i, -0.02265557786127289 + -0.04666236764801889i, 0.04714107177944449 + -0.02145286994554814i, 0.02236402470471761 + 0.04681842080535928i, -0.04729498903774116 + 0.02173939977631235i, 0.04751502815852508 + -0.04751844241467342i, 0.04677326569838006 + 0.04637319781115221i, -0.05375017628733026 + 0.03898807467706384i, -0.03911842466114791 + -0.05259262137187445i, 0.05413779425344784 + -0.03875376096967136i, 0.03958006572156292 + 0.05287875261011321i, -0.1959123579154331 + 0.1152702473266644i, -0.1904000292314342 + -0.1559287253901943i, 0.1743551383804932 + -0.1715825512913356i, 0.04518069085267876 + 0.2100360480210003i, 0.7886715491067711 + 0.0186250554204157i, -0.159902037239388 + -0.1839994360167489i, 0.194783821915703 + -0.1427731292914902i, -0.04489136591719057 + 0.02518482175782046i, -0.02362593873992407 + -0.045566508363115i, 0.04817586329477178 + -0.01777124955465136i, 0.01653978701700139 + 0.04853024119573128i, -0.04830118275675168 + 0.01746882704696061i, -0.01680705438952684 + -0.04870979060261692i, 0.04218108997208465 + 0.05143871799764225i, -0.0502412595392465 + 0.04155827180421484i, -0.0331700310279586 + -0.05674981134125896i, 0.05562218216996979 + -0.03341124502531951i, 0.03290143036696498 + 0.05710882815009482i, -0.05594904290774396 + 0.03383813393746703i, -0.134872821669944 + -0.1865348602265873i, 0.06560461647811704 + 0.1491213238028099i, -0.1543526472272893 + 0.04897346920098347i, -0.1690674139201165 + 0.1790415657178472i, -0.1599020372393835 + -0.1839994360167502i, -0.6509212324533222 + -0.2816135754998981i, 0.3016994794535154 + 0.3707113002644887i, -0.02966651580195464 + -0.04285958809110769i, 0.04368873002977677 + -0.02815919361329978i, 0.02250517339246269 + 0.04687706281816777i, -0.04735163176557099 + 0.02129739651325564i, -0.02221221478008752 + -0.04703221183381037i, 0.04750719252507979 + -0.02158395061777181i, -0.04784162760569108 + 0.0474253658015703i, -0.04668340062638959 + -0.04669370082388591i, 0.05405442458137942 + -0.03884650248912184i, 0.03898228472286762 + 0.05289460107391549i, -0.05444196019150098 + 0.03860990125543476i, -0.03944379647785421 + -0.05318347405036096i, 0.1952332425323231 + -0.117853710019941i, -0.1526109204735092 + 0.0525351549814783i, -0.03568594191998528 + -0.1564111202539286i, -0.162894873893359 + -0.181611940497136i, 0.1947838219157022 + -0.1427731292914828i, 0.3016994794535193 + 0.3707113002644895i, 0.6028673353657812 + 0.3801639461038105i, 0.04474964408033698 + -0.0257728811088443i, 0.02421666390976782 + 0.04544457143282204i, -0.04812846749371872 + 0.01837282063754986i, -0.01714143981691673 + -0.04849790414799845i, 0.0482576071657195 + -0.01807085131875214i, 0.01741159024350311 + 0.04867500304211311i, -0.04289740485081017 + -0.05112564375193637i, 0.0499314114347067 + -0.04225905653047613i, 0.03391768429683719 + 0.05655559458114036i, -0.05542164137547257 + 0.03414691314386125i, -0.03365228787972556 + -0.05691879427060867i, 0.05574471595422684 + -0.03457886033189477i, 0.1525968069538521 + 0.1918033025474123i, -0.03621012858184258 + -0.0373514054429251i, 0.04076554172092756 + -0.03181268282776888i, 0.02978157550996868 + 0.04262334705550333i, -0.04489136591718952 + 0.02518482175781987i, -0.02966651580195241 + -0.04285958809110606i, 0.0447496440803349 + -0.02577288110884331i, -0.744631950112659 + -0.1267370815520711i, 0.3773201896903062 + 0.2331266541325207i, -0.1027276726307461 + 0.1796959125817161i, -0.1768112120030458 + -0.1070125666319918i, -0.02008975443580694 + 0.1648766977762729i, -0.1650898958540091 + -0.02259578089641597i, -0.03000847267969491 + -0.04218203731418592i, 0.0412191954371112 + -0.02958871263073519i, 0.02272968977075173 + 0.04582549255435331i, -0.04496288246019536 + 0.02297701909673075i, -0.02250202910836154 + -0.04608992292373413i, 0.04523938999678508 + -0.02329102838245655i, -0.1961650118983897 + 0.1457842897593873i, 0.03841398593767489 + -0.03486153122936388i, 0.03036685699558397 + 0.04167221791250402i, -0.04345706462582219 + 0.02828162038081539i, -0.023625938739925 + -0.04556650836311486i, 0.04368873002977598 + -0.02815919361329991i, 0.0242166639097683 + 0.04544457143282122i, 0.3773201896902928 + 0.2331266541325197i, 0.7360719974054075 + 0.1733675459323022i, -0.1757161969412171 + -0.1082800162613657i, 0.1124557274617433 + -0.1727005695618684i, -0.1636600658313424 + -0.02543479722026158i, 0.02793935567490066 + -0.1637902704096546i, 0.04302469503647204 + -0.02852223916549212i, 0.02813550760964802 + 0.0420513256223061i, -0.04641687650997445 + 0.02114843462635754i, -0.0214232478724909 + -0.045565300841199i, 0.04667293929918584 + -0.02091286189966866i, 0.02172711841282066 + 0.04585118308076899i, -0.1202548022673231 + -0.2128927938006659i, 0.02981566718969052 + 0.04247674792409621i, -0.04515063395283176 + 0.02494794672353451i, -0.0226555778612723 + -0.04666236764801972i, 0.04817586329477125 + -0.01777124955465078i, 0.0225051733924614 + 0.04687706281816717i, -0.04812846749371848 + 0.0183728206375488i, -0.1027276726307522 + 0.1796959125817195i, -0.1757161969412099 + -0.1082800162613697i, -0.842943086785175 + 0.04172302067471873i, 0.04215973556516767 + 0.1540635211758709i, -0.1303058380650114 + 0.1598863184091735i, -0.1584504518803855 + -0.1330794653328389i, 0.02294834941767331 + 0.04626317934512574i, -0.04524877475564974 + 0.02268581269465888i, -0.0152060267984059 + -0.04871110755238831i, 0.04790006929419029 + -0.01558497441518787i, 0.0149402691467909 + 0.04893591519884621i, -0.04822172544703765 + 0.01585096559937093i, 0.2154972257377104 + -0.1147367306145593i, -0.04314129552047759 + 0.02870406590745015i, -0.02377868003346755 + -0.04568915089366452i, 0.04714107177944647 + -0.02145286994554885i, 0.01653978701700145 + 0.04853024119573253i, -0.04735163176557165 + 0.0212973965132549i, -0.01714143981691689 + -0.04849790414799892i, -0.1768112120030535 + -0.1070125666319965i, 0.1124557274617467 + -0.1727005695618709i, 0.04215973556515717 + 0.1540635211758769i, 0.8444451154781266 + 0.0002226071600751614i, -0.1563513897944216 + -0.1340473677415275i, 0.136780196312724 + -0.1548491121335388i, -0.04674989868530611 + 0.0217550424115373i, -0.02151823192704622 + -0.04573081126737644i, 0.04900066383863656 + -0.01396592444846201i, 0.01436436606220182 + 0.04820053723633964i, -0.04921844686769145 + 0.01369505667348232i, -0.01462186623230235 + -0.04852822649499765i, 0.1189288951605393 + 0.213712383149512i, -0.02955257758696452 + -0.04267736000877264i, 0.04532090989257979 + -0.02466660000838256i, 0.02236402470471674 + 0.04681842080535871i, -0.04830118275675066 + 0.01746882704695929i, -0.02221221478008605 + -0.0470322118338094i, 0.04825760716571876 + -0.01807085131875099i, -0.02008975443579509 + 0.1648766977762693i, -0.1636600658313508 + -0.02543479722025991i, -0.1303058380650149 + 0.1598863184091836i, -0.1563513897944193 + -0.1340473677415302i, -0.6658612035353044 + -0.3588440681493852i, 0.2876129526730812 + 0.3372282125261508i, -0.02265943033602364 + -0.04642100236243787i, 0.04540467013096164 + -0.02240328925444514i, 0.01489956602856806 + 0.04882022837048633i, -0.04801139987899267 + 0.01528377529354699i, -0.0146323094854446 + -0.04904339973284587i, 0.04833483133853983 + -0.01554778449019651i, -0.2161374308547739 + 0.1161298649238348i, 0.04324764199660035 + -0.02901160838497489i, 0.02407130670199345 + 0.04582696551264388i, -0.047294989037741 + 0.02173939977631128i, -0.01680705438952393 + -0.04870979060261676i, 0.04750719252507992 + -0.02158395061776957i, 0.01741159024350125 + 0.04867500304211325i, -0.1650898958539998 + -0.0225957808964029i, 0.02793935567490025 + -0.1637902704096557i, -0.1584504518803901 + -0.1330794653328474i, 0.1367801963127285 + -0.1548491121335362i, 0.2876129526730919 + 0.3372282125261562i, 0.6527807061759556 + 0.3801712114266651i, 0.04690071740676486 + -0.02204161134042052i, 0.02179974031791457 + 0.0458774390245941i, -0.0491923637398904 + 0.01422213787412679i, -0.01461952611164343 + -0.04838672330098354i, 0.04941226550648602 + -0.01395074528688192i, 0.01487955965922728 + 0.0487150701841322i, -0.1448292478234282 + -0.1821851739600287i, -0.03935537304529047 + 0.0545086558278562i, -0.04968364134529177 + -0.04469249669008771i, 0.04751502815852619 + -0.0475184424146767i, 0.04218108997208571 + 0.0514387179976397i, -0.04784162760568744 + 0.04742536580157171i, -0.0428974048508124 + -0.05112564375193359i, -0.03000847267969542 + -0.04218203731418973i, 0.04302469503647469 + -0.02852223916549118i, 0.02294834941767308 + 0.04626317934512832i, -0.04674989868530759 + 0.02175504241153602i, -0.02265943033602384 + -0.04642100236244033i, 0.04690071740676631 + -0.02204161134042121i, -0.6701272864339438 + -0.1870225420659977i, 0.3843563368792602 + 0.3011195320066258i, -0.1824964773967196 + 0.1742903661607326i, -0.1741897745490597 + -0.1780706256020184i, 0.04145180472311114 + 0.1607871533979154i, -0.1586016957179824 + 0.04450631308229829i, 0.1779630932908195 + -0.1427140458171124i, -0.05359061798842831 + -0.03834638056056874i, 0.04359767205298964 + -0.04888368946936066i, 0.04677326569838389 + 0.04637319781115375i, -0.05024125953924481 + 0.04155827180421641i, -0.04668340062639131 + -0.04669370082388278i, 0.0499314114347045 + -0.04225905653047844i, 0.04121919543711396 + -0.02958871263073597i, 0.02813550760964702 + 0.04205132562230712i, -0.04524877475565058 + 0.02268581269465832i, -0.02151823192704453 + -0.04573081126737701i, 0.04540467013096237 + -0.02240328925444569i, 0.02179974031791553 + 0.04587743902459408i, 0.3843563368792531 + 0.3011195320066235i, 0.6837040722541522 + 0.1658450212218595i, -0.1715914811907637 + -0.1781452287184037i, 0.1738076769438393 + -0.1714744110208234i, -0.157421507382712 + 0.04129994718467979i, -0.04428468167200321 + -0.1552666931592779i, 0.113024300770295 + 0.200286394208954i, 0.046873711810578 + -0.04707719629199332i, 0.04153950746586119 + 0.05133223835087421i, -0.05375017628732514 + 0.03898807467707156i, -0.03317003102795694 + -0.05674981134125159i, 0.05405442458137368 + -0.03884650248911956i, 0.03391768429683644 + 0.0565555945811351i, 0.02272968977075264 + 0.04582549255435554i, -0.04641687650997563 + 0.02114843462635614i, -0.01520602679840524 + -0.04871110755238914i, 0.04900066383863577 + -0.01396592444846114i, 0.01489956602856886 + 0.0488202283704865i, -0.04919236373988991 + 0.0142221378741285i, -0.1824964773967086 + 0.1742903661607282i, -0.1715914811907568 + -0.1781452287183965i, -0.7792243702971343 + -0.005997934715073036i, 0.03405196567006005 + 0.21720786159751i, -0.2070208765569836 + 0.1414777930049319i, -0.1446922980933196 + -0.202348031663166i, -0.1963544656890722 + 0.1139279967111949i, 0.04702140954937468 + 0.04570961469855871i, -0.05017584657756854 + 0.04160820276080084i, -0.03911842466115796 + -0.05259262137187071i, 0.05562218216996381 + -0.03341124502532055i, 0.0389822847228673 + 0.052894601073911i, -0.05542164137546849 + 0.03414691314386271i, -0.04496288246019942 + 0.0229770190967338i, -0.0214232478724913 + -0.0455653008412019i, 0.04790006929419356 + -0.01558497441518884i, 0.01436436606220225 + 0.04820053723634139i, -0.04801139987899555 + 0.01528377529354882i, -0.014619526111646 + -0.04838672330098591i, -0.1741897745490583 + -0.178070625602017i, 0.1738076769438413 + -0.171474411020817i, 0.03405196567006365 + 0.2172078615975175i, 0.7847686602105013 + -0.01362221963894329i, -0.1420924035275209 + -0.202665222370419i, 0.1980149205140884 + -0.1452099976583627i, -0.1120655810248029 + -0.2015280958573828i, -0.04729391562068966 + 0.04690755733334535i, -0.04132715783868043 + -0.05172947371213033i, 0.05413779425344873 + -0.03875376096967614i, 0.03290143036696783 + 0.05710882815009167i, -0.05444196019149714 + 0.03860990125543774i, -0.03365228787973033 + -0.05691879427060537i, -0.02250202910836354 + -0.04608992292373616i, 0.04667293929918737 + -0.02091286189966866i, 0.01494026914679084 + 0.04893591519884874i, -0.04921844686769183 + 0.01369505667348243i, -0.01463230948544588 + -0.04904339973284705i, 0.04941226550648659 + -0.01395074528688357i, 0.04145180472311153 + 0.1607871533979112i, -0.1574215073827062 + 0.04129994718467827i, -0.2070208765569861 + 0.1414777930049228i, -0.1420924035275116 + -0.2026652223704133i, -0.5794305523609922 + -0.3908983188618896i, 0.2735489758434193 + 0.4038664448678665i, 0.1975206021934352 + -0.1154039852576325i, -0.04752298655269267 + -0.04591863332607179i, 0.05043577037151492 + -0.04208166894894561i, 0.03958006572156556 + 0.05287875261011599i, -0.05594904290774268 + 0.03383813393746758i, -0.03944379647785543 + -0.0531834740503587i, 0.05574471595422555 + -0.03457886033189627i, 0.04523938999678728 + -0.02329102838245709i, 0.02172711841281893 + 0.04585118308077092i, -0.04822172544704037 + 0.01585096559936905i, -0.01462186623230064 + -0.04852822649499815i, 0.04833483133854115 + -0.015547784490196i, 0.01487955965922729 + 0.04871507018413282i, -0.1586016957179807 + 0.04450631308230069i, -0.04428468167200348 + -0.1552666931592789i, -0.144692298093329 + -0.2023480316631716i, 0.1980149205140862 + -0.145209997658369i, 0.2735489758434222 + 0.4038664448678754i, 0.6018651612875096 + 0.3614132990747745i;
];
Z = [
4.363283813063463e-013 + -532.1632843965498i, 2.329962659190866e-013 + -10.63786138830847i, 1.619818859784582e-013 + -124.8855568353299i, -5.592224168931314e-014 + 0.5510444922264284i, -3.306862552343881e-013 + 200.4786346930177i, -1.695255585749083e-013 + 4.571023525277397i, -4.775625478430648e-013 + -80.23298977660252i, -1.688616011816683e-013 + 13.25706275765797i, -1.4656243277198e-013 + 133.7314790489662i, -4.221705074631176e-014 + -3.574561379171352i, -2.739876211987219e-013 + -187.4830588455853i, 2.727815688802671e-013 + -4.075090715769153i, -4.329774481753175e-013 + 84.45234280184289i, 2.320793309709073e-014 + -14.61932423245595i, 6.41217731664116e-014 + -133.6656546900736i, -1.212295269874115e-014 + 4.666261729091753i, 9.042139654996431e-013 + 184.9353967592521i, -1.548794090509875e-013 + 3.334435595193464i, 4.205870897315203e-013 + -80.90347728632065i, 7.445885033941389e-014 + -10.63786138830946i, -2.50658059129949e-014 + 17.53927798927075i, 6.980695264903195e-014 + 59.52721452040361i, -1.463723085129272e-014 + -0.07568935122882825i, -4.982457638676713e-014 + -27.53697389044975i, 3.783224166288523e-014 + 1.075753596103027i, 5.966427620287075e-015 + -18.8821883766484i, 6.281777123395728e-015 + 0.3602832133560757i, 8.656007555231023e-014 + 3.634380245412034i, -9.90325023165268e-015 + -0.09714478229516109i, 1.329588832406584e-013 + -5.095170848572716i, 3.287878339196788e-015 + -0.1107475058361595i, 4.215167038776413e-014 + 2.295135735539021i, -4.21060183564476e-015 + -0.2110101453261967i, 7.134413291166229e-014 + -1.929282692742676i, -5.586547803613917e-015 + 0.067351158878632i, -1.508781201547071e-014 + 2.669291611261947i, -6.904235498130502e-015 + 0.04812808397917829i, 5.169439312392535e-014 + -1.167732660561928i, -2.833347014106744e-014 + -124.8855568353296i, -2.08966316295424e-013 + 59.52721452040377i, 5.797012748666208e-013 + -76.3792928142732i, -8.141643842365868e-014 + -0.8885720944791316i, 1.617175247194817e-013 + -323.2764737512222i, 1.179530193713543e-013 + 12.62905032909823i, 1.888080587158573e-013 + -221.6716803885504i, 4.164585538485512e-016 + 4.229625492943873i, 1.63219895186261e-013 + 42.66662105584969i, 2.144845403795997e-014 + -1.140452933886895i, 2.7650339130324e-013 + -59.81589958432728i, -6.921896650291711e-014 + -1.300145154144829i, -1.865605360206879e-013 + 26.94426011796815i, 9.346744393899875e-014 + -2.477200871023991i, 2.449833082188211e-013 + -22.64924636455791i, 6.275872496422066e-015 + 0.7906840174951608i, 4.245794408354456e-014 + 31.33674683848639i, -5.774911969936526e-014 + 0.5650104234075857i, 1.110396016725226e-013 + -13.70885916123005i, 9.338383709740105e-014 + 0.5510444922283169i, 3.673222844380641e-014 + -0.07568935122788859i, -6.995256242733119e-014 + -0.8885720944789993i, 5.070317182072558e-014 + 6.742502184050672i, -2.043989646263059e-014 + 64.18848600425947i, -2.976918770452049e-014 + 0.03177522733389196i, -1.650287275005222e-013 + -0.5577353683997793i, -1.490736958048754e-014 + -0.01866278190662393i, -6.151438746952188e-014 + -0.188262021019717i, -9.761799675198084e-015 + 0.005032129774081249i, 8.605413552452522e-016 + 0.263931426117825i, 2.068782637234824e-014 + 0.005736755070481308i, -1.342157719418348e-015 + -0.1188887412759427i, -5.300156500303967e-015 + 0.01093038641621063i, 2.297147263984269e-015 + 0.09993740019085493i, -9.938452925351371e-015 + -0.003488809463459823i, 4.503106348349671e-014 + -0.138270075661899i, 1.605414535571188e-015 + -0.002493048640337656i, 2.976673055980665e-014 + 0.06048888883578439i, -3.475061415545109e-013 + 200.4786346930251i, -2.296544201777761e-013 + -27.53697389043803i, 2.939517090174611e-013 + -323.2764737512287i, -3.060065399944962e-013 + 64.18848600426367i, 2.938434597993327e-013 + -513.9143477723495i, 4.209908571682946e-014 + 11.56032640342779i, 2.327957460831549e-013 + -202.9128804550632i, 4.805122407836098e-014 + -6.789812962438887i, -1.210503233396795e-014 + -68.49267794318014i, 6.999861414932179e-014 + 1.83076778904653i, 5.084644564515901e-013 + 96.02239511858468i, -2.861549987816926e-013 + 2.087121530692691i, 1.39225418962218e-013 + -43.25358998686869i, 7.972458948080593e-014 + 3.976646126210532i, -5.341665450665851e-014 + 36.35879466648639i, 3.356112530730126e-014 + -1.269283638254406i, -5.340932157621414e-013 + -50.30482361869346i, -3.800508794692086e-014 + -0.9070102227462826i, -2.343986021083332e-013 + 22.00680710708511i, 6.369983334018441e-014 + 4.571023525276945i, -4.209059466573008e-014 + 1.075753596103017i, 4.929532386210954e-014 + 12.62905032909919i, -2.761557844909962e-015 + 0.03177522733370328i, 2.384969324522019e-013 + 11.56032640342697i, -4.177186715922801e-014 + -1.945682856102992i, 1.464519352979163e-014 + 64.17172314213967i, -9.679117134641668e-015 + -0.1548114793823767i, 3.37139867994119e-014 + -1.561670823320588i, -1.578505068972955e-014 + 0.04174251506363809i, 8.128378262826197e-014 + 2.189363554823087i, -1.557297489524042e-015 + 0.0475875216700238i, 4.308755360671821e-014 + -0.9862057003087889i, -1.931159667157967e-014 + 0.09066975092781382i, -4.972225776596161e-014 + 0.8290008091788153i, 2.600275034753916e-015 + -0.0289403753088865i, -6.558530167579029e-014 + -1.146978052419418i, 2.458826312073248e-014 + -0.02068033923907484i, -2.172585780909948e-014 + 0.501767483023603i, 3.438180730822889e-013 + -80.23298977659938i, -1.362161804204325e-013 + -18.88218837664457i, 1.861811622044543e-013 + -221.6716803885493i, -2.423595327836546e-013 + -0.557735368397604i, 1.970604000429884e-013 + -202.9128804550423i, -1.218293830217718e-013 + 64.17172314213985i, 1.544080969982391e-013 + 18.330430152564i, -2.300550203953941e-014 + 2.717331856618082i, -6.56205359711264e-014 + 27.41126106855772i, 3.93307966158969e-015 + -0.7326864029093704i, 2.67922125685703e-013 + -38.4288513808799i, -4.872652267482558e-014 + -0.8352810084067754i, -2.351135111367897e-013 + 17.31039698941438i, 8.065452697813956e-014 + -1.591482773827792i, -3.283760352238673e-014 + -14.55105472107982i, 3.944962701600301e-014 + 0.5079765666157236i, -2.928446813366209e-013 + 20.13235719415079i, -1.983041940362934e-014 + 0.3629921039725813i, -7.080401513451747e-014 + -8.807284651467141i, -8.021854045604839e-014 + 13.25706275765922i, 9.97857388738412e-015 + 0.3602832133560507i, -4.71297548865705e-014 + 4.229625492944351i, -4.772283276659665e-016 + -0.01866278190660899i, 5.755640992037373e-015 + -6.789812962439598i, 1.476261344247513e-014 + -0.1548114793823961i, -3.958458933342667e-014 + 2.71733185661835i, 6.257066396650347e-015 + 14.26635048863314i, -5.728243939823773e-014 + 59.87418521301328i, 8.144399173569609e-015 + -0.545533431482837i, -5.125806176241055e-015 + -28.61281863360653i, 4.038758058729553e-014 + 0.9029206607595904i, -1.205344290520208e-013 + -18.71216384711491i, -9.909585822059717e-015 + 0.4895546641725251i, -7.967551694032362e-014 + 4.476037582361582i, 7.618072802464381e-015 + -0.1562582616649219i, -9.717896262417237e-015 + -6.192898147515159i, 6.265124001759656e-015 + -0.111659641572908i, -3.805722383912255e-014 + 2.709200108352525i, -6.367538102770695e-013 + 133.7314790489693i, 5.938066928994723e-014 + 3.634380245412054i, 7.854921699004633e-014 + 42.66662105585012i, 2.327538476840136e-014 + -0.1882620210194196i, 3.552108395130183e-013 + -68.49267794318499i, 1.32073668234806e-014 + -1.561670823320485i, -1.683804906125398e-014 + 27.41126106855633i, 1.95967656989238e-013 + 59.87418521301371i, 1.92232633280694e-013 + -37.62704803999463i, 1.580744314555088e-013 + -5.503103816922493i, -2.661960600488077e-013 + -288.6336608327091i, 1.4416665022693e-013 + 9.108270635399093i, -1.331009059581568e-013 + -188.7601645421291i, -5.707450923657895e-014 + 4.938414376690961i, -2.211321169291584e-013 + 45.15231896463442i, 1.007942504108855e-014 + -1.576265333284181i, 0 + -62.47126109351394i, 8.549144264804284e-014 + -1.126373865055138i, -1.141711272004443e-013 + 27.32923153134184i, -1.671287076782427e-013 + -3.57456137917237i, -1.107478585522094e-014 + -0.09714478229515187i, -8.431430424533813e-014 + -1.140452933887055i, 8.592530891810953e-015 + 0.005032129774059469i, 6.607906099902239e-014 + 1.830767789047065i, 7.579699368014241e-015 + 0.04174251506363186i, -5.683294508645784e-014 + -0.7326864029095664i, -3.047173944269208e-014 + -0.5455334314828556i, 2.14146717534576e-013 + -5.503103816924651i, 7.641156159500602e-015 + 3.839770321487506i, 3.051175976895434e-013 + 55.62372730436658i, -1.05654341538108e-013 + 0.16736717832491i, 9.987438449001482e-014 + -3.46852409027046i, -5.230721579535492e-015 + -0.1320008191374973i, -3.962882780864505e-014 + -1.206894082728287i, -2.133455728359078e-015 + 0.04213261571561208i, 1.130427535559943e-013 + 1.669818894918173i, -8.126228553104431e-015 + 0.03010728980737039i, -4.507455878793864e-014 + -0.7304937724975082i, -4.915738896267436e-013 + -187.4830588455837i, 1.46188646454684e-013 + -5.09517084857212i, -4.751034861317244e-014 + -59.81589958432844i, 6.103397165200503e-014 + 0.2639314261174354i, 2.484241525043167e-013 + 96.02239511858456i, -1.051767450183829e-013 + 2.189363554823203i, -2.66885154195591e-013 + -38.42885138088032i, 2.330185763622222e-013 + -28.61281863360555i, -2.889588883412997e-013 + -288.6336608327144i, 4.082500553192039e-013 + 55.62372730437161i, -1.502743149310175e-012 + -404.6276692723608i, 5.191948228700036e-014 + 8.778282763742276i, -8.669859346366059e-013 + -181.9214827062299i, 2.168652764152353e-014 + -6.92334379434964i, 1.310527117828357e-013 + -63.30068792521325i, -6.395814961138863e-014 + 2.209824040944732i, 1.076609600645954e-012 + 87.58074653900556i, -1.062109336523486e-014 + 1.579104731429469i, 3.718410506368382e-013 + -38.31384956305168i, -9.812426152309245e-014 + -4.075090715770208i, -6.542543724258667e-015 + -0.1107475058361727i, -2.05208591090966e-014 + -1.300145154145171i, -2.132587672432167e-014 + 0.005736755070477695i, 8.133198383274527e-014 + 2.087121530693205i, -6.595864643202453e-015 + 0.04758752167004381i, -1.205115987490244e-014 + -0.835281008407027i, -8.642512974923514e-014 + 0.9029206607595244i, -2.629581233777621e-013 + 9.108270635400359i, -1.166553407376856e-014 + 0.1673671783250519i, -2.364346847420607e-013 + 8.77828276374451i, 1.071191787308954e-013 + -5.095120141147505i, 5.573615550654175e-014 + 61.54805376601927i, -1.019161588339304e-014 + -0.1504843256342553i, -2.004954083815726e-014 + -1.375890265970141i, -7.561762920807625e-017 + 0.04803226451867922i, 8.942400413486587e-014 + 1.903636446414146i, -6.457123910805613e-015 + 0.03432308403035877i, -9.42136090458648e-015 + -0.8327816717315532i, -2.948051091117952e-013 + 84.45234280184448i, 3.471328744670189e-014 + 2.295135735539075i, 8.814390172926094e-014 + 26.94426011796767i, 7.577488068083802e-015 + -0.118888741275487i, 1.136921835964862e-013 + -43.25358998686867i, -7.463200244779136e-015 + -0.9862057003087874i, 1.943658109682839e-014 + 17.31039698941397i, -1.967499280249327e-014 + -18.71216384711233i, -2.32306953338144e-013 + -188.7601645421323i, 1.762620347167725e-013 + -3.468524090269034i, -4.946406522487501e-013 + -181.9214827062273i, 1.723474175694595e-013 + 61.54805376601906i, -1.405259489385782e-013 + 28.75765904176782i, -3.8335829458838e-014 + 3.118643176601613i, -6.521794257775781e-014 + 28.51400484160072i, -3.73057574221128e-016 + -0.9954225688724192i, 2.75530941930939e-014 + -39.45103776966259i, 5.480192553281879e-014 + -0.7113129647226804i, 5.051207174145056e-015 + 17.25860063629894i, 9.734027406783037e-015 + -14.61932423245701i, -4.336280651427716e-015 + -0.211010145326096i, 7.155309044135375e-014 + -2.477200871024246i, -4.759091334718549e-016 + 0.01093038641633122i, -1.507069370374158e-014 + 3.976646126211094i, -1.747038584516907e-014 + 0.09066975092784795i, 2.399359140665921e-014 + -1.591482773828008i, 1.902025160535436e-014 + 0.4895546641725542i, 6.478585400877472e-014 + 4.938414376691285i, -7.415384918831913e-015 + -0.132000819137452i, 1.318275229316359e-013 + -6.923343794350348i, -5.809158815610287e-015 + -0.150484325634236i, 4.671522729012558e-017 + 3.118643176601862i, 4.120824250124282e-014 + 12.46450987678707i, -2.850831362687906e-013 + 56.80898004533997i, 6.084416901945145e-015 + -0.7010668700531286i, -5.990320578826665e-013 + -27.78499949107484i, 8.765183218828216e-014 + 0.8725259883403096i, -2.121261653546414e-013 + -21.17011544055669i, -5.679088350065818e-013 + -133.6656546900869i, 3.935735896867761e-014 + -1.929282692742147i, 6.876089079635272e-014 + -22.64924636456232i, -4.559505076158809e-014 + 0.09993740019196835i, 2.460009846739691e-013 + 36.35879466649505i, -3.05328179859331e-014 + 0.8290008091792467i, -1.677405761128324e-013 + -14.55105472108231i, 5.455614841863887e-014 + 4.476037582361664i, 5.245906203847741e-014 + 45.152318964638i, -1.654636401377634e-014 + -1.206894082728071i, -3.577127785704872e-013 + -63.30068792521908i, 4.843258482152134e-014 + -1.375890265969933i, -2.694827826164276e-013 + 28.51400484160357i, -1.075210577796265e-013 + 56.80898004534027i, -9.119282390154928e-013 + -46.7791910257972i, 1.250614852155241e-013 + -6.409910655040485i, -4.48887656048578e-013 + -254.0404801536838i, 2.00416682994973e-013 + 7.977575133510864i, -2.057322974987377e-013 + -193.5600644209762i, 1.402704505900933e-013 + 4.666261729092639i, -2.683435756122944e-015 + 0.06735115887863814i, 2.32610403447443e-014 + 0.7906840174952405i, -1.761757435996971e-015 + -0.003488809463405718i, -4.23711754782818e-014 + -1.269283638254576i, 1.115452362248871e-015 + -0.02894037530890375i, 3.378630612217303e-014 + 0.507976566615893i, -2.152872446627717e-015 + -0.1562582616649366i, -2.842889908323028e-014 + -1.576265333284377i, 7.304271657874172e-015 + 0.0421326157156026i, -2.284167000806289e-015 + 2.209824040945134i, -4.097137783233552e-015 + 0.04803226451866956i, -3.714857126751316e-014 + -0.9954225688724434i, 2.660078861056343e-014 + -0.7010668700530172i, 8.180626868135265e-014 + -6.409910655039305i, -7.512440502849517e-015 + 2.074970180739836i, 1.235175797249998e-013 + 51.11971301488756i, -3.006332476889779e-015 + 0.1608616904458509i, 1.429001246391787e-013 + -3.902990400523795i, -1.060958733570967e-013 + 184.9353967592407i, -4.827416395607585e-014 + 2.669291611261879i, 1.968788804508876e-013 + 31.33674683848403i, -1.217892583359862e-014 + -0.1382700756603088i, 3.582885357598304e-013 + -50.30482361868131i, 4.511446259871478e-014 + -1.146978052419344i, 1.691749366311052e-013 + 20.13235719415196i, 9.481435394501944e-014 + -6.192898147514034i, -3.891619158860145e-013 + -62.47126109350791i, 1.777838640030405e-014 + 1.669818894917584i, -2.195501620020092e-013 + 87.58074653900087i, -7.608201441293187e-014 + 1.903636446413492i, -2.519891534335908e-013 + -39.45103776965986i, -4.108212241932234e-013 + -27.78499949107327i, -5.283937872223748e-013 + -254.0404801536892i, -4.89111292376415e-014 + 51.11971301488838i, -7.435632412902892e-013 + -411.3456406055272i, -5.257098270508714e-014 + 6.375343320434718i, -3.032839188576244e-013 + -154.6850820151223i, 6.10680146284556e-014 + 3.334435595196589i, -8.578603807780707e-015 + 0.04812808397925465i, -1.877546482245356e-015 + 0.5650104234083347i, 8.834650858716582e-015 + -0.002493048640372282i, 1.993898795382541e-014 + -0.9070102227473433i, -2.326828077514149e-015 + -0.02068033923914454i, 1.289459696885226e-014 + 0.3629921039728878i, 8.979485847134724e-015 + -0.1116596415729992i, 5.277938419318458e-014 + -1.126373865056137i, -7.784266097004589e-015 + 0.03010728980739785i, 9.85984612505475e-014 + 1.57910473143083i, -8.092702097923541e-015 + 0.03432308403037814i, 3.944571526363821e-014 + -0.7113129647233383i, 6.759116617406269e-014 + 0.8725259883402143i, 2.02996106747738e-013 + 7.977575133509109i, -1.386331857204506e-014 + 0.1608616904457078i, -5.988610411092888e-014 + 6.375343320430387i, 1.159760336409635e-014 + -6.247218076954501i, 3.179481263308867e-014 + 60.01155730204726i, 1.600770700679122e-013 + -80.903477286328i, 1.322493104428388e-014 + -1.167732660561587i, -7.326722953793497e-014 + -13.70885916123286i, -8.799968371138964e-015 + 0.06048888883650651i, -1.516073180325614e-013 + 22.00680710709038i, 1.614142735902762e-014 + 0.5017674830238411i, -1.448880301696603e-013 + -8.807284651468615i, 1.358242964741159e-014 + 2.709200108352422i, -2.059059975864087e-013 + 27.32923153134426i, -9.757633198314148e-015 + -0.7304937724975587i, -1.297328470456194e-013 + -38.31384956305455i, 4.521344687320923e-014 + -0.832781671731454i, -3.032779666601796e-013 + 17.25860063630076i, -1.465438649918282e-013 + -21.17011544055544i, 1.573237418185356e-013 + -193.5600644209757i, 9.671111401851995e-015 + -3.902990400524719i, 2.119925480402818e-013 + -154.6850820151186i, -3.98143653540778e-014 + 60.01155730204955i, -5.034788239564561e-014 + 30.2916897217036i;
];
