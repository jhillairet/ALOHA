% #####################################################################
% #                                                                   # 
% #                         ALOHA                                     #
% #                                                                   #
% #              (Advanced LOwer Hybrid Antenna)                      #
% #                                                                   #
% #####################################################################
%  
%  
%  Code de calcul du couplage des antennes à la fréquence hybride basse 
%  avec un plasma de tokamak ALOHA (Advanced LOwer Hybrid Antenna).
%  
%  ALOHA est un code de calcul base sur la theorie lineaire du couplage
%  des ondes electromagnétiques à la fréquence hybride basse. Il permet 
%  de modéliser le couplage d'une antenne de type 'grill', c'est-a-dire
%  un réseau de guides d'ondes déphasés, avec le plasma de bord d'un 
%  tokamak.
%  
%  Le coeur du code a ete developpé entre 1994 et 1996 par S.Berio et
%  Ph.Bibet. Le projet a ensuite ete repris par D.Voyer et Ph.Bibet. 
%  Quelques modifications ont ete apportees lors du stage d'O.Izacard.
%  Depuis 2008, le code fortran 1D a ete porte en fortran 90 et n'utilise
%  plus les librairies proprietaires NAG.
%  
%  Le plasma de bord est modelise dans ALOHA comme :
%   - un plasma 'froid' ; 
%   - un plasma inhomogène en densite suivant l'axe radial 
%     (pour chacune des lignes toroidales consituant l'antenne)
%  
%  La propagation des ondes electromagnetiques dans les guides de l'antenne
%  est modelise sous la forme d'une somme de modes TE et TM.
%  Le propagation des ondes electromagnetiques dans le plasma est 
%  est modelise sous la forme d'une spectre d'ondes planes.
%  
%  En regime harmonique, les hypotheses realisees sur le plasma permettent 
%  de deduire des equations de Maxwell Faraday et de Maxwell Ampere, 
%  un systeme de 4 equations differentielles couplées, reliant les 
%  composantes tangentielles dans le domaine spectrale des champs E et H. 
%  
%  Au niveau de l'interface antenne/plasma, les conditions aux limites 
%  sont exprimees dans le domaine spectrale grace a un terme matriciel 
%  d'admittance Ys, reliant les composantes tangentielles du champ E 
%  aux composantes tangentielles du champ H.
%  
%  Dans l'hypothese 1D, on suppose que la propagation des ondes EM
%  dans le plasma est reduite au plan xOz (direction radiale et toroidale).
%  Dans ce cas, l'admittance de surface Ys est reduite a un scalaire et 
%  on obtient alors une equation differentielle decouplee du second ordre en E.
%  Sa resolution analytique depend ensuite de l'hypothese que l'on fait 
%  sur le profil de densite electronique devant l'antenne. [Brambilla1978]
%  
%  Dans l'hypothese 2D, les termes de l'admittance de surface sont obtenus
%  selon l'hypothese utilisee pour modeliser le profil de densite 
%  electronique devant l'antenne. [Bers1983, Voyer200x]
%   
%  
%  NB : au prealable, le repertoire contenant les librairies
%  d'ALOHA doit etre charge dans le PATH matlab. Pour cela, 
%  on doit effectuer dans matlab : 
%  addpath(genpath([chemin_du_code_ALOHA_V2/libaloha'])); 
%  Cette commande peut etre rajoutee dans le fichier startup.m de 
%  matlab pour effectuer cette operation a chaque demarrage de matlab
%  (cf help matlab)
%  
%  AUTHORS: Ph.Bibet, S.Berio, O.Izacard, D.Voyer, J.Hillairet
%  
%  REFERENCES:
%  Coupling theory 2D/3D:
%   - M.Brambilla, Lower Hybrid waves propagation in cold plasma, 1978
%   - A.Bers, K.S.Theilaber, Three dimensional theory of waveguide plasma coupling, 1983
%  SWAN:
%   - D.Moreau, T.K.Bguyen, Couplage de l'onde lente au voisinage de la fréquence 
%     hybride basse dans les grands tokamaks, 1984
%  ALOHA:
%   - S.Berio, Developpement de coupleurs a la frequence hybride pour la generation
%     non inductive de courant dans un tokamak, 1996
%   - D.Voyer, Advanced modeling of the coupling of LH antenna to tokamak plasmas, 200x  
%     