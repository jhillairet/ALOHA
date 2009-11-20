C ********************************************************************* 
C Programme de calcul de matrice de scattering
C Couplage Grill Plasma avec differents cas poloidaux
C ********************************************************************* 

C      Subroutine SMinhomop(nbgp,Pgg,ag,bg,zg,k0g,X0g,D0g,Nmhg,Nmeg,
C     &Nmhegn,prheg,prhg,preg,RZhegg,Sg)

      Common/err/epsabs
      Common/nlim/nlimit
      Common/dim/b,z
      Common/idem/bcte
      Common/onde/k0
      Common/plasma/X0,D0
C      Common/dimg/bg,zg
C      Common/ondeg/k0g
C      Common/plasmag/X0g,D0g

      Integer *4 nbgp,i,j,q,intg,ni,nlimit
      Integer *4 Pg,Nmh,Nme,Nmhgn,Nmegn,Nmheg,Nmhegn,prhe,prhegn
      Integer *4 Pgg(5),Nmhg(5),Nmeg(5),prheg(5),Nmhg1(5)
      Integer *4 prh(60),pre(60),prhg(5,60),preg(5,60)

      Real *8 pi,k0,f,cl,me,Eps0,qe,epsrel,epsabs
      Real *8 a,ac,b(60),bc,b0,ei,e(60),z(60),pc
      Real *8 n0,nc,Dn,D0,X0
      Real *8 k0g(5),D0g(5),X0g(5)
      Real *8 ag(5),bg(5,60),eg(5,60),zg(5,60)
      Real *8 SMAmpli(2*10/2*60*5,2*10/2*60*5)
      Real *8 SMPhas(2*10/2*60*5,2*10/2*60*5)
      Real *8 SMAmplig(2*10*60*5,2*10*60*5)

      Complex *16 S(2*10*60,2*10*60),Spr(2*10/2*60,2*10/2*60)
      Complex *16 Sg(2*10*60*5,2*10*60*5),Sprg(10*60*5,10*60*5)
      Complex *16 RZheg(2*10*60,2*10*60),RZhegpr(2*10/2*60,2*10/2*60)
      Complex *16 RZhegg(2*10*60*5,2*10*60*5)
      Complex *16 RZhegprg(10*60*5,10*60*5)

      Character *1 repspc,repgril,reppr,sortie,bcte,rep

C Declaration des constantes

      pi=4.*ATAN(1.)
      cl=3.E8
      me=9.1091E-31
      Eps0=8.854E-12
      qe=1.6021E-19
      pc=1.*1.E-3

      Print *,'voulez vous une sortie sur ecran ("o/n")'
      Read *,sortie

      Write (*,*)  'Les grills sont repere du bas vers le haut'
      Print *,'Nombre de grill juxtaposes poloidalement (max. 5)'
      Read *,nbgp
      


C oblige de considere les deux cas sinon
C erreur en chargement du fichier de sauvegarde ds les prg l'utilant
C "input statement requires too much data"
      Print *,'voulez vous uniquement les modes propageants ("o/n")'
      Read *,reppr

C Domaine d etude
      Print *,'frequence d etude (GHz) : '
      Read *,f
C      f=3.7
      f=f*1.E9
      k0=2.*pi*f/cl
      Do 100 q=1,nbgp
         k0g(q)=k0
 100  Continue
C      Write (*,*) 'k0g',k0g

C Caracteristique des Plasmas

      If (nbgp.EQ.1) Then
         repspc='n'
         Goto 1
      EndIf
      Print *,'plusieurs plasma differents ("o/n")'
      Read *,repspc
C      repspc='n'

   1  If (repspc.EQ.'n'.OR.repspc.EQ.'N') Then
   5     Print *,'densite n0 (cm-3) a la bouche des grills '
         Read *,n0
C         n0=1.E12
         n0=n0*1.E6
         nc=(2.*pi*f)**2*me*Eps0/qe**2
         X0=n0/nc
         If (X0.LT.1) Then
            Print *,'n0 inferieure a la densite de coupure nc= ',nc
C            Goto 5
         EndIf
  10     Print *,'Gradient de densite Dn (cm-4) : '
         Read *,Dn
C         Dn=1.E12
         Dn=Dn*1.E8
         D0=k0*n0/Dn
         Do 110 q=1,nbgp
            X0g(q)=X0
            D0g(q)=D0
 110     Continue
      Else   
         Do 120 q=1,nbgp
  15        Print *,'densite  n0 (cm-3) a la bouche du grill  ',q
            Read *,n0
C            n0=1.E12
            n0=n0*1.E6
            nc=(2.*pi*f)**2*me*Eps0/qe**2
            X0=n0/nc
            If (X0.LT.1) Then
            Print *,'n0 inferieure a la densite de coupure nc= ',nc
C               Goto 15
            EndIf
  20        Print *,'Gradient de densite Dn (cm-4) : '
            Read *,Dn
C            Dn=1.E12
            Dn=Dn*1.E8
            D0=k0*n0/Dn
            X0g(q)=X0
            D0g(q)=D0
 120     Continue
      EndIf


C Caracteristique des grills

      Print *,'les guides des grills sont tous identiques ("o/n")'
      Read *,repgril
C      repgril='o'

      If (repgril.EQ.'o'.OR.repgril.EQ.'O') Then
C precision
  25     Print *,'hauteur (grand cote des guides) a=(cte) (mm) : '
         Read *,a
C         a=72.
         a=a*1.E-3
         Print *,'nombre de guide composant les grills : '
         Read *,Pg
C         Pg=6
         Print *,'largeurs des guides toutes identiques ("o/n")'
         Read *,bcte
C         bcte='o'
         If (bcte.EQ.'o'.OR.bcte.EQ.'O') Then
  30        Print *,'largeur (petit cote des guides) b=(cte) (mm) : '
            Read *,b0
C            b0=8.
            b0=b0*1.E-3
            bc=cl/f*0.5
            bc=DNINT(bc/pc)*pc
            Do 320 j=1,Pg  
               b(j)=b0
 320        Continue
         Else
            Do 325 j=1,Pg 
               Print *,'largeur (petit cote) b(mm) du guide ',j
               Read *,b0
               b0=b0*1.E-3
               b(j)=b0
 325        Continue
         EndIf
         Print *,'parois separant les guides toutes identiques ("o/n")'
         Read *,rep
C         rep='o'
         If (rep.EQ.'o'.OR.rep.EQ.'O') Then
  35        Print *,'largeur d une paroie e=(cte) (mm) : '
            Read *,e0
C            e0=2.
            e0=e0*1.E-3
            Do 330 j=1,Pg-1
               e(j)=e0
 330        Continue
         Else
            Do 335 j=1,Pg-1
               Print *,'largeur de paroi e(mm) a droite du guide ',j
               Read *,e0
               e0=e0*1.E-3
               e(j)=e0
 335        Continue
         EndIf

C Calcul des coordonnees en z
         z(1)=0
         Do 340 j=2,Pg
            z(j)=z(j-1)+b(j-1)+e(j-1)
            z(j)=DNINT(z(j)/pc)*pc
 340     Continue
         Do 300 q=1,nbgp
            Pgg(q)=Pg
            ag(q)=a
            Do 310 i=1,Pg
               bg(q,i)=b(i)
               zg(q,i)=z(i)
 310        Continue
 300     Continue
      Else   
         Do 400 q=1,nbgp
  40        Print *,'grand cote des guides a=(cte) (mm) du grill ',q
            Read *,a
C            a=72.
            a=a*1.E-3
            Print *,'nombre de guides composant le grill ',q
            Read *,Pg
C            Pg=6
            Print *,'largeurs des guides toutes identiques ("o/n")'
            Read *,bcte
C            bcte='o'
            If (bcte.EQ.'o'.OR.bcte.EQ.'O') Then
  45           Print *,'petit cote des guides b=(cte) (mm) du grill ',q
               Read *,b0
C               b0=8.
               b0=b0*1.E-3
               bc=cl/f*0.5
               bc=DNINT(bc/pc)*pc
               Do 420 j=1,Pg  
                  b(j)=b0
 420           Continue
            Else
               Do 425 j=1,Pg 
                  Print *,'petit cote b(mm) du guide ',j,'du grill',q
                  Read *,b0
                  b0=b0*1.E-3
                  b(j)=b0
 425           Continue
            EndIf
            Print *,'parois entre les guides toutes identiques ("o/n")'
            Read *,rep
C            rep='o'
            If (rep.EQ.'o'.OR.rep.EQ.'O') Then
 505           Print *,'epaisseur de paroie e=(cte) (mm) du grill ',q
               Read *,e0
C               e0=2.
               e0=e0*1.E-3
               Do 430 j=1,Pg-1
                  e(j)=e0
 430           Continue
            Else
               Do 435 j=1,Pg-1
               Print *,'paroi e(mm) a droite du guide ',j,'du grill',q
                  Read *,e0
                  e0=e0*1.E-3
                  e(j)=e0
 435           Continue
            EndIf
C Calcul des coordonnees en z
            z(1)=0
            Do 440 j=2,Pg
               z(j)=z(j-1)+b(j-1)+e(j-1)
               z(j)=DNINT(z(j)/pc)*pc
 440        Continue
            Pgg(q)=Pg
            ag(q)=a
            Do 410 i=1,Pg
               bg(q,i)=b(i)
               zg(q,i)=z(i)
 410        Continue
 400     Continue
      EndIf

      Nmhgn=0
      Nmegn=0
      Nmhg(0)=0
      Nmeg(0)=0
      Pgg(0)=0
      prhegn=0
      prheg(0)=0
      prhe=0

C Nombre de mode
      If (repgril.EQ.'o'.OR.repgril.EQ.'O') Then
         Print *,'nombre de mode TEmn des guides des grills '
         Read *,Nmh
         Print *,'nombre de mode TMmn des guides des grills '
         Read *,Nme
C         Nme=0
         Do 500 q=1,nbgp
            Nmhg(q)=Nmh
C            Nmhg1(q)=Nmh
            Nmeg(q)=Nme
C         write (*,*) 'q,Nmh',q,Nmh
 500        Continue
      Else
         Do 510 q=1,nbgp
            Print *,'nombre de mode TEmn des guides du grill ',q
            Read *,Nmh
            Print *,'nombre de mode TMmn des guides du grill ',q
            Read *,Nme
C            Nme=0
            Nmhg(q)=Nmh
            Nmeg(q)=Nme
 510        Continue
      EndIf
C      write (*,*) 'Nmhg(q)=',(Nmhg(q),q=1,nbgp)

C Couplage grill/plasma
      Write (*,*) 'Methode d integration (V2): '
      Print *,'1:AKF, 2:ARF 3:BDF 4:AHF 5:BAF : '
      Read *,intg
      If (intg.EQ.5) Then
         Print *,'precision : ni=  '
         Read *,ni
      ElseIf (intg.EQ.4) Then
         Print *,'erreur relative : epsrel=  '
         Read *,epsrel
         Print *,'nb fct eval : nlimit=  '
         Read *,nlimit
      Else
         Print *,'erreur relative : epsrel=  '
         Read *,epsrel
         Print *,'erreur absolue : epsabs=  '
         Read *,epsabs
      EndIf


C      write (*,*) 'Nmhg(q)=',(Nmhg(q),q=1,nbgp)

      If (repspc.EQ.'n'.OR.repspc.EQ.'N') Then

C      write (*,*) 'boucle 0 Nmhg(q)=',(Nmhg(q),q=1,nbgp)

         If (repgril.EQ.'o'.OR.repgril.EQ.'O') Then
      Call SMGrill(intg,ni,epsrel,Nmh,Nme,Pg,a,prh,pre,S,RZheg,
     &Spr,RZhegpr)

C         Do 2 q=1,nbgp
C le programme decale le Nmhg de Nme
C semble etre arrange
C            Nmhg(q)=Nmhg1(q)
C            Nmeg(q)=Nme
C   2        Continue
C      write (*,*) 'boucle 1 Nmhg(q)=',(Nmhg(q),q=1,nbgp)
C      write (*,*) 'boucle 1 Nmhg1(q)=',(Nmhg1(q),q=1,nbgp)

            Do 60 q=1,nbgp
               prhe=0
               Do 61 i=1,Pg
                  prhg(q,i)=prh(i)
                  preg(q,i)=pre(i)
                  prhe=prhe+(prh(i)+pre(i))
  61           Continue
               prheg(q)=prhe
  60        Continue
            Do 600 q=1,nbgp
               Nmhgn=Nmhgn+Nmhg(q-1)*Pgg(q-1)
               Nmegn=Nmegn+Nmeg(q-1)*Pgg(q-1)
               Nmhegn=Nmhgn+Nmegn
               Nmheg=Nmhg(q)+Nmeg(q)
               Do 610 i=1,Nmheg*Pgg(q)
                  Do 620 j=1,Nmheg*Pgg(q)
                     Sg(Nmhegn+i,Nmhegn+j)=S(i,j)
                     RZhegg(Nmhegn+i,Nmhegn+j)=RZheg(i,j)
 620              Continue
 610           Continue
               prhegn=prhegn+prheg(q-1)
               Do 615 i=1,prheg(q)
                  Do 625 j=1,prheg(q)
                     Sprg(prhegn+i,prhegn+j)=Spr(i,j)
                     RZhegprg(prhegn+i,prhegn+j)=RZhegpr(i,j)
 625              Continue
 615           Continue
C le reste de la matrice est compose de zeros
               If (q.GT.1) Then
                   Do 630 i=1,Nmhegn
                      Do 640 j=1,Nmheg*Pgg(q)
                         Sg(i,Nmhegn+j)=(0.,0.)
                         Sg(Nmhegn+j,i)=(0.,0.)
                         RZhegg(i,Nmhegn+j)=(0.,0.)
                         RZhegg(Nmhegn+j,i)=(0.,0.)
 640                  Continue
 630              Continue
                   Do 635 i=1,prhegn
                      Do 645 j=1,prheg(q)
                         Sprg(i,prhegn+j)=(0.,0.)
                         Sprg(prhegn+j,i)=(0.,0.)
                         RZhegprg(i,prhegn+j)=(0.,0.)
                         RZhegprg(prhegn+j,i)=(0.,0.)
 645                  Continue
 635              Continue
               EndIf
 600        Continue
C      write (*,*) 'boucle Nmhg(q)=',(Nmhg(q),q=1,nbgp)

         Else

            Do 700 q=1,nbgp
               Nmh=Nmhg(q)
               Nme= Nmeg(q)
               Pg=Pgg(q)
               a=ag(q)
               Do 705 i=1,Pgg(q)
                  b(i)=bg(q,i)
                  z(i)=zg(q,i)
 705           Continue
      Call SMGrill(intg,ni,epsrel,Nmh,Nme,Pg,a,prh,pre,S,RZheg,
     &Spr,RZhegpr)
               prhe=0
               Do 71 i=1,Pgg(q)
                  prhg(q,i)=prh(i)
                  preg(q,i)=pre(i)
                  prhe=prhe+(prh(i)+pre(i))
  71           Continue
               prheg(q)=prhe
               Nmhgn=Nmhgn+Nmhg(q-1)*Pgg(q-1)
               Nmegn=Nmegn+Nmeg(q-1)*Pgg(q-1)
               Nmhegn=Nmhgn+Nmegn
               Nmheg=Nmhg(q)+Nmeg(q)
               Do 710 i=1,Nmheg*Pgg(q)
                  Do 720 j=1,Nmheg*Pgg(q)
                     Sg(Nmhegn+i,Nmhegn+j)=S(i,j)
                     RZhegg(Nmhegn+i,Nmhegn+j)=RZheg(i,j)
 720              Continue
 710           Continue
               prhegn=prhegn+prheg(q-1)
               Do 715 i=1,prheg(q)
                  Do 725 j=1,prheg(q)
                     Sprg(prhegn+i,prhegn+j)=Spr(i,j)
                     RZhegprg(prhegn+i,prhegn+j)=RZhegpr(i,j)
 725              Continue
 715           Continue
C le reste de la matrice est compose de zeros
               If (q.GT.1) Then
                   Do 730 i=1,Nmhegn
                      Do 740 j=1,Nmheg*Pgg(q)
                         Sg(i,Nmhegn+j)=(0.,0.)
                         Sg(Nmhegn+j,i)=(0.,0.)
                         RZhegg(i,Nmhegn+j)=(0.,0.)
                         RZhegg(Nmhegn+j,i)=(0.,0.)
 740                  Continue
 730              Continue
                   Do 735 i=1,prhegn
                      Do 745 j=1,prheg(q)
                         Sprg(i,prhegn+j)=(0.,0.)
                         Sprg(prhegn+j,i)=(0.,0.)
                         RZhegprg(i,prhegn+j)=(0.,0.)
                         RZhegprg(prhegn+j,i)=(0.,0.)
 745                  Continue
 735              Continue
               EndIf
 700        Continue
         EndIf

      Else

         Do 800 q=1,nbgp
            X0=X0g(q)
            D0=D0g(q)
            If (repgril.EQ.'o'.OR.repgril.EQ.'O') Then
      Call SMGrill(intg,ni,epsrel,Nmh,Nme,Pg,a,prh,pre,S,RZheg,
     &Spr,RZhegpr)
               prhe=0
               Do 81 i=1,Pg
                  prhg(q,i)=prh(i)
                  preg(q,i)=pre(i)
                  prhe=prhe+(prh(i)+pre(i))
  81           Continue
               prheg(q)=prhe
               Nmhgn=Nmhgn+Nmhg(q-1)*Pgg(q-1)
               Nmegn=Nmegn+Nmeg(q-1)*Pgg(q-1)
               Nmhegn=Nmhgn+Nmegn
               Nmheg=Nmhg(q)+Nmeg(q)
               Do 810 i=1,Nmheg*Pgg(q)
                  Do 820 j=1,Nmheg*Pgg(q)
                     Sg(Nmhegn+i,Nmhegn+j)=S(i,j)
                     RZhegg(Nmhegn+i,Nmhegn+j)=RZheg(i,j)
 820              Continue
 810           Continue
               prhegn=prhegn+prheg(q-1)
               Do 815 i=1,prheg(q)
                  Do 825 j=1,prheg(q)
                     Sprg(prhegn+i,prhegn+j)=Spr(i,j)
                     RZhegprg(prhegn+i,prhegn+j)=RZhegpr(i,j)
 825              Continue
 815           Continue
C le reste de la matrice est composes de zeros
               If (q.GT.1) Then
                   Do 830 i=1,Nmhegn
                      Do 840 j=1,Nmheg*Pgg(q)
                         Sg(i,Nmhegn+j)=(0.,0.)
                         Sg(Nmhegn+j,i)=(0.,0.)
                         RZhegg(i,Nmhegn+j)=(0.,0.)
                         RZhegg(Nmhegn+j,i)=(0.,0.)
 840                  Continue
 830              Continue
                   Do 835 i=1,prhegn
                      Do 845 j=1,prheg(q)
                         Sprg(i,prhegn+j)=(0.,0.)
                         Sprg(prhegn+j,i)=(0.,0.)
                         RZhegprg(i,prhegn+j)=(0.,0.)
                         RZhegprg(prhegn+j,i)=(0.,0.)
 845                  Continue
 835              Continue
               EndIf
            Else
               Nmh=Nmhg(q)
               Nme= Nmeg(q)
               Pg=Pgg(q)
               a=ag(q)
               Do 905 i=1,Pgg(q)
                  b(i)=bg(q,i)
                  z(i)=zg(q,i)
 905           Continue
       Call SMGrill(intg,ni,epsrel,Nmh,Nme,Pg,a,prh,pre,S,RZheg,
     &Spr,RZhegpr)
               prhe=0
               Do 91 i=1,Pgg(q)
                  prhg(q,i)=prh(i)
                  preg(q,i)=pre(i)
                  prhe=prhe+(prh(i)+pre(i))
  91           Continue
               prheg(q)=prhe
               Nmhgn=Nmhgn+Nmhg(q-1)*Pgg(q-1)
               Nmegn=Nmegn+Nmeg(q-1)*Pgg(q-1)
               Nmhegn=Nmhgn+Nmegn
               Nmheg=Nmhg(q)+Nmeg(q)
               Do 910 i=1,Nmheg*Pgg(q)
                  Do 920 j=1,Nmheg*Pgg(q)
                     Sg(Nmhegn+i,Nmhegn+j)=S(i,j)
                     RZhegg(Nmhegn+i,Nmhegn+j)=RZheg(i,j)
 920              Continue
 910           Continue
               prhegn=prhegn+prheg(q-1)
               Do 915 i=1,prheg(q)
                  Do 925 j=1,prheg(q)
                     Sprg(prhegn+i,prhegn+j)=Spr(i,j)
                     RZhegprg(prhegn+i,prhegn+j)=RZhegpr(i,j)
 925              Continue
 915           Continue
C le reste de la matrice est composes de zeros
               If (q.GT.1) Then
                   Do 930 i=1,Nmhegn
                      Do 940 j=1,Nmheg*Pgg(q)
                         Sg(i,Nmhegn+j)=(0.,0.)
                         Sg(Nmhegn+j,i)=(0.,0.)
                         RZhegg(i,Nmhegn+j)=(0.,0.)
                         RZhegg(Nmhegn+j,i)=(0.,0.)
 940                  Continue
 930               Continue
                   Do 935 i=1,prhegn
                      Do 945 j=1,prheg(q)
                         Sprg(i,prhegn+j)=(0.,0.)
                         Sprg(prhegn+j,i)=(0.,0.)
                         RZhegprg(i,prhegn+j)=(0.,0.)
                         RZhegprg(prhegn+j,i)=(0.,0.)
 945                  Continue
 935              Continue
               EndIf
            EndIf

 800     Continue

      EndIf

      write (*,*) 'Nmhg(q)=',(Nmhg(q),q=1,nbgp)




C Calcul de la matrice amplitude et phase

      prhegn=prhegn+prheg(nbgp)
      Nmhegn=Nmhegn+(Nmhg(nbgp)+Nmeg(nbgp))*Pgg(nbgp)
      Write (*,*) 'prhegn,Nmhegn=',prhegn,Nmhegn


      If (sortie.EQ.'o'.OR.sortie.EQ.'O') Then
 
      If (intg.EQ.1) Then
         Write (*,*) 'AKF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.2) Then
         Write (*,*) 'ARF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.3) Then
         Write (*,*) 'BDF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.4) Then
         Write (*,*) 'AHF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      Else
         Write (*,*) 'BAF/ni,Nmh,Nme=',ni,Nmh,Nme
      EndIf




CC Matrice de scattering
      Write (*,*) 'MATRICE DE SCATTERING'
      Write (*,*) 'La matrice en Amplitude est :'
      Call CalMatAmplig(prhegn,Sprg,SMAmpli)
      Call ImpMatRAg(prhegn,prhegn,SMAmpli)
      Write (*,*) 'La matrice en Phase est :'
      Call CalMatPhasg(prhegn,Sprg,SMPhas)
      Call ImpMatRPg(prhegn,prhegn,SMPhas)
C      Write (*,*) 'MATRICE D IMPEDANCE'
C      Write (*,*) 'La matrice en Amplitude est :'
C      Call CalMatAmplig(prhegn,RZhegprg,SMAmpli)
C      Call ImpMatRAg(prhegn,prhegn,SMAmpli)
C      Write (*,*) 'La matrice en Phase est :'
C      Call CalMatPhasg(prhegn,RZhegprg,SMPhas)
C      Call ImpMatRPg(prhegn,prhegn,SMPhas)

      EndIf



C      Call ContSymg(Nmhegn,Sg,Statsymgg)
C      write (*,*) 'Statsym OK'
C      Call ContUnig(prhegn,Sprg,Statuni)
C      write (*,*) 'Statuni OK'

C      Call CalMatAmpligg(Nmhegn,Sg,SMAmplig)
C      Call ImpMatRAgg(Nmhegn,Nmhegn,SMAmplig)



C le programme perd le k0g il y a un decalage de nbgp
      Do 101 q=1,nbgp
         k0g(q)=k0
 101  Continue



C Sauvegarde des resultats pour trace du spectre par matlab\
      Open(10,File='ScattSpcg.don',Status='unknown',Form='unformatted')
          Write(10) nbgp,Nmhegn,prhegn
          Write(10) (Pgg(i),i=1,nbgp)
          Write(10) (ag(i),i=1,nbgp)
          Write(10) ((bg(i,j),j=1,Pgg(i)),i=1,nbgp)
          Write(10) ((zg(i,j),j=1,Pgg(i)),i=1,nbgp)
          Write(10) (Nmhg(i),i=1,nbgp)
          Write(10) (Nmeg(i),i=1,nbgp)
          Write(10) ((prhg(i,j),j=1,Pgg(i)),i=1,nbgp)
          Write(10) ((preg(i,j),j=1,Pgg(i)),i=1,nbgp)
C          Write(10) (prheg(i),i=1,nbgp)
          Write(10) (k0g(i),i=1,nbgp)
          Write(10) (X0g(i),i=1,nbgp)
          Write(10) (D0g(i),i=1,nbgp)
       If (reppr.EQ.'o'.OR.reppr.EQ.'O') then
          Write(10) ((Sprg(i,j),i=1,prhegn),j=1,prhegn)
          Write(10) ((RZhegprg(i,j),i=1,prhegn),j=1,prhegn)
       Else
          Write(10) ((Sg(i,j),i=1,Nmhegn),j=1,Nmhegn)
          Write(10) ((RZhegg(i,j),i=1,Nmhegn),j=1,Nmhegn)
       EndIf
      Rewind 10
      Close(10)
      Write (*,*) 'nom du fichier de sauvegarde : ScattSpcg.don '
C      Write (*,*) 'k0g',k0g

      Stop
      
                    
      
C      Return
      End

C Fin du programme principal


C ***
C Surotine de CONTrole de la SYMetrie du module d une mat. de scattering
C ***
      Subroutine ContSymg(Dim,Sg,Statsymgg)

      Integer *4 i,j,Dim
      Integer *4 Statsymgg(2*10*60*5,2*10*60*5)
      Real *8 A21r,A12r,Diffr
      Real *8 A21i,A12i,Diffi
      Complex *16 Sg(2*10*60*5,2*10*60*5)
      Character *1 rep

      Do 100 j=1,Dim
         Do 110 i=1,Dim
            If (i.NE.j) Then
               A12r=DBLE(Sg(i,j))
               A21r=DBLE(Sg(j,i))
               A12i=DIMAG(Sg(i,j))
               A21i=DIMAG(Sg(j,i))
               Diffr=DABS(A12r-A21r)
               Diffi=DABS(A12i-A21i)
               If ((Diffr.LT.0.001).AND.(Diffi.LT.0.001)) Then
                  Statsymgg(i,j)=0
               Else 
C               If ((Diffr.GT.0.001).OR.(Diffi.GT.0.001)) Then
                  Statsymgg(i,j)=-1
C                  Stop 'Probleme de Symetrie de la matice '
               EndIf
            EndIf 
 110     Continue
 100  Continue
  
   1  Return
      End

C ***
C Suroutine de CONTrole de l UNItarite d une matrice de scattering
C ***
      Subroutine ContUnig(Dim,Sprg,Statuni)

      Integer *4 i,j,Dim
      Integer *4 Statuni(10*60*5)
      Real *8 Diff,Somm,Som(10*60*5)
      Complex *16 Sprg(10*60*5,10*60*5)
      Character *1 rep

      Do 100 j=1,Dim
         Somm=0.
         Do 110 i=1,Dim
            Somm=Somm+(CDABS(Sprg(i,j)))**2
 110     Continue
         Som(j)=Somm
C         Diff=DABS(Somm-1)
C         If (Diff.LT.0.01) Then
         If (Somm.LT.1.) Then
            Statuni(i)=0
         Else 
             Statuni(i)=-1
C             Write (*,*) 'Probleme d Unitarite '
         EndIf
 100  Continue
  
C      Write (*,1) (Som(j),j=1,Dim)
C   1  Format(F8.4) 
      If (Dim.EQ.1) Then
         Write (*,1) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.2) Then
         Write (*,2) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.3) Then
         Write (*,3) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.4) Then
         Write (*,4) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.5) Then
         Write (*,5) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.6) Then
         Write (*,6) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.7) Then
         Write (*,7) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.8) Then
         Write (*,8) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.9) Then
         Write (*,9) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.10) Then
         Write (*,10) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.11) Then
         Write (*,11) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.12) Then
         Write (*,12) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.13) Then
         Write (*,13) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.14) Then
         Write (*,14) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.15) Then
         Write (*,15) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.16) Then
         Write (*,16) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.17) Then
         Write (*,17) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.18) Then
         Write (*,18) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.19) Then
         Write (*,19) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.20) Then
         Write (*,20) (Som(j),j=1,Dim)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4))
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4))  
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 

      Return
      End

C ***
C Subroutine de CALcul de la MATrice d AMPLItude d une matrice complexe
C ***
      Subroutine CalMatAmplig(prhe,Sprg,Matamp)

      Integer *4 i,j,prhe
      Real *8 Matamp(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matamp(i,j)=CDABS(Sprg(i,j))
 110     Continue
 100  Continue

      Return
      End
C ***
      Subroutine CalMatAmpligg(prhe,Sg,Matampg)

      Integer *4 i,j,prhe
      Real *8 Matampg(2*10*60*5,2*10*60*5)
      Complex *16 Sg(2*10*60*5,2*10*60*5)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matampg(i,j)=CDABS(Sg(i,j))
 110     Continue
 100  Continue

      Return
      End

C ***
C Subroutine de CALcul de la MATrice de PHASe d une matrice complexe

      Subroutine CalMatPhasg(prhe,Sprg,Matpha)

      Integer *4 i,j,prhe
      Real *8 pReel,pIma,phi,pi
      Real *8 Matpha(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

      pi=4.*ATAN(1.)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            pReel=DBLE(Sprg(i,j))
            pIma=DIMAG(Sprg(i,j))
            If (pReel.EQ.0) Then
               If (pIma.GT.0) Then
                  Matpha(i,j)=90.
               ElseIf (pIma.LT.0) Then
                  Matpha(i,j)=-90.
               Else
                  Matpha(i,j)=0.
               EndIf
            Else
               phi=pIma/pReel
               If (pReel.GT.0) Then
                  Matpha(i,j)=(180./pi)*DATAN(phi)
               Else
                  Matpha(i,j)=180.+(180./pi)*DATAN(phi)
                  If (Matpha(i,j).GT.180.) Then
                     Matpha(i,j)=Matpha(i,j)-360.
                  EndIf
               EndIf
            EndIf 
 110     Continue
 100  Continue

      Return
      End 

C *** 
C Subroutine pour IMPrimer une MATrice Reelle d Amplitude
C ***
      Subroutine ImpMatRAg(lprhe,cprhe,mat)

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

      If (cprhe.EQ.1) Then
         Write (*,1) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4)) 
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4)) 
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 
  21  Format(21(F8.4)) 
  22  Format(22(F8.4)) 
  23  Format(23(F8.4)) 
  24  Format(24(F8.4)) 
  25  Format(25(F8.4)) 
  26  Format(26(F8.4)) 
  27  Format(27(F8.4)) 
  28  Format(28(F8.4)) 
  29  Format(29(F8.4)) 
  30  Format(30(F8.4)) 

      Return
      End
C ***
      Subroutine ImpMatRAgg(lprhe,cprhe,matg)

      Integer *4 i,j,cprhe,lprhe
      Real *8 matg(2*10*60*5,2*10*60*5)

      If (cprhe.EQ.1) Then
         Write (*,1) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4)) 
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4)) 
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 
  21  Format(21(F8.4)) 
  22  Format(22(F8.4)) 
  23  Format(23(F8.4)) 
  24  Format(24(F8.4)) 
  25  Format(25(F8.4)) 
  26  Format(26(F8.4)) 
  27  Format(27(F8.4)) 
  28  Format(28(F8.4)) 
  29  Format(29(F8.4)) 
  30  Format(30(F8.4)) 

      Return
      End

C Subroutine pour IMPrimer une MATrice Reelle de Phase
C ***
      Subroutine ImpMatRPg(lprhe,cprhe,mat)

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

      If (cprhe.EQ.1) Then
         Write (*,1) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.2)) 
   2  Format(2(F8.2)) 
   3  Format(3(F8.2)) 
   4  Format(4(F8.2)) 
   5  Format(5(F8.2)) 
   6  Format(6(F8.2)) 
   7  Format(7(F8.2)) 
   8  Format(8(F8.2)) 
   9  Format(9(F8.2)) 
  10  Format(10(F8.2)) 
  11  Format(11(F8.2)) 
  12  Format(12(F8.2)) 
  13  Format(13(F8.2)) 
  14  Format(14(F8.2)) 
  15  Format(15(F8.2)) 
  16  Format(16(F8.2)) 
  17  Format(17(F8.2)) 
  18  Format(18(F8.2)) 
  19  Format(19(F8.2)) 
  20  Format(20(F8.2)) 
  21  Format(21(F8.2)) 
  22  Format(22(F8.2)) 
  23  Format(23(F8.2)) 
  24  Format(24(F8.2)) 
  25  Format(25(F8.2)) 
  26  Format(26(F8.2)) 
  27  Format(27(F8.2)) 
  28  Format(28(F8.2)) 
  29  Format(29(F8.2)) 
  30  Format(30(F8.2)) 

      Return
      End

