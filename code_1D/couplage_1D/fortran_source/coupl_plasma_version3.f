C ********************************************************************* 
C Programme de calcul de matrice de scattering
C Couplage Grill Plasma avec differents cas poloidaux
C ********************************************************************* 


      Common/err/epsabs
      Common/nlim/nlimit
      Common/dim/b,z
      Common/idem/bcte
      Common/onde/k0
      Common/plasma/X0,D0
      Common/integration/max_nz

cccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 nbgp,i,j,q,intg,ni,nlimit
      Integer *4 T_grill,D_guide_max
      Integer *4 Pg,Nmh,Nme,Nmhgn,Nmegn,Nmheg,Nmhegn,prhe,prhegn
      Integer *4 Pgg(5),Nmhg(5),Nmeg(5),prheg(5),Nmhg1(5)
      Integer *4 prh(60),pre(60),prhg(5,60),preg(5,60)

      Real *8 pi,k0,f,cl,me,Eps0,qe,epsrel,epsabs
      Real *8 a,ac,b(60),bc,b0,ei,e(60),z(60),pc
      Real *8 n0,nc,Dn,D0,X0,max_nz
      Real *8 k0g(5),D0g(5),X0g(5)
      Real *8 ag(5),bg(5,60),eg(5,60),zg(5,60)

      Character *1 repspc,repgril,reppr,sortie,bcte,rep


      Real *8 SMAmpli(2*10/2*60*5,2*10/2*60*5)
      Real *8 SMPhas(2*10/2*60*5,2*10/2*60*5)
      Real *8 SMAmplig(2*10*60*5,2*10*60*5)

      Complex *16 S(2*10*60,2*10*60),Spr(2*10/2*60,2*10/2*60)
      Complex *16 K_cpl(2*10*60,2*10*60)
      Complex *16 Sg(2*10*60*5,2*10*60*5),Sprg(10*60*5,10*60*5)
      Complex *16 RZheg(2*10*60,2*10*60),RZhegpr(2*10/2*60,2*10/2*60)
      Complex *16 RZhegg(2*10*60*5,2*10*60*5)
      Complex *16 RZhegprg(10*60*5,10*60*5)

c&&&      Integer *4 nbgp,i,j,q,intg,ni,nlimit
c&&&      Integer *4 T_grill,D_guide_max
c&&&      Integer *4 Pg,Nmh,Nme,Nmhgn,Nmegn,Nmheg,Nmhegn,prhe,prhegn
c&&&      Integer *4 Pgg(5),Nmhg(5),Nmeg(5),prheg(5),Nmhg1(5)
c&&&      Integer *4 prh(15),pre(15),prhg(5,15),preg(5,15)
c&&&
c&&&      Real *8 pi,k0,f,cl,me,Eps0,qe,epsrel,epsabs
c&&&      Real *8 a,ac,b(15),bc,b0,ei,e(15),z(15),pc
c&&&      Real *8 n0,nc,Dn,D0,X0,max_nz
c&&&      Real *8 k0g(5),D0g(5),X0g(5)
c&&&      Real *8 ag(5),bg(5,15),eg(5,15),zg(5,15)
c&&&
c&&&      Character *1 repspc,repgril,reppr,sortie,bcte,rep
c&&&
c&&&
c&&&      Real *8 SMAmpli(2*3/2*15*5,2*3/2*15*5)
c&&&      Real *8 SMPhas(2*3/2*15*5,2*3/2*15*5)
c&&&      Real *8 SMAmplig(2*3*15*5,2*3*15*5)
c&&&
c&&&      Complex *16 S(2*3*15,2*3*15),Spr(2*3/2*15,2*3/2*15)
c&&&      Complex *16 K_cpl(2*3*15,2*3*15)
c&&&      Complex *16 Sg(2*3*15*5,2*3*15*5),Sprg(3*15*5,3*15*5)
c&&&      Complex *16 RZheg(2*3*15,2*3*15),RZhegpr(2*3/2*15,2*3/2*15)
c&&&      Complex *16 RZhegg(2*3*15*5,2*3*15*5)
c&&&      Complex *16 RZhegprg(3*15*5,3*15*5)


ccccccccccccccccccccccccc



C Declaration des constantes

      pi=4.*ATAN(1.)
      cl=3.E8
      me=9.1091E-31
      Eps0=8.854E-12
      qe=1.6021E-19
      pc=1.*1.E-3

      Open(10,File='par_grill.dat',Status='unknown',Form='formatted')

      Read(10,*) Nmh
      Read(10,*) Nme
      Read(10,*) f
      Read(10,*) n0
      Read(10,*) Dn
      Read(10,*) Pg
      Read(10,*) a
      Read(10,*) (b(i),i=1,Pg)
      Read(10,*) (z(i),i=1,Pg)
      Read(10,*) T_grill
      Read(10,*) D_guide_max
      Read(10,*) epsrel
      Read(10,*) max_nz
      Rewind 10
      Close(10)

      k0=2.*pi*f/cl  
      D0=k0*n0/Dn
      nc=(2.*pi*f)**2*me*Eps0/qe**2
      X0=n0/nc

      
      
      intg=4
      nlimit=0
      ni = 32 
      epsabs=0.


      Call SMGrill(intg,ni,epsrel,T_grill,D_guide_max,Nmh,Nme,Pg,a,
     &prh,pre,S,RZheg,Spr,RZhegpr,K_cpl)
     

C Sauvegarde des resultats pour matlab

      Open(10,File='S_plasma.dat',Status='unknown',Form='unformatted') 
  
      Write(10) ((S(i,j),i=1,(Nme+Nmh)*Pg),j=1,(Nme+Nmh)*Pg)
      Write(10) ((RZheg(i,j),i=1,(Nme+Nmh)*Pg),j=1,(Nme+Nmh)*Pg)
      Write(10) ((K_cpl(i,j),i=1,(Nme+Nmh)*Pg),j=1,(Nme+Nmh)*Pg)

      Rewind 10
      Close(10)
    
      Stop
      
      End

C Fin du programme principal


C ***
C Surotine de CONTrole de la SYMetrie du module d une mat. de scattering
C ***
      Subroutine ContSymg(Dim,Sg,Statsymgg)


ccccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,Dim
      Integer *4 Statsymgg(2*10*60*5,2*10*60*5)
      Real *8 A21r,A12r,Diffr
      Real *8 A21i,A12i,Diffi
      Complex *16 Sg(2*10*60*5,2*10*60*5)
      Character *1 rep

c&&&      Integer *4 i,j,Dim
c&&&      Integer *4 Statsymgg(2*3*15*5,2*3*15*5)
c&&&      Real *8 A21r,A12r,Diffr
c&&&      Real *8 A21i,A12i,Diffi
c&&&      Complex *16 Sg(2*3*15*5,2*3*15*5)
c&&&      Character *1 rep

cccccccccccccccccccccccccccc

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

cccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,Dim
      Integer *4 Statuni(10*60*5)
      Real *8 Diff,Somm,Som(10*60*5)
      Complex *16 Sprg(10*60*5,10*60*5)
      Character *1 rep

c&&&      Integer *4 i,j,Dim
c&&&      Integer *4 Statuni(3*15*5)
c&&&      Real *8 Diff,Somm,Som(3*15*5)
c&&&      Complex *16 Sprg(3*15*5,3*15*5)
c&&&      Character *1 rep

ccccccccccccccccccccccccccccc

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


cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 Matamp(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

c&&&      Integer *4 i,j,prhe
c&&&      Real *8 Matamp(2*3/2*15*5,2*3/2*15*5)
c&&&      Complex *16 Sprg(2*3/2*15*5,2*3/2*15*5)

cccccccccccccccccccccccccccc

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matamp(i,j)=CDABS(Sprg(i,j))
 110     Continue
 100  Continue

      Return
      End
C ***
      Subroutine CalMatAmpligg(prhe,Sg,Matampg)

cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 Matampg(2*10*60*5,2*10*60*5)
      Complex *16 Sg(2*10*60*5,2*10*60*5)

c&&&      Integer *4 i,j,prhe
c&&&      Real *8 Matampg(2*3*15*5,2*3*15*5)
c&&&      Complex *16 Sg(2*3*15*5,2*3*15*5)


ccccccccccccccccccccccccccccc

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

cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 pReel,pIma,phi,pi
      Real *8 Matpha(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

c&&&      Integer *4 i,j,prhe
c&&&      Real *8 pReel,pIma,phi,pi
c&&&      Real *8 Matpha(2*3/2*15*5,2*3/2*15*5)
c&&&      Complex *16 Sprg(2*3/2*15*5,2*3/2*15*5)

cccccccccccccccccccccccccc

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


cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

c&&&      Integer *4 i,j,cprhe,lprhe
c&&&      Real *8 mat(2*3/2*15*5,2*3/2*15*5)

ccccccccccccccccccccccccccccc

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


cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 matg(2*10*60*5,2*10*60*5)

c&&&      Integer *4 i,j,cprhe,lprhe
c&&&      Real *8 matg(2*3*15*5,2*3*15*5)

ccccccccccccccccccccccccccccccc

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

cccccccccccccccccccccccccccccc
c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

c&&&      Integer *4 i,j,cprhe,lprhe
c&&&      Real *8 mat(2*3/2*15*5,2*3/2*15*5)

cccccccccccccccccccccccccccccc

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

