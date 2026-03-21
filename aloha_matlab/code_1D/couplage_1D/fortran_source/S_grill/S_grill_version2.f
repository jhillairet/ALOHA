C ********************************************************************* 
C Programme de calcul de matrice de scattering
C DIScontinuite GRILl PLasma tire du prg bersNSbR2.f
C valide (en spectre) %SWAN le 17 octdbre 1995
C ********************************************************************* 

      Subroutine SMGrill(intg,ni,epsrel,Nmh,Nme,Pg,a,prh,pre,S,RZheg,
     &Spr,RZhegpr)

C Declaration de variable

      Common/nbmod/Nmhm,Nmem
      Common/Airy/app
      Common/err/epsabs
      Common/nlim/nlimit
      Common/Nax/Nmax
      Common/dim/b,z
      Common/idem/bcte
      Common/onde/k0
      Common/plasma/X0,D0

C Nb. max. de guide du Grill = 60
      Integer *4 Pg,ni,m,n,i,j,Intg,k,l,Nmheg,dNmaxg,app
      Integer *4 Nmh,Nme,Nmax,Gmax,maxh,maxe,ifail,nlimit,Nmhm,Nmem
      Integer *4 mhi(2,10*10),mei(2,10*10)
      Integer *4 mhj(2,10*10),mej(2,10*10)
      Integer *4 prhi,prei,nprhi,nprei
      Integer *4 prh(60),pre(60),nprh(60),npre(60),prhe
      Integer *4 Statsym(2*10/2*60,2*10/2*60),Statuni(2*10/2*60)
      Integer *4 Statsymgg(2*10/2*60,2*10/2*60)

      Real *8 pi,k0,f,cl,me,Eps0,qe,dphas
      Real *8 a,ac,b(60),bc,b0,ei,e(60),z(60),pc
      Real *8 n0,nc,Dn,D0,X0
      Real *8 Wkspace(2*10*60),epsrel,epsabs,abserr
      Real *8 SMAmpli(2*10/2*60,2*10/2*60),SMPhas(2*10/2*60,2*10/2*60)

      Complex *16 zohe(10,10),zoeh(10,10)
      Complex *16 Phih(10,10),Phie(10,10)
      Complex *16 PhihI(10,10),PhieI(10,10)
      Complex *16 RZh(10,10),RYh(10,10)
      Complex *16 RZe(10,10),RYe(10,10)
      Complex *16 Id(2*10,2*10)
      Complex *16 PhiheI(2*10,2*10),Phiheg(2*10*60,2*10*60)
      Complex *16 RZhe(2*10,2*10),RZheg(2*10*60,2*10*60)
      Complex *16 RZhegpr(2*10/2*60,2*10/2*60)
      Complex *16 PhiRZhe(2*10,2*10),PhiRZheg(2*10*60,2*10*60)
      Complex *16 Y(2*10,2*10*60),H(2*10*60,2*10*60)
      Complex *16 YRZheg(2*10,2*10*60),RZYRZg(2*10,2*10*60)
      Complex *16 Idg(2*10*60,2*10*60),IdgpH(2*10*60,2*10*60)
      Complex *16 Idg12(2*10*60,2*10*60)
      Complex *16 IdgpHI(2*10*60,2*10*60),IdgmH(2*10*60,2*10*60)
      Complex *16 S(2*10*60,2*10*60),Spr(2*10/2*60,2*10/2*60)
      Complex *16 Alpha,Beta

      Character *1 rep,reps,bcte,TransA,TransB

C Declaration des constantes

      pi=4.*ATAN(1.)
      cl=3.E8
      me=9.1091E-31
      Eps0=8.854E-12
      qe=1.6021E-19
C nombre de mode max
      Nmax=10
C nombre de guide max
      Gmax=60
      Nmhm=Nmh
      Nmem=Nme

C Matrice Indentite Idg((Nmh+Nme)*Pg,(Nmh+Nme)*Pg)
      Do 500 k=1,(Nmh+Nme)*Pg
         Do 510 l=1,(Nmh+Nme)*Pg
            If (k.EQ.l) Then
               Idg(k,l)=(1.,0.)
            Else
               Idg(k,l)=(0.,0.)
            EndIf
 510     Continue
 500  Continue

      TransA='N'
      TransB='N'
      Alpha=(1.,0.)
      Beta=(0.,0.)
      ifail=0
      Nmheg=(Nmh+Nme)*Pg
      dNmaxg=2*Nmax*Gmax

C Calcul de la matrice d'impedance globale
      Call CalMatig(Pg,a,b,Nmh,Nme,RZheg)

C Boucle sur les differents guides
      Do 400 i=1,Pg
C         write (*,*) 'i= ',i

         Call CalRZh(i,a,b,Nmh,RZh,RYh)
         Call CalRZe(i,a,b,Nme,RZe,RYe)
         Call Zero(Nmh,Nme,zohe)
         Call Zero(Nme,Nmh,zoeh)
         Call Mpm(Nmh,Nmh,Nme,Nme,RZh,zohe,zoeh,RZe,RZhe)

         Call CalcouPgi(intg,i,ni,a,epsrel,Nmh,Nme,Pg,Y)
	 

C Calcul de RYhe^-1*(Y*RZheg)=RZhe*(Y*RZheg)=RZYRZg
         Call F06ZAF(TransA,TransB,Nmh+Nme,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,
     &Alpha,Y,2*Nmax,RZheg,2*Nmax*Gmax,Beta,YRZheg,2*Nmax)

         Call F06ZAF(TransA,TransB,Nmh+Nme,(Nmh+Nme)*Pg,Nmh+Nme,
     &Alpha,RZhe,2*Nmax,YRZheg,2*Nmax,Beta,RZYRZg,2*Nmax)


C Constrution de la matrice H
         Do 410 k=1,Nmh+Nme
            Do 420 l=1,(Nmh+Nme)*Pg
               H((Nmh+Nme)*(i-1)+k,l)=RZYRZg(k,l)
 420        Continue
 410     Continue

 400  Continue


C Calcul de la matrice de Scattering S=(I+H)^-1*(I+H)

      TransA='N'
      TransB='N'
      Alpha=(1.,0.)
      Beta=(0.,0.)
      ifail=0
C Calcul de I+H=IdgpH
      Beta=(1.,0.)
      Call F01CWF(TransA,TransB,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,Alpha,Idg,
     &2*Nmax*Gmax,Beta,H,2*Nmax*Gmax,IdgpH,2*Nmax*Gmax,ifail)
      Beta=(0.,0.)

C Calcul de (I-H)=IdgmH
      Beta=(-1.,0.)
      Call F01CWF(TransA,TransB,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,Alpha,Idg,
     &2*Nmax*Gmax,Beta,H,2*Nmax*Gmax,IdgmH,2*Nmax*Gmax,ifail)
      Beta=(0.,0.)
 
      Do 520 k=1,(Nmh+Nme)*Pg
         Do 530 l=1,(Nmh+Nme)*Pg
               Idg12(k,l)=IdgpH(k,l)
 530     Continue
 520  Continue

C Calcul de (I+H)^-1=IdgpHI
      Call F04ADF(IdgpH,2*Nmax*Gmax,Idg,2*Nmax*Gmax,(Nmh+Nme)*Pg,
     &(Nmh+Nme)*Pg,IdgpHI,2*Nmax*Gmax,Wkspace,ifail)

C Calcul de S=(I+H)^-1*(I-H)
      Call F06ZAF(TransA,TransB,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,
     &Alpha,IdgpHI,2*Nmax*Gmax,IdgmH,2*Nmax*Gmax,Beta,S,2*Nmax*Gmax)
         
C Calcul du nombre de modes propageants - construction de la matrice
C le nombre de mode propageant est propre a chque guide
C Verification de l'unitarite et de la symetrie

      Do 600 i=1,Pg
         Call CalNbmpr(i,Nmh,Nme,a,b,mhi,mei,prhi,prei,nprhi,nprei)
         prh(i)=prhi
         pre(i)=prei
         nprh(i)=nprhi
         npre(i)=nprei
 600  Continue

      prhe=0
      Call ConstMatpr(Pg,prh,pre,Nmh,Nme,S,Spr)
      Call ConstMatpr(Pg,prh,pre,Nmh,Nme,RZheg,RZhegpr)
      Do 700 i=1,Pg
         prhe=prhe+(prh(i)+pre(i))
 700  Continue

C Calcul de la matrice amplitude et phase

C Integration 
      If (intg.EQ.1) Then
         Write (*,*) 'AKF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.2) Then
         Write (*,*) 'ARF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.3) Then
         Write (*,*) 'BDF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      ElseIf (intg.EQ.4) Then
C         Write (*,*) 'AHF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
      Else
         Write (*,*) 'BAF/ni,Nmh,Nme=',ni,Nmh,Nme
      EndIf
C      Write (*,*) 'mhd=',(mhd(1,i),i=1,nprhd)
C      Write (*,*) 'nprh=',(nprh(i),i=1,Pg)
C      Write (*,*) 'prh=',(prh(i),i=1,Pg)
C      Write (*,*) 'pre=',(pre(i),i=1,Pg)
C      Write (*,*) 'mh=',(mhi(1,i),i=1,Nmh)
C      Write (*,*) 'nh=',(mhi(2,i),i=1,Nmh)
C      Write (*,*) 'me=',(mei(1,i),i=1,Nme)
C      Write (*,*) 'ne=',(mei(2,i),i=1,Nme)
C      Write (*,*) 'mh=',(mhi(1,i),i=1,nprh(1))
C      Write (*,*) 'nh=',(mhi(2,i),i=1,nprh(1))
C      Write (*,*) 'La matrice en Amplitude calculee :'
C      Call CalMatAmpli(prhe,Spr,SMAmpli)
C      Call ImpMatRA(prhe,prhe,SMAmpli)
C      Write (*,*) 'La matrice en Phase calculee :'
C      Call CalMatPhas(prhe,Spr,SMPhas)
C      Call ImpMatRP(prhe,prhe,SMPhas)

      Call ContSymgg((Nmh+Nme)*Pg,S,Statsymgg)
      Call ContSym1(prhe,Spr,Statsym)
C      write (*,*) 'Statsym OK'
C      Call ContUni1(prhe,Spr,Statuni)
C      write (*,*) 'Statuni OK'

C      Stop
      Return
      End

C Fin du programme principal


C CARACTERISTIQUE DU GRILL

C ***
C Subroutine CLASSant les MODEs par ordre croissant
C ***

      Subroutine ClassModeTE(Nmax,p,a,b,mh,kch)

      Common/nbmod/Nmhm,Nmem
      Integer *4 Nmax,i,ri,m,n,p,Nmhm,Nmem
      Integer *4 mh(2,10*10),IFAIL,IRANK(10*10)
      Real *8 a,b(60),kch(10*10),pi

      pi=4.*ATAN(1.)

      i=1

C      Do 100 m=0,Nmax-1
C         Do 110 n=0,Nmax-1
C            If (m.EQ.0.AND.n.EQ.0) Then
C               Goto 110
C            EndIf

C on force la polarisation en TEm0
      n=0
      Do 100 m=1,Nmhm
         kch(i)=DSQRT(m**2*b(p)/a+n**2*a/b(p))*pi/DSQRT(a*b(p))
         mh(1,i)=m
         mh(2,i)=n
         i=i+1
 100  Continue 
C      ri=i-1
 
C classement des modes se referant au rang de kch
C la troncature des Nmh premier terme ce fera ds les subroutines
C      IFAIL=0
C      Call M01DAF(kch,1,ri,'Ascending',IRANK,IFAIL)
C      Call M01EAF(mh,1,ri,IRANK,IFAIL)
C classement des coupures par ordre croissant
C      Call M01CAF(kch,1,ri,'Ascending',IFAIL) 

      Return
      End

C ***
      Subroutine ClassModeTM(Nmax,p,a,b,me,kce)

      Common/nbmod/Nmhm,Nmem
      Integer *4 Nmax,i,ri,m,n,p,Nmhm,Nmem
      Integer *4 me(2,10*10),IFAIL,IRANK(10*10)
      Real *8 a,b(60),kce(10*10),pi,mn

      pi=4.*ATAN(1.)

      i=1
      mn=Nmem/Nmhm

C      Do 100 m=0,Nmax-1
C         Do 110 n=0,Nmax-1
C            If (m.EQ.0.OR.n.EQ.0) Then
C               Goto 110
C            EndIf

C on ne prend que les TM avec m identique au m du TE
      Do 100 m=1,Nmhm
C         Do 110 n=m,(m+IDNINT(mn)-1)
C En multimode le n peut etre < a m car "a" augmente Ez/Ey#a/b.n/m
         Do 110 n=1,IDNINT(mn)
            kce(i)=DSQRT(m**2*b(p)/a+n**2*a/b(p))*pi/DSQRT(a*b(p))
            me(1,i)=m
            me(2,i)=n
            i=i+1
 110     Continue
 100  Continue 
C      ri=i-1

C classement des modes se referant au rang de kce
C la troncature des Nme premier terme ce fera ds les subroutines
C      IFAIL=0
C      Call M01DAF(kce,1,ri,'Ascending',IRANK,IFAIL)
C      Call M01EAF(me,1,ri,IRANK,IFAIL)
C classement des coupures par ordre croissant
C      Call M01CAF(kce,1,ri,'Ascending',IFAIL) 
       
      Return
      End


C ***
C Subroutine CALculant la Racine des Zimpedances carcteristiques 
C ***
      Subroutine CalRZh(p,a,b,Nmh,RZhi,RYhi)
     
      Common/Nax/Nmax
      Common/onde/k0
      Integer *4 Nmh,Nmax,i,p
      Integer *4 mh(2,10*10)
      Real *8 a,b(60),k0,kch(10*10),pi,Z0
      Complex *16 RZh(10),RYh(10),tampon,tamponbis
      Complex *16 RZhi(10,10),RYhi(10,10)

      pi=4.*ATAN(1.)
      Z0=120.*pi
      cl=3.E8
      Eps0=8.854E-12
      
      Call ClassModeTE(Nmax,p,a,b,mh,kch)

      Do 100 i=1,Nmh
C         RZh(i)=CDSQRT(Z0/CDSQRT((1.,0.)*(1-(kch(i)/k0)**2)))
         tampon=(0.,1.)/((cl**2*Eps0)*CDSQRT((1.,0.)*(kch(i)**2-k0**2)))
	 tamponbis=k0*cl*tampon
	 RZh(i)=CDSQRT(tamponbis)
         RYh(i)=1/RZh(i)
 100  Continue
      
      Call Diago(RZh,Nmh,RZhi)
      Call Diago(RYh,Nmh,RYhi)

      Return
      End

C ***
      Subroutine CalRZe(p,a,b,Nme,RZei,RYei)

      Common/Nax/Nmax
      Common/onde/k0
      Integer *4 Nme,Nmax,i,p
      Integer *4 me(2,10*10)
      Real *8 a,b(60),k0,kce(10*10),pi,Z0
      Complex *16 RZe(10),RYe(10),tampon
      Complex *16 RZei(10,10),RYei(10,10)

      pi=4.*ATAN(1.)
      Z0=120.*pi
      cl=3.E8
      Eps0=8.854E-12

      Call ClassModeTM(Nmax,p,a,b,me,kce)

      Do 100 i=1,Nme
C         RZe(i)=CDSQRT(Z0*CDSQRT((1.,0.)*(1-(kce(i)/k0)**2)))
         tampon=CDSQRT((1.,0.)*(kce(i)**2-k0**2))/((0.,1.)*k0*cl*Eps0) 
	 RZe(i)=CDSQRT(tampon)
         RYe(i)=1/RZe(i)
 100  Continue

      Call Diago(RZe,Nme,RZei)
      Call Diago(RYe,Nme,RYei)

      Return
      End

C ***
C CALcul de la MATrice d'Impedance Global
C ***
      Subroutine CalMatig(Pg,a,b,Nmh,Nme,RZheg)

      Integer *4 Nmh,Nme,p,d,Pg,i,j
      Real *8 a,b(60)
      Complex *16 RZh(10,10),RYh(10,10) 
      Complex *16 RZe(10,10),RYe(10,10) 
      Complex *16 RZheg(2*10*60,2*10*60)

      d=0
      Do 100 p=1,Pg
         Call CalRZh(p,a,b,Nmh,RZh,RYh)
         Call CalRZe(p,a,b,Nme,RZe,RYe)
         Do 200 j=1,Nmh
            Do 250 i=1,Nmh
               If (i.EQ.j) Then
                  RZheg(i+d,j+d)=RZh(i,j)
               EndIf
 250        Continue
 200     Continue
         Do 300 j=1,Nme
            Do 350 i=1,Nme
               If (i.EQ.j) Then
                  RZheg(i+d+Nmh,j+d+Nmh)=RZe(i,j)
               EndIf
 350        Continue
 300     Continue
         d=(Nmh+Nme)*p
 100  Continue
      Do 400 j=1,Pg*(Nmh+Nme)
         Do 450 i=1,Pg*(Nmh+Nme)
C            If (i.NE.j) Then
C               RZheg(i,j)=(0.,0.)
C            EndIf
            If (i.EQ.j) Then
               RZheg(i,j)=RZheg(i,j)
            Else
               RZheg(i,j)=(0.,0.)
            EndIf
 450      Continue
 400   Continue

       Return
       End


C ***
C Sudroutine donnant la normalisation des fct. de forme sin et cos
C ***

C coeff. de normalisation en modE TE 

      Subroutine CE(p,a,b,m,n,E,sg)

      Integer *4 m,n,p
      Real *8 a,b(60),Epsm,Epsn,E
      Character *3 sg

      Call Eps(m,Epsm)
      Call Eps(n,Epsn)

      If (m.EQ.0.AND.n.EQ.0) Then
         Write (*,*) 'Pb m=0 et n=0 en Subroutine CE'
         Stop
      EndIf
      E=(DSQRT(Epsm*Epsn)/DSQRT(m**2*b(p)/a+n**2*a/b(p)))*m*m/a

C      If (sg.EQ.'sgp') Then
C         E=E
C      Else
      If (sg.EQ.'sgm') Then
         E=-1.*E
      EndIf

      Return 
      End

C coeff. de normalisation en moDe TM

      Subroutine CM(p,a,b,m,n,D,sg)

      Integer *4 m,n,p
      Real *8 a,b(60),D
      Character *3 sg

      If (m.EQ.0.OR.n.EQ.0) Then
         Write (*,*) 'Pb m=0 ou n=0 en Subroutine CM'
         Stop
      EndIf
      D=(2./DSQRT(m**2*b(p)/a+n**2*a/b(p)))*m*n/b(p)

C      If (sg.EQ.'sgp') Then
C         D=D
C      Else
      If (sg.EQ.'sgm') Then
         D=-1.*D
      EndIf

      Return
      End

C ***
      Subroutine Eps(m,Epsm)

      Integer *4 m
      Real *8 Epsm

      If (m.EQ.0) Then
         Epsm=1.
      Else 
         Epsm=2.
      EndIf
  
      Return
      End

C ***
C Couplage au plasma
C ***

C ***
C Couplage au plasma
C ***

C integration de 0 a 1
      Complex *16 function f011(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/idem/bcte
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,mn,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f11
      Character *1 scale,bcte

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
       
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)

      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)

      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)

      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)

      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f011=(0.,0.)
C         f011=f11
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

      coupl=(2.,0.)

      f11=(-1.,0.)*Airy*mode*coupl/(t*t)
      f011=f11
 
  10  EndIf

      Return
      End

      Complex *16 function f012(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/idem/bcte
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,mn,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f12
      Character *1 scale,bcte

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
       
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)
       
      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)
       
      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)
       
      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)

      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f012=(0.,0.)
C         f012=f12
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

C      coupl=(2.,0.)*(-(-1.)**m*CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*b(i)))
      coupl=(2.,0.)*(-(-1.)**n*CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*b(i)))

      f12=(-1.,0.)*Airy*mode*coupl/(t*t)
      f012=f12

  10   EndIf

      Return
      End


      Complex *16 function f0131(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f131
      Character *1 scale

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
       
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)

      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)

      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)

      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)

      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f0131=(0.,0.)
C         f0131=f131
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

      coupl=CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*DABS(z(j)-z(i)))

      f131=(-1.,0.)*Airy*mode*coupl/(t*t)
      f0131=f131

  10  EndIf

      Return
      End

      Complex *16 function f0132(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f132
      Character *1 scale

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
       
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)

      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)

      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)

      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)
        
      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f0132=(0.,0.)
C         f0132=f132
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

      coupl=CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*DABS(z(j)-z(i)+b(j)-b(i)))

      f132=(-1.,0.)*Airy*mode*coupl/(t*t)
      f0132=f132

  10  EndIf

      Return
      End

      Complex *16 function f0133(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f133
      Character *1 scale

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
       
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)

      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)

      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)

      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)
        
      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f0133=(0.,0.)
C         f0133=f133
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

      coupl=CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*DABS(z(j)-z(i)+b(j)))

      f133=(-1.,0.)*Airy*mode*coupl/(t*t)
      f0133=f133

  10  EndIf

      Return
      End

      Complex *16 function f0134(t)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,nz,m,n,ifail
      Real *8 k0,X0,D0,v,z(60),b(60),pi,t,pc
      Complex *16 y,ye,Aie,Aied,Ai,Aid,Airy,mode,coupl,f134
      Character *1 scale

      v=(1.-t)/t

      scale='u'
      pi=4.*ATAN(1.)
      y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)
     &*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)
      ye=y*CDEXP((0.,1.)*2.*pi/3.)
   
      ifail=1
      Call S17DGF('d',ye,scale,Aied,nz,ifail)

      ifail=1
      Call S17DGF('f',ye,scale,Aie,nz,ifail)

      ifail=1
      Call S17DGF('d',y,scale,Aid,nz,ifail)

      ifail=1
      Call S17DGF('f',y,scale,Ai,nz,ifail)
        
      pc=1.e-300
      If (CDABS(Ai).LT.pc.OR.CDABS(Aie).LT.pc) Then
         f0134=(0.,0.) 
C         f0134=f134
         GoTo 10
      Else
    
      Airy=(Aid/Ai)-CDEXP((0.,-1.)*4.*pi/3.)*(Aied/Aie)

      mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.)/(
     &((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*
     &((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))
     &*CDEXP((0.,1.)*pi/6.)*(0.,1.)

      coupl=CDEXP((0.,1.)*k0*(-1+(0.,1.)*v)*DABS(z(j)-z(i)-b(i)))

      f134=(-1.,0.)*Airy*mode*coupl/(t*t)
      f0134=f134

  10  EndIf

      Return
      End

C ***
C Suroutine de CALcul de RESIDUs
      Subroutine Calresidu(residu)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 i,j,m,n,ifail,nz
      Real *8 pi,z(60),b(60),k0,X0,D0,pc
      Complex *16 y,Ai,Aid,residu
      Character *1 scale

      pi=4.*ATAN(1.)
      scale='u'

      If ((m.EQ.n).AND.(i.EQ.j)) Then

         y=((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)*
     &(D0/X0)**(2./3.)*(X0-1.)*CDEXP((0.,-1.)*pi/3.)
C      y=(0.,0.)

         ifail=0
         Call S17DGF('d',y,scale,Aid,nz,ifail)

         ifail=0
         Call S17DGF('f',y,scale,Ai,nz,ifail)

         If (n.EQ.0) Then
    
            residu=-2.*pi*k0*b(i)*Aid/Ai/
     &((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(2./3.)
     &*CDEXP((0.,1.)*pi/6.)
         Else
            residu=-pi*k0*b(i)*Aid/Ai/
     &((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(2./3.)
     &*CDEXP((0.,1.)*pi/6.)
         EndIf  

      Else
         residu=(0.,0.)
      EndIf

      Return
      End

    
C ***
C separation partie reelle , partie imaginaire
C integration de 0 a 1
      Real *8 function rf011(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f011
      External f011

      rf011=DBLE(f011(t))

      Return
      End

      Real *8 function if011(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f011
      External f011

      if011=DIMAG(f011(t))

      Return
      End

      Real *8 function rf012(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f012
      External f012

      rf012=DBLE(f012(t))

      Return
      End

      Real *8 function if012(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f012
      External f012

      if012=DIMAG(f012(t))

      Return
      End


      Real *8 function rf0131(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0131
      External f0131

      rf0131=DBLE(f0131(t))

      Return
      End

      Real *8 function if0131(t)

      Real *8 t
      Complex *16 f0131
      External f0131

      if0131=DIMAG(f0131(t))

      Return
      End

      Real *8 function rf0132(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0132
      External f0132

      rf0132=DBLE(f0132(t))

      Return
      End

      Real *8 function if0132(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0132
      External f0132

      if0132=DIMAG(f0132(t))

      Return
      End

      Real *8 function rf0133(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0133
      External f0133

      rf0133=DBLE(f0133(t))

      Return
      End

      Real *8 function if0133(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0133
      External f0133

      if0133=DIMAG(f0133(t))

      Return
      End

      Real *8 function rf0134(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0134
      External f0134

      rf0134=DBLE(f0134(t))

      Return
      End

      Real *8 function if0134(t)
      Common/telnum/knout

      Integer *4 knout
      Real *8 t
      Complex *16 f0134
      External f0134

      if0134=DIMAG(f0134(t))

      Return
      End

   
C ***
C Subroutine d'INTegration


C INTegration de 0 a 1 routine NAG D01AKF
      Subroutine Int01AKF(epsrel,resultr,resulti,absrm)

      Common/err/epsabs
      Common/telnum/knout
      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j

      Integer *4 lw,liw
      parameter (lw=4*500,liw=lw/4)
      Integer *4 iw(liw)
      Integer *4 ifail,knout
      Integer *4 ct,i,j,m,n
      Real *8 w(lw)
      Real *8 d,f,abserr,absr(12),absrm,epsabs,epsrel
      Real *8 b(60),z(60),mn
      Real *8 resultr0112,resultr013,resultr
      Real *8 resulti0112,resulti013,resulti
      Real *8 resultr011,resultr012
      Real *8 resulti011,resulti012
      Real *8 resultr0131,resultr0132,resultr0133,resultr0134
      Real *8 resulti0131,resulti0132,resulti0133,resulti0134
      Complex *16 residu

      Real *8 rf011,if011,rf012,if012
      Real *8 rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134
      External rf011,if011,rf012,if012
      External rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134

C      epsabs=epsrel
C      epsabs=0.
      d=0.
      f=1.
      knout=0
      ifail=-1

      If (i.EQ.j) Then
         mn=(m+n)/2.
         If (mn.EQ.DINT(mn)) Then
         Call D01AKF(rf011,d,f,epsabs,epsrel,resultr011,abserr,
     &w,lw,iw,liw,ifail)
         absr(1)=abserr
         Call D01AKF(if011,d,f,epsabs,epsrel,resulti011,abserr,
     &w,lw,iw,liw,ifail)
         absr(2)=abserr
         Call D01AKF(rf012,d,f,epsabs,epsrel,resultr012,abserr,
     &w,lw,iw,liw,ifail)
         absr(3)=abserr
         Call D01AKF(if012,d,f,epsabs,epsrel,resulti012,abserr,
     &w,lw,iw,liw,ifail)
         absr(4)=abserr
            resultr0112=resultr011+resultr012
            resulti0112=resulti011+resulti012
         Else
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
         EndIf
         resultr013=(0.,0.)
         resulti013=(0.,0.)
         absr(5)=0.
         absr(6)=0.
         absr(7)=0.
         absr(8)=0.
         absr(9)=0.
         absr(10)=0.
         absr(11)=0.
         absr(12)=0.
      Else
         resultr0112=(0.,0.)
         resulti0112=(0.,0.)
         absr(1)=0.
         absr(2)=0.
         absr(3)=0.
         absr(4)=0.
         Call D01AKF(rf0131,d,f,epsabs,epsrel,resultr0131,abserr,
     &w,lw,iw,liw,ifail)
         absr(5)=abserr
         Call D01AKF(if0131,d,f,epsabs,epsrel,resulti0131,abserr,
     &w,lw,iw,liw,ifail)
         absr(6)=abserr
         Call D01AKF(rf0132,d,f,epsabs,epsrel,resultr0132,abserr,
     &w,lw,iw,liw,ifail)
         absr(7)=abserr
         Call D01AKF(if0132,d,f,epsabs,epsrel,resulti0132,abserr,
     &w,lw,iw,liw,ifail)
         absr(8)=abserr
         Call D01AKF(rf0133,d,f,epsabs,epsrel,resultr0133,abserr,
     &w,lw,iw,liw,ifail)
         absr(9)=abserr
         Call D01AKF(if0133,d,f,epsabs,epsrel,resulti0133,abserr,
     &w,lw,iw,liw,ifail)
         absr(10)=abserr
         Call D01AKF(rf0134,d,f,epsabs,epsrel,resultr0134,abserr,
     &w,lw,iw,liw,ifail)   
         absr(11)=abserr
         Call D01AKF(if0134,d,f,epsabs,epsrel,resulti0134,abserr,
     &w,lw,iw,liw,ifail)
         absr(12)=abserr

         resultr013=resultr0131+(-1.)**(m+n)*resultr0132
     &-(-1.)**m*resultr0133-(-1.)**n*resultr0134
         resulti013=resulti0131+(-1.)**(m+n)*resulti0132
     &-(-1.)**m*resulti0133-(-1.)**n*resulti0134

      EndIf

      If ((m.EQ.n).AND.(i.EQ.j)) Then
         Call Calresidu(residu)
      Else
         residu=(0.,0.)
      EndIf

      If (i.EQ.j) Then
         resultr=DBLE(residu)-resultr0112
         resulti=DIMAG(residu)-resulti0112
      Else
         resultr=DBLE(residu)-resultr013
         resulti=DIMAG(residu)-resulti013
      EndIf
   
      absrm=absr(1)
      Do 10 ct=1,12
         If (absrm.LT.absr(ct+1)) Then
            absrm=absr(ct+1)
         EndIf
  10  Continue

      Return
      End

C INTegration de 0 a 1 routine NAG D01BDF
      Subroutine Int01BDF(epsrel,resultr,resulti,absrm)

      Common/err/epsabs
      Common/telnum/knout
      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j

      Integer *4 ifail,knout
      Integer *4 ct,i,j,m,n
      Real *8 d,f,abserr,absr(12),absrm,epsabs,epsrel
      Real *8 b(60),z(60),mn
      Real *8 resultr0112,resultr013,resultr
      Real *8 resulti0112,resulti013,resulti
      Real *8 resultr011,resultr012
      Real *8 resulti011,resulti012
      Real *8 resultr0131,resultr0132,resultr0133,resultr0134
      Real *8 resulti0131,resulti0132,resulti0133,resulti0134
      Complex *16 residu

      Real *8 rf011,if011,rf012,if012
      Real *8 rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134
      External rf011,if011,rf012,if012
      External rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134

C      epsabs=epsrel
C      epsabs=0.
      d=0.
      f=1.
      knout=0
      ifail=-1

      If (i.EQ.j) Then
         mn=(m+n)/2.
         If (mn.EQ.DINT(mn)) Then
      Call D01BDF(rf011,d,f,epsabs,epsrel,resultr011,abserr)
      absr(1)=abserr
      Call D01BDF(if011,d,f,epsabs,epsrel,resulti011,abserr)
      absr(2)=abserr
      Call D01BDF(rf012,d,f,epsabs,epsrel,resultr012,abserr)
      absr(3)=abserr
      Call D01BDF(if012,d,f,epsabs,epsrel,resulti012,abserr)
      absr(4)=abserr
           resultr0112=resultr011+resultr012
            resulti0112=resulti011+resulti012
         Else
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
         EndIf
         resultr013=(0.,0.)
         resulti013=(0.,0.)
         absr(5)=0.
         absr(6)=0.
         absr(7)=0.
         absr(8)=0.
         absr(9)=0.
         absr(10)=0.
         absr(11)=0.
         absr(12)=0.
      Else
         resultr0112=(0.,0.)
         resulti0112=(0.,0.)
         absr(1)=0.
         absr(2)=0.
         absr(3)=0.
         absr(4)=0.
      Call D01BDF(rf0131,d,f,epsabs,epsrel,resultr0131,abserr)
      absr(5)=abserr
      Call D01BDF(if0131,d,f,epsabs,epsrel,resulti0131,abserr)
      absr(6)=abserr
      Call D01BDF(rf0132,d,f,epsabs,epsrel,resultr0132,abserr)
      absr(7)=abserr
      Call D01BDF(if0132,d,f,epsabs,epsrel,resulti0132,abserr)
      absr(8)=abserr
      Call D01BDF(rf0133,d,f,epsabs,epsrel,resultr0133,abserr)
      absr(9)=abserr
      Call D01BDF(if0133,d,f,epsabs,epsrel,resulti0133,abserr)
      absr(10)=abserr
      Call D01BDF(rf0134,d,f,epsabs,epsrel,resultr0134,abserr)
      absr(11)=abserr
      Call D01BDF(if0134,d,f,epsabs,epsrel,resulti0134,abserr)
      absr(12)=abserr

          resultr013=resultr0131+(-1.)**(m+n)*resultr0132
     &-(-1.)**m*resultr0133-(-1.)**n*resultr0134
         resulti013=resulti0131+(-1.)**(m+n)*resulti0132
     &-(-1.)**m*resulti0133-(-1.)**n*resulti0134

      EndIf

      If ((m.EQ.n).AND.(i.EQ.j)) Then
         Call Calresidu(residu)
      Else
         residu=(0.,0.)
      EndIf

      If (i.EQ.j) Then
         resultr=DBLE(residu)-resultr0112
         resulti=DIMAG(residu)-resulti0112
      Else
         resultr=DBLE(residu)-resultr013
         resulti=DIMAG(residu)-resulti013
      EndIf
   
      absrm=absr(1)
      Do 10 ct=1,12
         If (absrm.LT.absr(ct+1)) Then
            absrm=absr(ct+1)
         EndIf
  10  Continue

      Return
      End

C INTegration de 0 a 1 routine NAG D01ARF
      Subroutine Int01ARF(relacc,resultr,resulti,absrm)
      
      Common/err/epsabs
      Common/telnum/knout
      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j

      Integer *4 maxrul
      parameter (maxrul=0)
      Integer *4 ifail,knout,iparm,ni
      Integer *4 ct,i,j,m,n
      Real *8 d,f,absacc,acc,absr(12),absrm,relacc
      Real *8 b(60),z(60),mn,alpha(390)
      Real *8 resultr0112,resultr013,resultr
      Real *8 resulti0112,resulti013,resulti
      Real *8 resultr011,resultr012
      Real *8 resulti011,resulti012
      Real *8 resultr0131,resultr0132,resultr0133,resultr0134
      Real *8 resulti0131,resulti0132,resulti0133,resulti0134
      Complex *16 residu

      Real *8 rf011,if011,rf012,if012
      Real *8 rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134
      External rf011,if011,rf012,if012
      External rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134

C      relacc=0.
C      absacc=0.
C      absacc=relacc
      absacc=epsabs
      d=0.
      f=1.
      iparm=0
      ifail=-1

      If (i.EQ.j) Then
         mn=(m+n)/2.
         If (mn.EQ.DINT(mn)) Then
      Call D01ARF(d,f,rf011,relacc,absacc,maxrul,iparm,acc,resultr011
     &,ni,alpha,ifail)
      absr(1)=acc 
      Call D01ARF(d,f,if011,relacc,absacc,maxrul,iparm,acc,resulti011
     &,ni,alpha,ifail)
      absr(2)=acc 
      Call D01ARF(d,f,rf012,relacc,absacc,maxrul,iparm,acc,resultr012
     &,ni,alpha,ifail)
      absr(3)=acc 
      Call D01ARF(d,f,if012,relacc,absacc,maxrul,iparm,acc,resulti012
     &,ni,alpha,ifail)
      absr(4)=acc 
            resultr0112=resultr011+resultr012
            resulti0112=resulti011+resulti012
         Else
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
         EndIf
         resultr013=(0.,0.)
         resulti013=(0.,0.)
         absr(5)=0.
         absr(6)=0.
         absr(7)=0.
         absr(8)=0.
         absr(9)=0.
         absr(10)=0.
         absr(11)=0.
         absr(12)=0.
      Else
         resultr0112=(0.,0.)
         resulti0112=(0.,0.)
         absr(1)=0.
         absr(2)=0.
         absr(3)=0.
         absr(4)=0.
      Call D01ARF(d,f,rf0131,relacc,absacc,maxrul,iparm,acc,resultr0131
     &,ni,alpha,ifail)
      absr(5)=acc 
      Call D01ARF(d,f,if0131,relacc,absacc,maxrul,iparm,acc,resulti0131
     &,ni,alpha,ifail)
      absr(6)=acc 
      Call D01ARF(d,f,rf0132,relacc,absacc,maxrul,iparm,acc,resultr0132
     &,ni,alpha,ifail)
      absr(7)=acc 
      Call D01ARF(d,f,if0132,relacc,absacc,maxrul,iparm,acc,resulti0132
     &,ni,alpha,ifail)
      absr(8)=acc 
      Call D01ARF(d,f,rf0133,relacc,absacc,maxrul,iparm,acc,resultr0133
     &,ni,alpha,ifail)
      absr(9)=acc 
      Call D01ARF(d,f,if0133,relacc,absacc,maxrul,iparm,acc,resulti0133
     &,ni,alpha,ifail)
      absr(10)=acc 
      Call D01ARF(d,f,rf0134,relacc,absacc,maxrul,iparm,acc,resultr0134
     &,ni,alpha,ifail)
      absr(11)=acc 
      Call D01ARF(d,f,if0134,relacc,absacc,maxrul,iparm,acc,resulti0134
     &,ni,alpha,ifail)
      absr(12)=acc 

         resultr013=resultr0131+(-1.)**(m+n)*resultr0132
     &-(-1.)**m*resultr0133-(-1.)**n*resultr0134
         resulti013=resulti0131+(-1.)**(m+n)*resulti0132
     &-(-1.)**m*resulti0133-(-1.)**n*resulti0134

      EndIf

      If ((m.EQ.n).AND.(i.EQ.j)) Then
         Call Calresidu(residu)
      Else
         residu=(0.,0.)
      EndIf

      If (i.EQ.j) Then
         resultr=DBLE(residu)-resultr0112
         resulti=DIMAG(residu)-resulti0112
      Else
         resultr=DBLE(residu)-resultr013
         resulti=DIMAG(residu)-resulti013
      EndIf
  
      absrm=absr(1)
      Do 10 ct=1,12
         If (absrm.LT.absr(ct+1)) Then
            absrm=absr(ct+1)
         EndIf
  10  Continue

      Return
      End


C INTegration de 0 a 1 routine NAG D01AHF
      Subroutine Int01AHF(epsrel,resultr,resulti,absrm)

      Common/nlim/nlimit
      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j

      Integer *4 ct,i,j,m,n,ifail,knout,ni,nlimit
      Real *8 b(60),z(60),mn
      Real *8 d,f,relr,absr(12),absrm,epsrel
      Real *8 resultr0112,resultr013,resultr
      Real *8 resulti0112,resulti013,resulti
      Real *8 resultr011,resultr012
      Real *8 resulti011,resulti012
      Real *8 resultr0131,resultr0132,resultr0133,resultr0134
      Real *8 resulti0131,resulti0132,resulti0133,resulti0134
      Complex *16 residu

      Real *8 D01AHF,rf011,if011,rf012,if012
      Real *8 rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134
      External rf011,if011,rf012,if012
      External rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134

      d=0.
      f=1.
C      nlimit=0
      ifail=-1

      If (i.EQ.j) Then
         mn=(m+n)/2.
         If (mn.EQ.DINT(mn)) Then
      resultr011=D01AHF(d,f,epsrel,ni,relr,rf011,nlimit,ifail)
      absr(1)=relr
      resulti011=D01AHF(d,f,epsrel,ni,relr,if011,nlimit,ifail)
      absr(2)=relr
      resultr012=D01AHF(d,f,epsrel,ni,relr,rf012,nlimit,ifail)
      absr(3)=relr
      resulti012=D01AHF(d,f,epsrel,ni,relr,if012,nlimit,ifail)
      absr(4)=relr
            resultr0112=resultr011+resultr012
            resulti0112=resulti011+resulti012
         Else
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
         EndIf
         resultr013=(0.,0.)
         resulti013=(0.,0.)
         absr(5)=0.
         absr(6)=0.
         absr(7)=0.
         absr(8)=0.
         absr(9)=0.
         absr(10)=0.
         absr(11)=0.
         absr(12)=0.
      Else
         resultr0112=(0.,0.)
         resulti0112=(0.,0.)
         absr(1)=0.
         absr(2)=0.
         absr(3)=0.
         absr(4)=0.
      resultr0131=D01AHF(d,f,epsrel,ni,relr,rf0131,nlimit,ifail)
      absr(5)=relr
      resulti0131=D01AHF(d,f,epsrel,ni,relr,if0131,nlimit,ifail)
      absr(6)=relr
      resultr0132=D01AHF(d,f,epsrel,ni,relr,rf0132,nlimit,ifail)
      absr(7)=relr
      resulti0132=D01AHF(d,f,epsrel,ni,relr,if0132,nlimit,ifail)
      absr(8)=relr
      resultr0133=D01AHF(d,f,epsrel,ni,relr,rf0133,nlimit,ifail)
      absr(9)=relr
      resulti0133=D01AHF(d,f,epsrel,ni,relr,if0133,nlimit,ifail)
      absr(10)=relr
      resultr0134=D01AHF(d,f,epsrel,ni,relr,rf0134,nlimit,ifail)
      absr(11)=relr
      resulti0134=D01AHF(d,f,epsrel,ni,relr,if0134,nlimit,ifail)
      absr(12)=relr

         resultr013=resultr0131+(-1.)**(m+n)*resultr0132
     &-(-1.)**m*resultr0133-(-1.)**n*resultr0134
         resulti013=resulti0131+(-1.)**(m+n)*resulti0132
     &-(-1.)**m*resulti0133-(-1.)**n*resulti0134

      EndIf

      If ((m.EQ.n).AND.(i.EQ.j)) Then
         Call Calresidu(residu)
      Else
         residu=(0.,0.)
      EndIf

      If (i.EQ.j) Then
         resultr=DBLE(residu)-resultr0112
         resulti=DIMAG(residu)-resulti0112
      Else
         resultr=DBLE(residu)-resultr013
         resulti=DIMAG(residu)-resulti013
      EndIf

      absrm=absr(1)
      Do 10 ct=1,12
         If (absrm.LT.absr(ct+1)) Then
            absrm=absr(ct+1)
         EndIf
  10  Continue

      Return
      End

C INTegration de 0 a 1 routine NAG D01BAF
      Subroutine Int01BAF(ni,resultr,resulti)

      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j

      Integer *4 i,j,m,n,ifail,knout,ni
      Real *8 b(60),z(60),mn
      Real *8 d,f
      Real *8 resultr0112,resultr013,resultr
      Real *8 resulti0112,resulti013,resulti
      Real *8 resultr011,resultr012
      Real *8 resulti011,resulti012
      Real *8 resultr0131,resultr0132,resultr0133,resultr0134
      Real *8 resulti0131,resulti0132,resulti0133,resulti0134
      Complex *16 residu

      Real *8 D01BAF,rf011,if011,rf012,if012
      Real *8 rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134
      External D01BAZ,rf011,if011,rf012,if012
      External rf0131,if0131,rf0132,if0132,rf0133,if0133,rf0134,if0134

C      write (*,*) 'ni= ',ni
      d=0.
      f=1.
      ifail=-1

      If (i.EQ.j) Then
         mn=(m+n)/2.
         If (mn.EQ.DINT(mn)) Then
            resultr011=D01BAF(D01BAZ,d,f,ni,rf011,ifail)
            resulti011=D01BAF(D01BAZ,d,f,ni,if011,ifail)
            resultr012=D01BAF(D01BAZ,d,f,ni,rf012,ifail)
            resulti012=D01BAF(D01BAZ,d,f,ni,if012,ifail)
            resultr0112=resultr011+resultr012
            resulti0112=resulti011+resulti012
         Else
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
         EndIf
         resultr013=(0.,0.)
         resulti013=(0.,0.)
      Else
         resultr0112=(0.,0.)
         resulti0112=(0.,0.)
         resultr0131=D01BAF(D01BAZ,d,f,ni,rf0131,ifail)
         resulti0131=D01BAF(D01BAZ,d,f,ni,if0131,ifail)
         resultr0132=D01BAF(D01BAZ,d,f,ni,rf0132,ifail)
         resulti0132=D01BAF(D01BAZ,d,f,ni,if0132,ifail)
         resultr0133=D01BAF(D01BAZ,d,f,ni,rf0133,ifail)
         resulti0133=D01BAF(D01BAZ,d,f,ni,if0133,ifail)
         resultr0134=D01BAF(D01BAZ,d,f,ni,rf0134,ifail)
         resulti0134=D01BAF(D01BAZ,d,f,ni,if0134,ifail)

         resultr013=resultr0131+(-1.)**(m+n)*resultr0132
     &-(-1.)**m*resultr0133-(-1.)**n*resultr0134
         resulti013=resulti0131+(-1.)**(m+n)*resulti0132
     &-(-1.)**m*resulti0133-(-1.)**n*resulti0134
      EndIf

      If ((m.EQ.n).AND.(i.EQ.j)) Then
         Call Calresidu(residu)
      Else
         residu=(0.,0.)
      EndIf

      If (i.EQ.j) Then
         resultr=DBLE(residu)-resultr0112
         resulti=DIMAG(residu)-resulti0112
      Else
         resultr=DBLE(residu)-resultr013
         resulti=DIMAG(residu)-resulti013
      EndIf
   
      Return
      End


C ***
C Suroutine de CALcul de RESIDUs en dimension y
      Subroutine Calresiduy(a,mi,mj,residuy)

      Common/guid/i,j
      Common/onde/k0

      Integer *4 i,j,mi,mj
      Real *8 pi,a,k0,residuy

      pi=4.*ATAN(1.)

      If (mi.EQ.mj) Then

      If (mi.NE.0) Then
         residuy=pi*k0*a/((mi*pi/(k0*a))**2.)
      Else
         write (*,*) 'pb mi=mj=0'
         stop
C         residuy=pi*(k0*a)**3./3.
      EndIf

      Else
         residuy=(0.,0.)
      EndIf
   
      Return
      End
         
C ***
C Subroutine de CALcul de la matrice de COUplage Plasma/Guidei
C ***
      Subroutine CalcouPgi(intg,p,ni,a,epsrel,Nmh,Nme,Pg,Y)

      Common/telnum/knout
      Common/Nax/Nmax
      Common/dim/b,z
      Common/mod/m,n
      Common/guid/i,j
      Common/onde/k0
      Common/plasma/X0,D0

      Integer *4 knout,ct,intg
      Integer *4 Nmax,i,j,p,k,l,Pg
      Integer *4 Nmh,Nme,m,n,mi,mj,ma,na,ni
      Integer *4 mhi(2,10*10),mhj(2,10*10)
      Integer *4 mei(2,10*10),mej(2,10*10)
      Real *8 a,b(60),z(60),k0,X0,D0,Z0,Y0
      Real *8 pi,kch(10*10),kce(10*10)
      Real *8 absrr(10*10*4),absrm,epsrel
      Real *8 resultrhh,resultihh,resultrhe,resultihe
      Real *8 resultreh,resultieh,resultree,resultiee
      Real *8 cte,cofm,Di,Dj,Ei,Ej
      Complex *16 Y(2*10,2*10*60)
      Character *3 sgp,sgm

      pi=4.*ATAN(1.)
      Z0=120.*pi
      Y0=1./Z0
      cte=-Y0*(X0/D0)**(1./3.)*(pi/a)**2./k0**4./((2*pi)**2.)

      sgp='sgp'
      sgm='sgm'                 

      j=p
      Do 100 i=1,Pg
cc         write (*,*) 'j,i= ',j,i
         ct=1

         Call ClassModeTE(Nmax,i,a,b,mhi,kch)
         Call ClassModeTE(Nmax,j,a,b,mhj,kch)
         Call ClassModeTM(Nmax,i,a,b,mei,kce)
         Call ClassModeTM(Nmax,j,a,b,mej,kce)

         Do 200 k=1,Nmh
            Do 210 l=1,Nmh
               Call CE(j,a,b,mhj(1,k),mhj(2,k),Ej,sgp)
               Call CE(i,a,b,mhi(1,l),mhi(2,l),Ei,sgm)
               mj=mhj(1,k)
               mi=mhi(1,l)
               If (mj.EQ.mi) Then
                  Call Calresiduy(a,mi,mj,cofm)
               Else
                  cofm=0.
               EndIf
               m=mhj(2,k)
               n=mhi(2,l)
               If (cofm.NE.0.) Then
                  If (Intg.EQ.1) Then
                     Call Int01AKF(epsrel,resultrhh,resultihh,absrm)
                  ElseIf (intg.EQ.2) Then
                     Call Int01ARF(epsrel,resultrhh,resultihh,absrm)
                  ElseIf (intg.EQ.3) Then
                    Call Int01BDF(epsrel,resultrhh,resultihh,absrm)
                  ElseIf (intg.EQ.4) Then
                     Call Int01AHF(epsrel,resultrhh,resultihh,absrm)
                  Else
                     Call Int01BAF(ni,resultrhh,resultihh)
                  EndIf
                 EndIf
               Y(k,l+(Nmh+Nme)*(i-1))=
     &Ej*Ei*cofm*cte*(resultrhh+(0.,1.)*resultihh)
                absrr(ct)=absrm
                ct=ct+1
 210        Continue
 200     Continue

         Do 300 k=1,Nmh
            Do 310 l=1,Nme
               Call CE(j,a,b,mhj(1,k),mhj(2,k),Ej,sgp)
               Call CM(i,a,b,mei(1,l),mei(2,l),Di,sgm)
               mj=mhj(1,k)
               mi=mei(1,l)
               If (mj.EQ.mi) Then
                  Call Calresiduy(a,mi,mj,cofm)
               Else
                  cofm=0.
               EndIf
               m=mhj(2,k)
               n=mei(2,l)
               If (cofm.NE.0.) Then
                  If (Intg.EQ.1) Then
                     Call Int01AKF(epsrel,resultrhe,resultihe,absrm)
                  ElseIf (intg.EQ.2) Then
                     Call Int01ARF(epsrel,resultrhe,resultihe,absrm)
                  ElseIf (intg.EQ.3) Then
                    Call Int01BDF(epsrel,resultrhe,resultihe,absrm)
                  ElseIf (intg.EQ.4) Then
                     Call Int01AHF(epsrel,resultrhe,resultihe,absrm)
                  Else
                     Call Int01BAF(ni,resultrhe,resultihe)
                  EndIf
               EndIf
               Y(k,l+Nmh+(Nmh+Nme)*(i-1))=
     &Ej*Di*cofm*cte*(resultrhe+(0.,1.)*resultihe)
                absrr(ct)=absrm
                ct=ct+1
 310        Continue
 300     Continue

         Do 400 l=1,Nme
            Do 410 k=1,Nmh
               Call CM(j,a,b,mej(1,l),mej(2,l),Dj,sgp)
               Call CE(i,a,b,mhi(1,k),mhi(2,k),Ei,sgm)
               mj=mej(1,l)
               mi=mhi(1,k)
               If (mj.EQ.mi) Then
                  Call Calresiduy(a,mi,mj,cofm)
               Else
                  cofm=0.
               EndIf
               m=mej(2,l)
               n=mhi(2,k)
               If (cofm.NE.0.) Then
                  If (Intg.EQ.1) Then
                     Call Int01AKF(epsrel,resultreh,resultieh,absrm)
                  ElseIf (intg.EQ.2) Then
                     Call Int01ARF(epsrel,resultreh,resultieh,absrm)
                  ElseIf (intg.EQ.3) Then
                    Call Int01BDF(epsrel,resultreh,resultieh,absrm)
                  ElseIf (intg.EQ.4) Then
                     Call Int01AHF(epsrel,resultreh,resultieh,absrm)
                  Else
                     Call Int01BAF(ni,resultreh,resultieh)
                  EndIf
               EndIf
               Y(Nmh+l,k+(Nmh+Nme)*(i-1))=
     &Dj*Ei*cofm*cte*(resultreh+(0.,1.)*resultieh)
                absrr(ct)=absrm
                ct=ct+1
 410        Continue
 400     Continue

         Do 500 l=1,Nme
            Do 510 k=1,Nme
               Call CM(j,a,b,mej(1,l),mej(2,l),Dj,sgp)
               Call CM(i,a,b,mei(1,k),mei(2,k),Di,sgm)
               mj=mej(1,l)
               mi=mei(1,k)
               If (mj.EQ.mi) Then
                  Call Calresiduy(a,mi,mj,cofm)
               Else
                  cofm=0.
               EndIf
               m=mej(2,l)
               n=mei(2,k)
               If (cofm.NE.0.) Then
                  If (Intg.EQ.1) Then
                     Call Int01AKF(epsrel,resultree,resultiee,absrm)
                  ElseIf (intg.EQ.2) Then
                     Call Int01ARF(epsrel,resultree,resultiee,absrm)
                  ElseIf (intg.EQ.3) Then
                    Call Int01BDF(epsrel,resultree,resultiee,absrm)
                  ElseIf (intg.EQ.4) Then
                     Call Int01AHF(epsrel,resultree,resultiee,absrm)
                  Else
                     Call Int01BAF(ni,resultree,resultiee)
                  EndIf
               EndIf
               Y(Nmh+l,k+Nmh+(Nmh+Nme)*(i-1))=
     &Dj*Di*cofm*cte*(resultree+(0.,1.)*resultiee)
                absrr(ct)=absrm
                ct=ct+1
 510        Continue
 500     Continue

C      ct=ct-1
      Call tri(ct,absrr,absrm)
CC      Write (*,*) 'absrm=',absrm

 100  Continue

      Return
      End


      Subroutine tri(ct,absrr,absrm)

      Integer *4 i,ct
      Real *8 absrr(10*10*4),absrm

      absrm=absrr(1)
      Do 10 i=1,ct
         If (absrm.LT.absrr(i+1)) Then
            absrm=absrr(i+1)
         EndIf
  10  Continue

      Return
      End



C ***
C Subroutine creant une matrice DIAGOnale a partir d un vect. colonne
C ***
      Subroutine Diago(Mc,Nm,Md)

      Integer *4 Nm,i,j
      Complex *16 Mc(10),Md(10,10)

      Do 100 j=1,Nm
         Do 110 i=1,Nm
            If (i.EQ.j) Then
               Md(i,j)=Mc(i)
            Else
               Md(i,j)=(0.,0.)
            EndIf
 110     Continue 
 100  Continue

      Return
      End

C ***
C Subroutine donnant une matrice de ZEROs 
C ***
      Subroutine Zero(Mm,Nm,zox)

      Integer *4 i,j,Mm,Nm
      Complex *16 zox(10,10)

      Do 100 j=1,Nm
         Do 110 i=1,Mm
            zox(i,j)=(0.,0.)
 110     Continue
 100  Continue

      Return
      End


C ***
C Subroutine de construction d une Matrice Par Morceau 
C ***
      Subroutine Mpm(Nm11,Nm12,Nm21,Nm22,Mor11,Mor12,Mor21,Mor22,Mat)

      Integer *4 i,j,Nm11,Nm12,Nm21,Nm22
      Complex *16 Mor11(10,10),Mor12(10,10)
      Complex *16 Mor21(10,10),Mor22(10,10)
      Complex *16 Mat(2*10,2*10)

      Do 100 j=1,Nm11
         Do 110 i=1,Nm12
            Mat(i,j)=Mor11(i,j)
 110     Continue
 100  Continue

      Do 120 j=1,Nm22
         Do 130 i=1,Nm12 
            Mat(i,Nm11+j)=Mor12(i,j)
 130     Continue
 120  Continue

      Do 140 j=1,Nm11
         Do 150 i=1,Nm21
            Mat(Nm12+i,j)=Mor21(i,j)
 150     Continue
 140  Continue

      Do 160 j=1,Nm22
         Do 170 i=1,Nm21
            Mat(Nm12+i,Nm11+j)=Mor22(i,j)
 170     Continue
 160  Continue

      Return
      End

C ***
C Subroutine de CALcul  du NomBre de Mode PRopageant
C ***
      Subroutine CalNbmpr(i,Nmh,Nme,a,b,mhi,mei,prhi,prei,nprhi,nprei)

      Common/Nax/Nmax
      Common/onde/k0
      Integer *4 i,l,Nmax,Nmh,Nme,prhi,prei,nprhi,nprei
      Integer *4 mhi(2,10*10),mei(2,10*10)
      Real *8 a,b(60)
      Real *8 k0,kch(10*10),kce(10*10)

      Call ClassModeTE(Nmax,i,a,b,mhi,kch)
      prhi=0
      Do 100 l=1,Nmh
         If (k0.GT.kch(l)) Then
            prhi=prhi+1
         EndIf
 100  Continue
      nprhi=0
      Do 200 l=1,Nmax
         If (k0.GT.kch(l)) Then
            nprhi=nprhi+1
         EndIf
 200  Continue

      Call ClassModeTM(Nmax,i,a,b,mei,kce)
      prei=0
      Do 110 l=1,Nme
         If (k0.GT.kce(l)) Then
            prei=prei+1
         EndIf
 110  Continue
      nprei=0
      Do 210 l=1,Nmax
         If (k0.GT.kce(l)) Then
            nprei=nprei+1
         EndIf
 210  Continue

      Return
      End

C ***
C Subroutine de CONSTruction de la MATrice des modes PRopageants
C ***
      Subroutine ConstMatpr(Pg,prh,pre,Nmh,Nme,S,Spr)

      Integer *4 p,q,i,j,Pg,deci,Ndeci,decj,Ndecj
      Integer *4 Nmh,Nme,prh(60),pre(60)
      Complex *16 S(2*10*60,2*10*60),Spr(2*10/2*60,2*10/2*60)

C Construction de Sprhh
      deci=0
      Ndeci=0
      Do 10 p=1,Pg
         decj=0
         Ndecj=0
         Do 20 q=1,Pg
            Do 30 i=1,prh(p)
               Do 40 j=1,prh(q)
                  Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
  40           Continue
  30        Continue
            decj=decj+prh(q)+pre(q)
            Ndecj=Ndecj+Nmh+Nme
  20     Continue
         deci=deci+prh(p)+pre(p)
         Ndeci=Ndeci+Nmh+Nme
  10  Continue

C Construction de Sprhe
      deci=0
      Ndeci=0
      Do 11 p=1,Pg
         decj=prh(1)
         Ndecj=Nmh
         Do 21 q=1,Pg
            Do 31 i=1,prh(p)
               Do 41 j=1,pre(q)
                  Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
  41           Continue
  31        Continue
            decj=decj+pre(q)+prh(q+1)
            Ndecj=Ndecj+Nme+Nmh
  21     Continue
         deci=deci+prh(p)+pre(p)
         Ndeci=Ndeci+Nmh+Nme
  11  Continue

C Construction de Spreh
      deci=prh(1)
      Ndeci=Nmh
      Do 12 p=1,Pg
         decj=0
         Ndecj=0
         Do 22 q=1,Pg
            Do 32 i=1,pre(p)
               Do 42 j=1,prh(q)
                  Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
  42           Continue
  32        Continue
            decj=decj+prh(q)+pre(q)
            Ndecj=Ndecj+Nmh+Nme
  22     Continue
         deci=deci+pre(p)+prh(p+1)
         Ndeci=Ndeci+Nme+Nmh
  12  Continue

C Construction de Spree
      deci=prh(1)
      Ndeci=Nmh
      Do 13 p=1,Pg
         decj=prh(1)
         Ndecj=Nmh
         Do 23 q=1,Pg
            Do 33 i=1,pre(p)
               Do 43 j=1,pre(q)
                  Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
  43           Continue
  33        Continue
            decj=decj+pre(q)+prh(q+1)
            Ndecj=Ndecj+Nme+Nmh
  23     Continue
         deci=deci+pre(p)+prh(p+1)
         Ndeci=Ndeci+Nme+Nmh
  13  Continue

      Return
      End

C ***
C Surotine de CONTrole de la SYMetrie du module d une mat. de scattering
C ***
      Subroutine ContSym1(Dim,Spr,Statsym)

      Integer *4 i,j,Dim
      Integer *4 Statsym(2*10/2*60,2*10/2*60)
      Real *8 A21,A12,Diff
      Real *8 A21r,A12r,Diffr
      Real *8 A21i,A12i,Diffi
      Complex *16 Spr(2*10/2*60,2*10/2*60)
      Character *1 rep

      Do 100 j=1,Dim
         Do 110 i=1,Dim
            If (i.NE.j) Then
C               A12=CDABS(Spr(i,j))
C               A21=CDABS(Spr(j,i))
C               Diff=DABS(A12-A21)
               A12r=DBLE(Spr(i,j))
               A21r=DBLE(Spr(j,i))
               A12i=DIMAG(Spr(i,j))
               A21i=DIMAG(Spr(j,i))
               Diffr=DABS(A12r-A21r)
               Diffi=DABS(A12i-A21i)
C               If (Diff.LT.0.001) Then
               If ((Diffr.LT.0.001).AND.(Diffi.LT.0.001)) Then
                  Statsym(i,j)=0
               Else 
C               If ((Diffr.GT.0.001).OR.(Diffi.GT.0.001)) Then
C               If (Diff.GT.0.001) Then
                  Statsym(i,j)=-1
C                  Stop 'Probleme de Symetrie de la matice '
               EndIf
            EndIf 
 110     Continue
 100  Continue
  
   1  Return
      End

C ***
      Subroutine ContSymgg(Dim,S,Statsymgg)

      Integer *4 i,j,Dim
      Integer *4 Statsymgg(2*10*60,2*10*60)
      Real *8 A21r,A12r,Diffr
      Real *8 A21i,A12i,Diffi
      Complex *16 S(2*10*60,2*10*60)
      Character *1 rep

      Do 100 j=1,Dim
         Do 110 i=1,Dim
            If (i.NE.j) Then
               A12r=DBLE(S(i,j))
               A21r=DBLE(S(j,i))
               A12i=DIMAG(S(i,j))
               A21i=DIMAG(S(j,i))
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
      Subroutine ContUni1(Dim,Spr,Statuni)

      Integer *4 i,j,Dim
      Integer *4 Statuni(2*10/2*60)
      Real *8 Diff,Somm,Som(2*10/2*60)
      Complex *16 Spr(2*10/2*60,2*10/2*60)
      Character *1 rep

      Do 100 j=1,Dim
         Somm=0.
         Do 110 i=1,Dim
            Somm=Somm+(CDABS(Spr(i,j)))**2
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

      Return
      End

C ***
C Subroutine de CALcul de la MATrice d AMPLItude d une matrice complexe
C ***
      Subroutine CalMatAmpli(prhe,Spr,Matamp)

      Integer *4 i,j,prhe
      Real *8 Matamp(2*10/2*60,2*10/2*60)
      Complex *16 Spr(2*10/2*60,2*10/2*60)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matamp(i,j)=CDABS(Spr(i,j))
 110     Continue
 100  Continue

      Return
      End

C ***
C Subroutine de CALcul de la MATrice de PHASe d une matrice complexe

      Subroutine CalMatPhas(prhe,Spr,Matpha)

      Integer *4 i,j,prhe
      Real *8 pReel,pIma,phi,pi
      Real *8 Matpha(2*10/2*60,2*10/2*60)
      Complex *16 Spr(2*10/2*60,2*10/2*60)

      pi=4.*ATAN(1.)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            pReel=DBLE(Spr(i,j))
            pIma=DIMAG(Spr(i,j))
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
      Subroutine ImpMatRA(lprhe,cprhe,mat)

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60,2*10/2*60)

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

      Return
      End

C Subroutine pour IMPrimer une MATrice Reelle de Phase
C ***
      Subroutine ImpMatRP(lprhe,cprhe,mat)

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60,2*10/2*60)

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

      Return
      End


C fin
