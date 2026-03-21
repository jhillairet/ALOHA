
c **************************************************************************************************
c **************************************************************************************************
c Declaration de variable


      Common/erreur_APF/err_absAPF,err_relAPF  
      Common/lim/nlimit 
      Common/comptage/kount  
      Common/dim/a,b,z
      Common/mode/m,n
      Common/onde/k0,Y0
      Common/plasma/X0,D0
      Common/pi/pi


      Integer *4 m,n,ifail
      Integer *4 nlimit,ni,nbre
      Integer *4 knout
      Real *8 pi,k0,Y0,f,cl,me,Eps0,qe
      Real *8 a,b(2),z(2)
      Real *8 n0,nc,Dn,D0,X0
      Real *8 err_absAPF,err_relAPF
      Complex *16 K




c **************************************************************************************************
c **************************************************************************************************
c constantes universelles

      pi=4.*ATAN(1.)
      cl=3.E8
      me=9.1091E-31
      Eps0=8.854E-12
      qe=1.6021E-19
      Y0=1/(120.*pi)
 
 
c **************************************************************************************************
c **************************************************************************************************      
c variables du programme


      f=3.7E9
      n0=1.E12*1.E6
      Dn=1.E12*1.E8
      
      a=70.E-3
      b(1)=8.E-3
      b(2)=8.E-3
      z(1)=0.
      z(2)=300.E-3
      
      m=2
      n=2




c ***********
c NAG D01APF

c     on fixe LW ds le prog
      
      err_absAPF=0.
      err_relAPF=1.E-4
      
      





c **************************************************************************************************
c **************************************************************************************************
c programme

      k0=2.*pi*f/cl
 
      nc=(2.*pi*f)**2*me*Eps0/qe**2
      X0=n0/nc

      D0=k0*n0/Dn  
      
      
      kount=0
      Call Int01APF(K,ifail,nbre)

      Write(*,*) 'intégrale avec APF'
      If (ifail.EQ.0) Then
          Write(*,*) 'convergence ok'
      Else
          Write(*,*) 'erreur convergence'
      Endif
      Write(*,*) 'après ',kount,'évaluations'
      Write(*,*) 'après ',nbre,'subdivisions'	 
      Write (*,*) 'K=',K



         
      End


c **************************************************************************************************
c **************************************************************************************************
c expression de la fonction


c ******************************

      Complex *16 function f_nz(nz)

      Common/dim/a,b,z
      Common/mode/m,n
      Common/onde/k0,Y0
      Common/plasma/X0,D0
      Common/pi/pi
      Common/comptage/kount 

      Integer *4 kount
      Integer *4 m,n,v,ifail
      Real *8 Y0,k0,X0,D0,nz,a,z(2),b(2),pi
      Real *8 coeff_mode_m,coeff_mode_n
      Real *8 PClimit
      Real *8 neta,S17AGF,S17AHF,S17AJF,S17AKF
      Complex *16 Ai,Aid,Bi,Bid,ys_mod

      kount=kount+1
      
      PClimit=1.e-300
 
      neta=(DABS(nz**2.-1.))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)

     
      If (DABS(nz).LT.1) Then
      
         ifail=1
         Aid=(1.,0.)*S17AJF(neta,ifail)

         ifail=1
	 Ai=(1.,0.)*S17AGF(neta,ifail)

         If (CDABS(Ai).LT.PClimit) Then
              f_nz=(0.,0.)
	      goto 20
	 Else      
              ys_mod=(0.,1.)*(Aid/Ai)*(X0/D0)**(1./3.)
	      goto 10
	 Endif
	 
      Elseif(DABS(nz).GT.1) Then
      
         ifail=1
         Aid=(1.,0.)*S17AJF(-neta,ifail)

         ifail=1
	 Ai=(1.,0.)*S17AGF(-neta,ifail)
		  
	 ifail=1
         Bid=(1.,0.)*S17AKF(-neta,ifail)
	 
         ifail=1
	 Bi=(1.,0.)*S17AHF(-neta,ifail) 
	 
         If (CDABS((Ai+(0.,-1.)*Bi)).LT.PClimit) Then
              f_nz=(0.,0.)
	      goto 20
	 Else    
              ys_mod=(0.,1.)*(Aid+(0.,-1.)*Bid)/(Ai+(0.,-1.)*Bi)*
     &           (X0/D0)**(1./3.)
	     goto 10
	 Endif
      
      Else
         f_nz=(0.,0.)
	 goto 20
            
  10  EndIf
  

      f_nz=Y0/(2.*pi*k0)*nz**2.*ys_mod*
     &(1.-(-1.)**n*CDEXP((0.,1.)*k0*nz*b(1)))*
     &1/(nz**2-(n*pi/(k0*b(1)))**2)*
     &(1.-(-1.)**m*CDEXP((0.,-1.)*k0*nz*b(2)))*
     &1/(nz**2-(m*pi/(k0*b(2)))**2)*
     &CDEXP((0.,1.)*k0*nz*(z(1)-z(2)))

    
      If (m.EQ.0) Then  
          coeff_mode_m=DSQRT(2./(a*b(1)))
      Else
          coeff_mode_m=2.*m/(b(1)*DSQRT((b(1)/a)+m**2*(a/b(1))))
      Endif
      
      If (n.EQ.0) Then  
          coeff_mode_n=DSQRT(2./(a*b(2)))
      Else
         coeff_mode_n=2.*n/(b(2)*DSQRT((b(2)/a)+n**2*(a/b(2))))
      Endif
        
      f_nz=f_nz*(a/2.)*coeff_mode_m*coeff_mode_n

         
               
  20  Return
      End
      
      

c ********************************************

      Complex *16 function f_t(t)
      
      Real *8 t,nz
      Complex *16 f_nz,f_neg,f_pos
      External f_nz
      
      If (t.EQ.0) Then
          f_t=(0.,0.)
	  goto 30
      Else  
          
c ***************
c ******* nz < 0  
   
          nz=(t-1.)/t
          f_neg=f_nz(nz)/(2*t)**(2./3.)	
	  
c ***************
c ******* nz > 0

	  nz=(1.-t)/t
          f_pos=f_nz(nz)/(2*t)**(2./3.)	  
	  	  
	  f_t=f_neg+f_pos	
	    
      Endif
      
  30  Return
      End
      


   
       
c ********************************************
c separation des parties reelle et imaginaire

      
      Real *8 function rf_t(t)
      Real *8 t
      Complex *16 f_t
      External f_t
      rf_t=DBLE(f_t(t))
      Return
      End
      
      Real *8 function if_t(t)   
      Real *8 t
      Complex *16 f_t
      External f_t
      if_t=DIMAG(f_t(t))
      Return
      End
      

      



 
 
c **************************************************************************************************
c **************************************************************************************************
c intégration


c ******************************************************
c subroutine d'integration de 0 a 1 routine NAG D01APF


      Subroutine Int01APF(K,ifail,nbre)


      Common/erreur_APF/err_absAPF,err_relAPF            
      Common/comptage/kount 

      Integer *4 kount           
      Integer *4 m,n
         
         
      Integer *4 LW,LIW,IW(50000)           
      Integer *4 ifail,ifail_tamp,nbre 
      Real *8 W(200000)
      Real *8 err_absAPF,err_relAPF,err_r
      Real *8 min,max,tamp_K
      Complex *16 K
      Real *8 rf_t,if_t
      External rf_t,if_t
      
      LW=200000
      LIW=LW/4
 
 
      nbre=0
      ifail=0
 
      
      min=0.
      max=1./2.




      ifail_tamp=-1
      Call D01APF(rf_t,min,max,0.,-2./3.,1,err_absAPF,err_relAPF,
     &tamp_K,err_r,W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=(1.,0.)*tamp_K

      
      ifail_tamp=-1
      Call D01APF(if_t,min,max,0.,-2./3.,1,err_absAPF,err_relAPF,
     &tamp_K,err_r,W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=K+(0.,1.)*tamp_K


     
      min=1./2.
      max=1.



      ifail_tamp=-1
      Call D01APF(rf_t,min,max,-2./3.,0.,1,err_absAPF,err_relAPF,
     &tamp_K,err_r,W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=K+(1.,0.)*tamp_K

      
      ifail_tamp=-1
      Call D01APF(if_t,min,max,-2./3.,0.,1,err_absAPF,err_relAPF,
     &tamp_K,err_r,W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=K+(0.,1.)*tamp_K
  
      
            

      Return
      End
      
      
