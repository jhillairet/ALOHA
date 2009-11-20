
c **************************************************************************************************
c **************************************************************************************************
c Declaration de variable


      Common/erreur/errAHF,err_absAJF,err_relAJF
      Common/erreur_bis/err_absAMF,err_relAMF 
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
      Real *8 errAHF,err_absAJF,err_relAJF
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
      n0=0.6E12*1.E6
      Dn=0.6E12*1.E8
      
      a=70.E-3
      b(1)=8.E-3
      b(2)=8.E-3
      z(1)=0.
      z(2)=0.
      
      m=1
      n=1

      
c ***********
c NAG D01AMF

c     on fixe LW ds le prog
      
      err_absAMF=0.
      err_relAMF=1.E-6


      
      
c ***********
c NAG D01AJF

c     on fixe LW ds le prog
      
      err_absAJF=0.
      err_relAJF=1.E-6
      
      

c ***********
c NAG D01AHF
      
      nlimit=1.E5
      errAHF=1.E-9




c **************************************************************************************************
c **************************************************************************************************
c programme

      k0=2.*pi*f/cl
 
      nc=(2.*pi*f)**2*me*Eps0/qe**2
      X0=n0/nc

      D0=k0*n0/Dn  
   
      kount=0
      Call Int01AMF(K,ifail,nbre)
      Write(*,*) 'intégrale avec AMF'
      If (ifail.EQ.0) Then
          Write(*,*) 'convergence ok'
      Else
          Write(*,*) 'erreur convergence'
      Endif
      Write(*,*) 'après ',kount,'évaluations'
      Write(*,*) 'après ',nbre,'subdivisions'	 
      Write (*,*) 'K=',K
   

           
      kount=0
      Call Int01AJF(K,ifail,nbre)

      Write(*,*) 'intégrale  avec AJF'
      If (ifail.EQ.0) Then
          Write(*,*) 'convergence ok'
      Else
          Write(*,*) 'erreur convergence'
      Endif
      Write(*,*) 'après ',kount,'évaluations'
      Write(*,*) 'après ',nbre,'subdivisions'  	 
      Write (*,*) 'K=',K



      Call Int01AHF(K,ifail,ni)

      Write(*,*) 'intégrale avec AHF'
      If (ifail.EQ.0) Then
         Write(*,*) 'convergence ok'
      Else
          Write(*,*) 'erreur convergence'
      Endif
      Write(*,*) 'après ',ni,'évaluations'	 
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
      Real *8 PC_limit,pertes
      Complex *16 neta,tampon,Ai,Aid,Bi,Bid,ys,alpha
      
     
      External S17DGF,S17DHF

      kount=kount+1
      
      PC_limit=1.e-300
      
      pertes=1.0E-5
      alpha=(1.,0.)-(0.,1.)*pertes
      
  
      
      
      If (DABS(nz).EQ.1) Then
          f_nz=(0.,0.)
	  goto 10	  
      Else
    

          neta=((1.,0.)*nz**2-alpha)**(1./3.)*(D0/X0)**(2./3.)*
     &         (X0-1.)*CDEXP(-(0.,1.)*pi/3.)


          ifail=1
          Call S17DGF('f',neta,'u',Ai,v,ifail)

          ifail=1
          Call S17DGF('d',neta,'u',Aid,v,ifail)

          If (CDABS(Ai).LT.PC_limit) Then
             f_nz=(0.,0.)
	     goto 10	  
          Else      
            ys=-alpha*(Aid/Ai)*(X0/D0)**(1./3.)/
     &          ((1.,0.)*nz**2-alpha)**(2./3.)*CDEXP((0.,1.)*pi/6.)
  	             
          EndIf

      Endif
  

      f_nz=Y0/(2.*pi*k0)*nz**2.*ys*
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

         
               
  10  Return
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
          f_neg=f_nz(nz)/t**2.	
	  
c ***************
c ******* nz > 0

	  nz=(1.-t)/t
          f_pos=f_nz(nz)/t**2.	  
	  	  
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
      
      
      Real *8 function rf_nz(nz)
      Real *8 nz
      Complex *16 f_nz
      External f_nz
      If (abs(nz).GT.1.E4)Then
          rf_nz=0.
      Else
          rf_nz=DBLE(f_nz(nz))
      Endif
      Return
      End
      
      Real *8 function if_nz(nz)   
      Real *8 nz
      Complex *16 f_nz
      External f_nz
      If (abs(nz).GT.1.E4) Then
          if_nz=0.
      Else
          if_nz=DIMAG(f_nz(nz))
      Endif
      Return
      End


 
 
c **************************************************************************************************
c **************************************************************************************************
c intégration


******************************************************
c subroutine d'integration de 0 a 1 routine NAG D01AHF

      Subroutine Int01AHF(K,ifail,ni)


      Common/erreur/errAHF,err_absAJF,err_relAJF
      Common/lim/nlimit  

      Integer *4 ifail,ifail_tamp
      Integer *4 nlimit,ni,ni_tamp
      Real *8 errAHF,err_r
      Real *8 min,max
      Real *8 Krf,Kif
      Complex *16 K
      Real *8 D01AHF,rf_t,if_t
      External rf_t,if_t
          
      min=0.
      max=1.
      
      ifail=0
      ni=0
      
            
      ifail_tamp=-1
      Krf=D01AHF(min,max,errAHF,ni_tamp,err_r,rf_t,
     &nlimit,ifail_tamp)
      ifail=ifail+ifail_tamp
      ni=ni+ni_tamp
      
      ifail_tamp=-1
      Kif=D01AHF(min,max,errAHF,ni_tamp,err_r,if_t,
     &nlimit,ifail_tamp)
      ifail=ifail+ifail_tamp
      ni=ni+ni_tamp
      

      
      K=(1.,0.)*Krf+(0.,1.)*Kif

       
      Return
      End



******************************************************
c subroutine d'integration de 0 a 1 routine NAG D01AJF


      Subroutine Int01AJF(K,ifail,nbre)


      Common/erreur/errAHF,err_absAJF,err_relAJF            
      Common/comptage/kount 

      Integer *4 kount           
      Integer *4 m,n
         
         
      Integer *4 LW,LIW,IW(5000)           
      Integer *4 ifail,ifail_tamp,nbre 
      Real *8 W(20000)
      Real *8 err_absAJF,err_relAJF,err_r
      Real *8 min,max,tamp_K
      Complex *16 K
      Real *8 rf_t,if_t
      External rf_t,if_t

      
      LW=20000
      LIW=LW/4
      
      min=1./50.
      max=1.

      nbre=0
      ifail=0
      
      
      ifail_tamp=-1
      Call D01AJF(rf_t,min,max,err_absAJF,err_relAJF,tamp_K,err_r,
     &W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=(1.,0.)*tamp_K
 
      Write(*,*) IW(1)
      Write(*,*) ifail_tamp
      
      Open(10,File='data_W.dat',Status='unknown',Form='unformatted')  
      Write(10) (W(i),i=1,4*IW(1))
      Rewind 10
      Close(10)      

            
      ifail_tamp=-1
      Call D01AJF(if_t,min,max,err_absAJF,err_relAJF,tamp_K,err_r,
     &W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=K+(0.,1.)*tamp_K
      
      

             

      Return
      End
      





******************************************************
c subroutine d'integration de -infini a +infini routine NAG D01AMF

      Subroutine Int01AMF(K,ifail,nbre)


      Common/erreur_bis/err_absAMF,err_relAMF            
      Common/comptage/kount 

      Integer *4 kount           
      Integer *4 m,n
         
         
      Integer *4 LW,LIW,IW(5000)           
      Integer *4 ifail,ifail_tamp
      Integer *4 nbre 
      Real *8 W(20000)
      Real *8 err_absAMF,err_relAMF,err_r
      Real *8 min,max,tamp_K
      Complex *16 K
      Real *8 rf_nz,if_nz
      External rf_nz,if_nz
      
      LW=20000
      LIW=LW/4
      

      nbre=0
      ifail=0
            
      ifail_tamp=-1
      Call D01AMF(rf_nz,0.,2,err_absAMF,err_relAMF,tamp_K,err_r,
     &W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=(1.,0.)*tamp_K
 
 
                     
      ifail_tamp=-1
      Call D01AMF(if_nz,0.,2,err_absAMF,err_relAMF,tamp_K,err_r,
     &W,LW,IW,LIW,ifail_tamp)
      nbre=nbre+IW(1)
      ifail=ifail+ifail_tamp
      K=K+(0.,1.)*tamp_K
	    


      Return
      End
      

      

      
      
