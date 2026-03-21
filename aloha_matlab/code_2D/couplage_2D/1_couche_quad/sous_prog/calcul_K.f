


      Subroutine Calcul_K(K,eyt_ny_nz,ezt_ny_nz,hyt_ny_nz,hzt_ny_nz)


c **************************************************************************************************
c **************************************************************************************************
c Declaration de variable

      implicit none


      Common/dim/a1,b1,y1,z1,a2,b2,y2,z2
      Common/mode/m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/integral/min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Common/admittance_surface/ys_yy_stok_sp,ys_yz_stok_sp,
     &ys_zy_stok_sp,ys_zz_stok_sp
      Common/tab_admittance/ tab_rempli
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
      Common/tamp/ind_tamp

      Integer *4 m1,n1,TE_TM_1,m2,n2,TE_TM_2,ind_tamp
      Integer *4 tab_rempli,ind_i,quad_key
      Real *8 pi,cl,Eps0,k0,Y0,f,un
      Real *8 a1,b1,y1,z1,a2,b2,y2,z2
      Real *8 S,D,dD,P,dP
      Real *8 min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz,ny_tampon
      Real *8 nz_quad,distance,lambda,error,omega_y,omega_z
      Complex *16 K
      Complex *16 eyt_ny_nz(30000),ezt_ny_nz(30000)
      Complex *16 hyt_ny_nz(30000),hzt_ny_nz(30000)
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)

      
c **************************************************************************************************
c **************************************************************************************************
c constantes universelles

      un=1.
      pi=4.*DATAN(un)
      cl=3.E8
      Eps0=8.854E-12
      Y0=1/(120.*pi)
  

c **************************************************************************************************
c **************************************************************************************************
c programme

      ind_tamp=0
 
      Call integral_discr(K) 

c      write(*,*) ind_tamp

      Call spect_discr(eyt_ny_nz,ezt_ny_nz,hyt_ny_nz,hzt_ny_nz)


      Return
      End


c **************************************************************************************************
c **************************************************************************************************
c expression de la fonction


c ******************************

      
      subroutine chps_spect(ny,nz,f_ny_nz,ey_ny_nz,
     &                      ez_ny_nz,hy_ny_nz,hz_ny_nz,int_sp,num)      

      implicit none

      Common/dim/a1,b1,y1,z1,a2,b2,y2,z2
      Common/mode/m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/admittance_surface/ys_yy_stok_sp,ys_yz_stok_sp,
     &ys_zy_stok_sp,ys_zz_stok_sp
      Common/tab_admittance/ tab_rempli
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
      Common/tamp/ind_tamp
           
      Integer *4 m1,n1,TE_TM_1,m2,n2,TE_TM_2,ind_tamp
      Integer *4 ifail,n_Ai
      Integer *4 tab_rempli
      Integer *4 int_sp,num,case_plasma,quad_key
      Real *8 k0,Y0
      Real *8 S,D,dD,P,dP
      Real *8 ny,nz,pi,a1,b1,y1,z1,a2,b2,y2,z2
      Real *8 coeff_mode_1y,coeff_mode_1z
      Real *8 coeff_mode_2y,coeff_mode_2z
      Real *8 PC_limit,pertes,un,deux
      Real *8 sign_y,sign_z
      Real *8 alpha_r,a_r,neta_r
      Real *8 U,dU,phi,Ai,dAi
      Real *8 nz_quad,error,omega_y,omega_z
      Complex *16 TF_cos_m1_cj,TF_sin_n1_cj,TF_sin_m1_cj,TF_cos_n1_cj
      Complex *16 TF_cos_m2,TF_sin_n2,TF_sin_m2,TF_cos_n2
      Complex *16 Hy1_cj,Hz1_cj,Ey2,Ez2
      Complex *16 f_ny_nz
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz  
      Complex *16 ys_yy,ys_yz,ys_zy,ys_zz          
      Complex *16 neta,alpha,Aiz,dAiz
      Complex *16 E_conj,dE_conj
      Complex *16 g_F,g_S,nxF,nxS
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)
 


      un=1.
      deux=2.
      pi=4.*DATAN(un)
      PC_limit=1.E-300

 
c coefficients de normalisation


      If (TE_TM_1.EQ.1) Then
      
          coeff_mode_1y=(n1/b1)/DSQRT(m1**2*b1/a1+n1**2*a1/b1)
	  coeff_mode_1z=-(m1/a1)/DSQRT(m1**2*b1/a1+n1**2*a1/b1)
          If (m1.GT.0) Then
	     coeff_mode_1y=coeff_mode_1y*DSQRT(deux)
	     coeff_mode_1z=coeff_mode_1z*DSQRT(deux)
	  Endif
	  If (n1.GT.0) Then
	     coeff_mode_1y=coeff_mode_1y*DSQRT(deux)
	     coeff_mode_1z=coeff_mode_1z*DSQRT(deux)
	  Endif
	  
      Else 
            
	  coeff_mode_1y=-2.*(m1/a1)/DSQRT(m1**2*b1/a1+n1**2*a1/b1)
          coeff_mode_1z=-2.*(n1/b1)/DSQRT(m1**2*b1/a1+n1**2*a1/b1)
	  
      Endif
      
      
      If (TE_TM_2.EQ.1) Then
      
          coeff_mode_2y=(n2/b2)/DSQRT(m2**2*b2/a2+n2**2*a2/b2)
	  coeff_mode_2z=-(m2/a2)/DSQRT(m2**2*b2/a2+n2**2*a2/b2)
          If (m2.GT.0) Then
	     coeff_mode_2y=coeff_mode_2y*DSQRT(deux)
	     coeff_mode_2z=coeff_mode_2z*DSQRT(deux)
	  Endif
	  If (n2.GT.0) Then
	     coeff_mode_2y=coeff_mode_2y*DSQRT(deux)
	     coeff_mode_2z=coeff_mode_2z*DSQRT(deux)
	  Endif
	  
      Else 
            
	  coeff_mode_2y=-2.*(m2/a2)/DSQRT(m2**2*b2/a2+n2**2*a2/b2)
          coeff_mode_2z=-2.*(n2/b2)/DSQRT(m2**2*b2/a2+n2**2*a2/b2)
	  
      Endif
 
 
c TF des cosinus et sinus

      If (((m1*m2+DABS(ny)).EQ.0).OR.((n1*n2+DABS(nz)).EQ.0)) Then
      
         f_ny_nz=(0.,0.)
         ey_ny_nz=(0.,0.)
         ez_ny_nz=(0.,0.)
         hy_ny_nz=(0.,0.)
         hz_ny_nz=(0.,0.)
	 goto 8
	 
      Endif


      TF_cos_m1_cj=((0.,1.)*ny/k0)*
     &(1.-(-1.)**m1*CDEXP((0.,1.)*k0*ny*a1))*
     &1/(ny**2-(m1*pi/(k0*a1))**2)*
     &CDEXP((0.,1.)*k0*ny*y1)
      TF_sin_n1_cj=(-n1*pi/(b1*k0**2))*
     &(1.-(-1.)**n1*CDEXP((0.,1.)*k0*nz*b1))*
     &1/(nz**2-(n1*pi/(k0*b1))**2)*
     &CDEXP((0.,1.)*k0*nz*z1)
      TF_sin_m1_cj=(-m1*pi/(a1*k0**2))*
     &(1.-(-1.)**m1*CDEXP((0.,1.)*k0*ny*a1))*
     &1/(ny**2-(m1*pi/(k0*a1))**2)*
     &CDEXP((0.,1.)*k0*ny*y1)
      TF_cos_n1_cj=((0.,1.)*nz/k0)*
     &(1.-(-1.)**n1*CDEXP((0.,1.)*k0*nz*b1))*
     &1/(nz**2-(n1*pi/(k0*b1))**2)*
     &CDEXP((0.,1.)*k0*nz*z1)
      TF_cos_m2=(-(0.,1.)*ny/k0)*
     &(1.-(-1.)**m2*CDEXP(-(0.,1.)*k0*ny*a2))*
     &1/(ny**2-(m2*pi/(k0*a2))**2)*
     &CDEXP(-(0.,1.)*k0*ny*y2)
      TF_sin_n2=(-n2*pi/(b2*k0**2))*
     &(1.-(-1.)**n2*CDEXP(-(0.,1.)*k0*nz*b2))*
     &1/(nz**2-(n2*pi/(k0*b2))**2)*
     &CDEXP(-(0.,1.)*k0*nz*z2)
      TF_sin_m2=(-m2*pi/(a2*k0**2))*
     &(1.-(-1.)**m2*CDEXP(-(0.,1.)*k0*ny*a2))*
     &1/(ny**2-(m2*pi/(k0*a2))**2)*
     &CDEXP(-(0.,1.)*k0*ny*y2)
      TF_cos_n2=(-(0.,1.)*nz/k0)*
     &(1.-(-1.)**n2*CDEXP(-(0.,1.)*k0*nz*b2))*
     &1/(nz**2-(n2*pi/(k0*b2))**2)*
     &CDEXP(-(0.,1.)*k0*nz*z2)
     
c expression des modes
     
      Hy1_cj=-coeff_mode_1z*TF_sin_m1_cj*TF_cos_n1_cj
      Hz1_cj=coeff_mode_1y*TF_cos_m1_cj*TF_sin_n1_cj
      Ey2=coeff_mode_2y*TF_cos_m2*TF_sin_n2
      Ez2=coeff_mode_2z*TF_sin_m2*TF_cos_n2


      If ((tab_rempli.EQ.1).AND.(int_sp.EQ.2)) Then

            ys_yy=ys_yy_stok_sp(num)
            ys_yz=ys_yz_stok_sp(num)
            ys_zy= ys_zy_stok_sp(num)
            ys_zz=ys_zz_stok_sp(num)

            goto 2

      Endif
    
 
c expression de l'admittance       
 

      case_plasma=4


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      If (case_plasma.EQ.1) Then

c vide


      If (DABS(1-ny**2-nz**2).LT.1.E-4) Then  

          ys_yy=(0.,0.)
          ys_yz=(0.,0.)
          ys_zy=(0.,0.)
          ys_zz=(0.,0.)
          goto 2

      Endif
    
      ys_yy=-(0.,1.)*ny*nz/CDSQRT(ny**2+nz**2-1.+1.E-6*(0,1.))
      ys_yz=(0.,1.)*(ny**2-1)/CDSQRT(ny**2+nz**2-1.+1.E-6*(0,1.))
      ys_zy=-(0.,1.)*(nz**2-1)/CDSQRT(ny**2+nz**2-1.+1.E-6*(0,1.))
      ys_zz=-ys_yy
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ElseIf (case_plasma.EQ.2) Then

c plasma 2D step 

      If ((S-ny**2-nz**2-D**2/(S-nz**2)).LT.0.) Then
         nxF=(0.,1.)*DSQRT(D**2/(S-nz**2)+ny**2+nz**2-S)
      Else    
         sign_z=DABS((S-nz**2)/(S-ny**2-nz**2))*
     &	                     (S-ny**2-nz**2)/(S-nz**2)
         nxF=-sign_z*DSQRT(S-ny**2-nz**2-D**2/(S-nz**2))     
      Endif

      If (((P/S)*(S-nz**2)-ny**2).LT.0.) Then
         nxS=(0.,1.)*DSQRT(DABS((P/S)*(S-nz**2)-ny**2))
      Else
         sign_z=DABS(S/(S-nz**2))*(S-nz**2)/S
         nxS=-sign_z*DSQRT((P/S)*(S-nz**2)-ny**2)    
      Endif
        
      ys_yy=nz*(nxF*ny+D*(0.,1.))/(S-ny**2-nz**2)
      ys_yz=(nxS*S*(S-nz**2-ny**2)+(nz**2)*ny*(nxF*ny+D*(0.,1.)))/
     &      ((S-nz**2-ny**2)*(S-nz**2))
      ys_zy=-(nxF*(S-nz**2)+ny*D*(0.,1.))/(S-nz**2-ny**2)
      ys_zz=-(ny*nz/(S-nz**2))*(nxF*(S-nz**2)+ny*D*(0.,1.))/
     &      (S-nz**2-ny**2)    
         


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ElseIf (case_plasma.EQ.3) Then

 
c plasma 1D lin�aire

      If (DABS(S-nz**2.).LT.1.E-2) Then
     
          ys_yy=(0.,0.)
          ys_yz=(0.,0.)
          ys_zy=(0.,0.)
          ys_zz=(0.,0.)
          goto 2
       
      Endif
    
      alpha_r=(DABS((S-nz**2.)*dP/S))**(1./3.)
      neta_r=alpha_r*((P/dP)-(S/dP)*ny**2/(S-nz**2))

      If (((S-nz**2.)*dP/S).GT.0.) Then
	      
          Call eval_airy_1(neta_r,Aiz,dAiz)
          g_S=alpha_r*dAiz/Aiz
      Else  

          Call eval_airy_2(neta_r,Ai,dAi)      
          g_S=alpha_r*dAi/Ai
      
      Endif


 
      ys_yy= (0.,0.)
      ys_yz=-(0.,1.)*S*g_S/(S-nz**2.)             
      ys_zy= (0.,0.)
      ys_zz= (0.,0.)




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ElseIf (case_plasma.EQ.4) Then


c plasma 2D lin�aire

      If ((DABS(S-ny**2.-nz**2.).LT.1.E-2).OR.
     &(DABS(S-nz**2.).LT.1.E-2)) Then
     
          ys_yy=(0.,0.)
          ys_yz=(0.,0.)
          ys_zy=(0.,0.)
          ys_zz=(0.,0.)
          goto 2
       
      Endif

      alpha_r=(DABS(nz**2.-S)/((2*dD)**2.))**(1./4.)	     
      a_r=(alpha_r**2)*(ny**2+nz**2-S+dD*ny/(S-nz**2))
      neta_r=(D/dD)/alpha_r 

      If ((nz**2-S).GT.0.) Then
	      
          Call eval_whittaker_E(a_r,neta_r,E_conj,dE_conj)
          g_F=(1./alpha_r)*dE_conj/E_conj
      Else  

          Call eval_parabolic_U(a_r,neta_r,U,dU)      
          g_F=((1.,0.)/alpha_r)*dU/U
      
      Endif
  

      alpha_r=(DABS((S-nz**2.)*dP/S))**(1./3.)
      neta_r=alpha_r*((P/dP)-(S/dP)*ny**2/(S-nz**2))

      If (((S-nz**2.)*dP/S).GT.0.) Then
	      
          Call eval_airy_1(neta_r,Aiz,dAiz)
          g_S=alpha_r*dAiz/Aiz
      Else  

          Call eval_airy_2(neta_r,Ai,dAi)      
          g_S=(1.,0.)*alpha_r*dAi/Ai
      
      Endif


      ys_yy= -(0.,1.)*ny*nz*g_F/(S-ny**2.-nz**2.)
     &       +(0.,1.)*nz*D/(S-ny**2.-nz**2.)

      ys_yz=-(0.,1.)*S*g_S/(S-nz**2.)
     &      -(0.,1.)*ny**2.*nz**2.*g_F/
     &       ((S-nz**2.)*(S-ny**2.-nz**2.))
     &      +(0.,1.)*ny*nz**2.*D/
     &       ((S-nz**2.)*(S-ny**2.-nz**2.))
               
      ys_zy=+(0.,1.)*(S-nz**2.)*g_F/(S-ny**2-nz**2.)
     &      -(0.,1.)*ny*D/(S-ny**2.-nz**2.) 
     &                         
 
      ys_zz=+(0.,1.)*ny*nz*g_F/(S-ny**2.-nz**2.)  
     &      -(0.,1.)*ny**2.*nz*D/
     &       ((S-nz**2.)*(S-ny**2.-nz**2.))


 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      EndIf 
    




      If (int_sp.EQ.2) Then
c         write(*,*) 'num=',num, size(ys_yy_stok_sp)
         ys_yy_stok_sp(num)=ys_yy
         ys_yz_stok_sp(num)=ys_yz
         ys_zy_stok_sp(num)=ys_zy
         ys_zz_stok_sp(num)=ys_zz

      Endif
  

c expression de l'int�grant  

 2    f_ny_nz=Y0*k0**2/(4.*pi**2)*
     &(Hy1_cj*(ys_yy*Ey2+ys_yz*Ez2)+
     &Hz1_cj*(ys_zy*Ey2+ys_zz*Ez2))


      If (CDABS(f_ny_nz)*0.EQ.0.) Then

          ey_ny_nz=Ey2
          ez_ny_nz=Ez2
          hy_ny_nz=Y0*ys_yy*Ey2+Y0*ys_yz*Ez2
          hz_ny_nz=Y0*ys_zy*Ey2+Y0*ys_zz*Ez2

      Else
     
          f_ny_nz=(0.,0.)
          ey_ny_nz=(0.,0.)
          ez_ny_nz=(0.,0.)
          hy_ny_nz=(0.,0.)
          hz_ny_nz=(0.,0.)

       
      Endif
      
      ind_tamp=ind_tamp+1
  
          
 8    Return
      End
      
      




 
 
c **************************************************************************************************
c **************************************************************************************************
c int�gration


      Real *8  function real_f(ny)

      implicit none

      Common/quad/nz_quad,error,omega_y,omega_z,quad_key

      Integer *4 int_sp,num,quad_key    
      Real *8 ny,nz,nz_quad,error,omega_y,omega_z
      Complex *16 f_ny_nz,f_ny_nz_conj
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz 

      int_sp=1
      num=1
      nz=nz_quad
      Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                hy_ny_nz,hz_ny_nz,int_sp,num)      
      f_ny_nz_conj=CONJG(f_ny_nz)
      real_f=(f_ny_nz+f_ny_nz_conj)/2.

      Return
      End


      Real *8 function imag_f(ny)

      implicit none

      Common/quad/nz_quad,error,omega_y,omega_z,quad_key

      Integer *4 int_sp,num,quad_key    
      Real *8 ny,nz,nz_quad,error,omega_y,omega_z
      Complex *16 f_ny_nz,f_ny_nz_conj
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz 

      int_sp=1
      num=1
      nz=nz_quad
      Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                hy_ny_nz,hz_ny_nz,int_sp,num)     
      f_ny_nz_conj=CONJG(f_ny_nz)
      imag_f=-(0.,1.)*(f_ny_nz-f_ny_nz_conj)/2.

      Return
      End




      Real *8 function real_F_nz(nz)

      implicit none

      Common/plasma/S,D,dD,P,dP    
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
 
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_ny,b_ny,nz,nz_quad,real_f,epsabs,epsrel,points(2)
      Real *8 abserr,work(2000),S,D,dD,P,dP,tamp,epsabsp,epsrelp 
      Real *8 borne_max,demi_intervalle,error,omega_y,omega_z

      External real_f

 
      nz_quad=nz

      npts2=2
      limit=400
      leniw=802
      lenw=2000

      borne_max=3.0
      demi_intervalle=0.2


      epsabs=error
      epsrel=error
      key=quad_key
      epsabsp=error
      epsrelp=error


      If (DSQRT(DABS(S-nz**2.)).LE.borne_max) Then

         If (DSQRT(DABS(S-nz**2.)).GE.(borne_max-demi_intervalle)) Then

            a_ny=-borne_max
            b_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=-DSQRT(DABS(S-nz**2.))
            Call dqagp(real_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            real_F_nz=tamp


            a_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(real_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

  
            a_ny=+DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=borne_max
            points(1)=DSQRT(DABS(S-nz**2.))
            Call dqagp(real_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

         Else

            a_ny=-borne_max
            b_ny=-DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(real_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            real_F_nz=tamp

            a_ny=-DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=-DSQRT(DABS(S-nz**2.))
            Call dqagp(real_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

            a_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(real_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

  
            a_ny=+DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=+DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=DSQRT(DABS(S-nz**2.))
            Call dqagp(real_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

            a_ny=DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=borne_max
            Call dqag(real_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            real_F_nz=real_F_nz+tamp

         Endif

      Else

         a_ny=-borne_max
         b_ny=borne_max
         Call dqag(real_f,a_ny,b_ny,epsabs,epsrel,key,real_F_nz,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)

      Endif


      Return
      End


      Real *8 function imag_F_nz(nz)

      implicit none
   
      Common/plasma/S,D,dD,P,dP  
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
 
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_ny,b_ny,nz,nz_quad,imag_f,epsabs,epsrel,points(2)
      Real *8 abserr,work(2000),S,D,dD,P,dP,tamp,epsabsp,epsrelp
      Real *8 borne_max,demi_intervalle,error,omega_y,omega_z

      External imag_f


      nz_quad=nz

      npts2=2
      limit=400
      leniw=802
      lenw=2000

      borne_max=3.0 
      demi_intervalle=0.2

      epsabs=error
      epsrel=error
      key=quad_key
      epsabsp=error
      epsrelp=error


      If (DSQRT(DABS(S-nz**2.)).LE.borne_max) Then

         If (DSQRT(DABS(S-nz**2.)).GE.(borne_max-demi_intervalle)) Then

            a_ny=-borne_max
            b_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=-DSQRT(DABS(S-nz**2.))
            Call dqagp(imag_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            imag_F_nz=tamp


            a_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(imag_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

  
            a_ny=+DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=borne_max
            points(1)=DSQRT(DABS(S-nz**2.))
            Call dqagp(imag_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

         Else

            a_ny=-borne_max
            b_ny=-DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(imag_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            imag_F_nz=tamp

            a_ny=-DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=-DSQRT(DABS(S-nz**2.))
            Call dqagp(imag_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

            a_ny=-DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=DSQRT(DABS(S-nz**2.))-demi_intervalle
            Call dqag(imag_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

  
            a_ny=+DSQRT(DABS(S-nz**2.))-demi_intervalle
            b_ny=+DSQRT(DABS(S-nz**2.))+demi_intervalle
            points(1)=DSQRT(DABS(S-nz**2.))
            Call dqagp(imag_f,a_ny,b_ny,npts2,points,epsabsp,
     &         epsrelp,tamp,abserr,neval,ier,leniw,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

            a_ny=DSQRT(DABS(S-nz**2.))+demi_intervalle
            b_ny=borne_max
            Call dqag(imag_f,a_ny,b_ny,epsabs,epsrel,key,tamp,abserr,
     &                neval,ier,limit,lenw,last,iwork,work)
            imag_F_nz=imag_F_nz+tamp

         Endif

      Else

         a_ny=-borne_max
         b_ny=borne_max
         Call dqag(imag_f,a_ny,b_ny,epsabs,epsrel,key,imag_F_nz,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)

      Endif

      Return
      End



      Subroutine integral_real_f(Kr)

      implicit none

      Common/plasma/S,D,dD,P,dP 
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
    
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      integer *4 quad_key
      Real *8 a_nz,b_nz,epsabs,epsrel,epsabsp,epsrelp,real_F_nz
      Real *8 Kr,points(2),abserr,work(2000),S,D,dD,P,dP,K_tamp 
      Real *8 borne_max,demi_intervalle,error,nz_quad,omega_y,omega_z

      External real_F_nz


      npts2=2
      limit=400
      leniw=802
      lenw=2000

      borne_max=24.
      demi_intervalle=0.2

      epsabs=error
      epsrel=error
      key=quad_key
      epsabsp=error
      epsrelp=error

      a_nz=-borne_max
      b_nz=-DSQRT(DABS(S))-demi_intervalle
      Call dqag(real_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Kr=K_tamp


      a_nz=-DSQRT(DABS(S))-demi_intervalle
      b_nz=-DSQRT(DABS(S))+demi_intervalle
      points(1)=-DSQRT(DABS(S))
      Call dqagp(real_F_nz,a_nz,b_nz,npts2,points,epsabsp,epsrelp,K_tamp
     &           ,abserr,neval,ier,leniw,lenw,last,iwork,work)
      Kr=Kr+K_tamp


      a_nz=-DSQRT(DABS(S))+demi_intervalle
      b_nz=DSQRT(DABS(S))-demi_intervalle
      Call dqag(real_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Kr=Kr+K_tamp

  
      a_nz=+DSQRT(DABS(S))-demi_intervalle
      b_nz=+DSQRT(DABS(S))+demi_intervalle
      points(1)=DSQRT(DABS(S))
      Call dqagp(real_F_nz,a_nz,b_nz,npts2,points,epsabsp,epsrelp,K_tamp
     &           ,abserr,neval,ier,leniw,lenw,last,iwork,work)
      Kr=Kr+K_tamp


      a_nz=+DSQRT(DABS(S))+demi_intervalle
      b_nz=borne_max
      Call dqag(real_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Kr=Kr+K_tamp

      Return
      End



      Subroutine integral_imag_f(Ki)

      implicit none

      Common/plasma/S,D,dD,P,dP 
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key
    
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_nz,b_nz,epsabs,epsrel,epsabsp,epsrelp,imag_F_nz
      Real *8 Ki,points(2),abserr,work(2000),S,D,dD,P,dP,K_tamp 
      Real *8 borne_max,demi_intervalle,error,nz_quad,omega_y,omega_z

      External imag_F_nz


      npts2=2
      limit=400
      leniw=802
      lenw=2000

      borne_max=24.
      demi_intervalle=0.2

      epsabs=error
      epsrel=error
      key=quad_key
      epsabsp=error
      epsrelp=error


      a_nz=-borne_max
      b_nz=-DSQRT(DABS(S))-demi_intervalle
      Call dqag(imag_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Ki=K_tamp


      a_nz=-DSQRT(DABS(S))-demi_intervalle
      b_nz=-DSQRT(DABS(S))+demi_intervalle
      points(1)=-DSQRT(DABS(S))
      Call dqagp(imag_F_nz,a_nz,b_nz,npts2,points,epsabsp,epsrelp,K_tamp
     &           ,abserr,neval,ier,leniw,lenw,last,iwork,work)
      Ki=Ki+K_tamp


      a_nz=-DSQRT(DABS(S))+demi_intervalle
      b_nz=DSQRT(DABS(S))-demi_intervalle
      Call dqag(imag_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Ki=Ki+K_tamp

  
      a_nz=+DSQRT(DABS(S))-demi_intervalle
      b_nz=+DSQRT(DABS(S))+demi_intervalle
      points(1)=DSQRT(DABS(S))
      Call dqagp(imag_F_nz,a_nz,b_nz,npts2,points,epsabsp,epsrelp,K_tamp
     &           ,abserr,neval,ier,leniw,lenw,last,iwork,work)
      Ki=Ki+K_tamp


      a_nz=+DSQRT(DABS(S))+demi_intervalle
      b_nz=borne_max
      Call dqag(imag_F_nz,a_nz,b_nz,epsabs,epsrel,key,K_tamp,abserr,
     &        neval,ier,limit,lenw,last,iwork,work)
      Ki=Ki+K_tamp


      Return
      End






      Subroutine integral_discr(K) 

      implicit none
 
      Common/dim/a1,b1,y1,z1,a2,b2,y2,z2
      Common/onde/k0,Y0
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key

      Integer *4 quad_key
      Real *8 a1,b1,y1,z1,a2,b2,y2,z2,k0,Y0
      Real *8 nz_quad,error,omega_y,omega_z
      Real *8 Kr,Ki
      Complex *16 K


      error=1.E-3
      quad_key=6
      Call integral_real_f(Kr)
      Call integral_imag_f(Ki)
      K = (1.,0.)*Kr + (0.,1.)*Ki
 

      Return 

      End




      
       

 
 
c **************************************************************************************************
c **************************************************************************************************
c donnees spectrales



      Subroutine spect_discr(eyt_ny_nz,ezt_ny_nz,hyt_ny_nz,hzt_ny_nz)

      implicit none

      Common/integral/min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Common/quad/nz_quad,error,omega_y,omega_z,quad_key

      Integer *4 nbre,n_1,n_2,incr,int_sp,quad_key
      Real *8 min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Real *8 pas_ny,pas_nz,ny,nz,nz_quad,error,omega_y,omega_z
      Complex *16 f_ny_nz
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz 
      Complex *16 eyt_ny_nz(30000),ezt_ny_nz(30000)
      Complex *16 hyt_ny_nz(30000),hzt_ny_nz(30000)

      int_sp=2

      pas_ny=(max_ny-min_ny)/nbre_ny
      pas_nz=(max_nz-min_nz)/nbre_nz

      incr=1
     
      Do 40 n_1=0,nbre_nz
      
         nz=min_nz+n_1*pas_nz
	 
	 Do 50 n_2=0,nbre_ny
	 
            ny=min_ny+n_2*pas_ny
            Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,
     &                      ez_ny_nz,hy_ny_nz,hz_ny_nz,int_sp,incr)  

            eyt_ny_nz(incr)=ey_ny_nz
            ezt_ny_nz(incr)=ez_ny_nz
            hyt_ny_nz(incr)=hy_ny_nz
            hzt_ny_nz(incr)=hz_ny_nz

            incr=incr+1
	    
   
 50      Continue  
 
 40   Continue  
        
      Return
      End


       

c *********************************************************************************
c *********************************************************************************
c fonctions sp�ciales



c ********************************************
c fonction parabolique (parfois � des facteurs pr�s ; voir abramovitz)
c le rapport d�riv�e/fonction est correct

      Subroutine eval_parabolic_U(a,x,U,dU)

      implicit none

      Real *8 a,x,U,dU
      Real *8 v_a,p,vr,vi,dvr,dvi,v1,dv1,x_moins
      Real *8 un,pi

      un=1.
      pi=4.*DATAN(un)

      If ((a.GT.(-70.)).AND.(a.LT.70.)) Then

	 v_a=-a-1./2.
         Call PARABOLIC_D(v_a,x,U,dU)      
	  
      Elseif (a.LT.(-70.)) Then

c valable pour |a|>>x� et a<0

         p=DSQRT(-a)

	 vr=x**2./(16.*p**2.)+x**4./(2.**7.*p**4.)-
     &      9.*x**2./(2.**8.*p**6.)+ x**6./(3*2.**8.*p**6)

         vi=-x**3./(24.*p)+x/(16.*p**3.)+
     &     x**5./(5.*2.**7.*p**3.)+x**3./(48.*p**5.)-
     &     x**7/(7*2.**10.*p**5.)
	  
	 dvr=x/(8.*p**2.)+x**3./(2.**5.*p**4.)-
     &      9.*x/(2.**7.*p**6.)+ x**5./(2.**7.*p**6)

         dvi=-x**2./(8.*p)+1./(16.*p**3.)+
     &      x**4./(2.**7.*p**3.)+x**2./(16.*p**5.)-
     &      x**6/(2.**10.*p**5.)


         U=DCOS(p*x+vi+pi/4+pi*a/2)
         dU=(dvr*DCOS(p*x+vi+pi/4+pi*a/2)-
     &                       (p+dvi)*DSIN(p*x+vi+pi/4+pi*a/2))

	  
      Else

c valable pour |a|>>x� et a>0

         p=DSQRT(a)

         v1=-x**3./(24.*p)-x**2./(16.*p**2)-x/(16.*p**3.)+
     &     x**5./(5*2.**7.*p**3.)+x**4./(2.**7.*p**4.)+
     &     x**3./(48.*p**5.)-x**7./(7*2.**10.*p**5.) 
            
         dv1=-x**2./(8.*p)-x/(8.*p**2)-1./(16.*p**3.)+
     &      x**4./(2.**7.*p**3.)+x**3./(2.**5.*p**4.)+
     &      x**2./(16.*p**5.)-x**6./(2.**10.*p**5.) 

         U=1.
         dU=(-p+dv1)


      Endif
   


      Return
      End








c ********************************************
c fonction Whittaker (parfois � des facteurs pr�s ; voir abramovitz)
c le rapport d�riv�e/fonction est correct

      Subroutine eval_whittaker_E(a,x,E_conj,dE_conj)

      implicit none

      Real *8 a,x,W1,dW1,W2,dW2,ke
      Real *8 v_a,p,vr,vi,dvr,dvi,v1,dv1,v2,dv2
      Real *8 un,pi
      Complex *16 E,dE,E_conj,dE_conj

      un=1.
      pi=4.*DATAN(un)


      If (a.LT.(-10.)) Then

c     valable pour |a|>>x� et a<0

         p=DSQRT(-a)
	  
	 vr=-x**2./(16.*p**2.)+x**4./(2.**7.*p**4.)-
     &      9.*x**2./(2.**8.*p**6.)- x**6./(3*2.**8.*p**6)

         vi=x**3./(24.*p)-x/(16.*p**3.)-
     &     x**5./(5.*2.**7.*p**3.)+x**3./(48.*p**5.)+
     &     x**7/(7*2.**10.*p**5.)

	 dvr=-x/(8.*p**2.)+x**3./(2.**5.*p**4.)-
     &      9.*x/(2.**7.*p**6.)- x**5./(2.**7.*p**6)

         dvi=x**2./(8.*p)-1./(16.*p**3.)-
     &      x**4./(2.**7.*p**3.)+x**2./(16.*p**5.)+
     &      x**6/(2.**10.*p**5.)

         W1=DCOS(p*x+vi+pi/4)
         dW1=dvr*DCOS(p*x+vi+pi/4)-(p+dvi)*DSIN(p*x+vi+pi/4)
         W2=DSIN(p*x+vi+pi/4)
         dW2=-(dvr*DSIN(p*x+vi+pi/4)+(p+dvi)*DCOS(p*x+vi+pi/4))

         ke=1.-DEXP(pi*a)        
         E_conj=(1.,0.)*W1-(0.,1.)*ke*W2
         dE_conj=(1.,0.)*dW1-(0.,1.)*ke*dW2 


      Elseif (((a.LT.60.).AND.(x.LT.1.5)).OR.
     &       ((a.LT.40.).AND.(x.GT.1.5).AND.(x.LT.2.)).OR.   
     &       ((a.LT.30.).AND.(x.GT.2.).AND.(x.LT.2.5)).OR. 
     &       ((a.LT.20.).AND.(x.GT.2.5).AND.(x.LT.3.)).OR. 
     &       ((a.LT.15.).AND.(x.GT.3.).AND.(x.LT.3.5)).OR. 
     &       ((a.LT.12.).AND.(x.GT.3.5).AND.(x.LT.4.)).OR.
     &       ((a.LT.11.).AND.(x.GT.4.).AND.(x.LT.4.5)).OR.       
     &       ((a.LT.10.5).AND.(x.GT.4.5).AND.(x.LT.5.)) ) Then

c valable a priori pour |a|<5 et |x|<5 => �tendu via des tests num�riques

          Call PARABOLIC_W(a,x,W1,dW1,W2,dW2) 

          If (a.GT.2.) Then

              ke=DEXP(-pi*a)/2.
	    
	  Elseif (a.LT.(-4.)) Then

              ke=1.-DEXP(pi*a)
	    
	  Else

              ke=DSQRT(1.+DEXP(2.*pi*a))-DEXP(pi*a)

          Endif
	    
          E_conj=(1.,0.)*W1-(0.,1.)*ke*W2
          dE_conj=(1.,0.)*dW1-(0.,1.)*ke*dW2   
	  
 
	  
      Else

c valable pour |a|>>x� et a>0

         p=DSQRT(a)

         v1=+x**3./(24.*p)+x**2./(16.*p**2)+x/(16.*p**3.)+
     &     x**5./(5*2.**7.*p**3.)+x**4./(2.**7.*p**4.)+
     &     x**3./(48.*p**5.)+x**7./(7*2.**10.*p**5.) 
            
         dv1=x**2./(8.*p)+x/(8.*p**2)+1./(16.*p**3.)+
     &      x**4./(2.**7.*p**3.)+x**3./(2.**5.*p**4.)+
     &      x**2./(16.*p**5.)+x**6./(2.**10.*p**5.) 

         v2=-x**3./(24.*p)+x**2./(16.*p**2)-x/(16.*p**3.)-
     &     x**5./(5*2.**7.*p**3.)+x**4./(2.**7.*p**4.)-
     &     x**3./(48.*p**5.)-x**7./(7*2.**10.*p**5.) 
            
         dv2=-x**2./(8.*p)+x/(8.*p**2)-1./(16.*p**3.)+
     &      x**4./(2.**7.*p**3.)+x**3./(2.**5.*p**4.)-
     &      x**2./(16.*p**5.)-x**6./(2.**10.*p**5.) 

         W1=DEXP(-p*x+v1)
         dW1=(-p+dv1)*W1
         W2=DEXP(p*x+v2)
         dW2=-(p+dv2)*W2

         E_conj=DEXP(v1)-DEXP((1-pi)*a)*DEXP(v2)
         dE_conj=(-p+v1)*DEXP(v1)-
     &                 (p+dv2)*DEXP((1-pi)*a)*DEXP(v2)
	    
	    



      Endif
   

      Return
      End





c ********************************************
c fonction airy (parfois � des facteurs pr�s ; voir abramovitz)
c le rapport d�riv�e/fonction est correct


      Subroutine eval_airy_1(x,Aiz,dAiz)

      implicit none

      Real *8 x,Ai,dAi,Bi,dBi,x_moins,c(6),d(6)
      Real *8 un,pi
      Complex *16 Aiz,dAiz,ksi

      un=1.
      pi=4.*DATAN(un)

      c(1)=1.
      c(2)=0.069444444444444
      c(3)=0.037133487654321
      c(4)=0.037993059127801
      c(5)=0.057649190412670
      c(6)=0.116099064025515


      d(1)=1.
      d(2)=-0.097222222222222
      d(3)=-0.043885030864198
      d(4)=-0.042462830789895
      d(5)=-0.062662163492032
      d(6)=-0.124105896027275


      If (DABS(x).LT.50.) Then

         x_moins=-x
         Call AIRYA(x_moins,Ai,Bi,dAi,dBi)
         Aiz=(1.,0.)*Ai-(0.,1.)*Bi
         dAiz=-((1.,0.)*dAi-(0.,1.)*dBi)
	  
      Else 

c valable pour |x|>>1

         If (x.GT.0.) Then

            ksi=(2./3.)*x**(3./2.)*CDEXP(-(0.,1.)*pi/2.)
            dAiz=-CDEXP(-(0.,1.)*pi/3)*x**(1./4.)*CDEXP(-(0.,1.)*pi/12.)
     &        *(d(1)-d(2)*ksi**(-1.)+d(3)*ksi**(-2.)-d(4)*ksi**(-3.)+
     &          d(5)*ksi**(-4.)-d(6)*ksi**(-5.))
            Aiz=x**(-1./4.)*CDEXP(+(0.,1.)*pi/12.)*
     &         (c(1)-c(2)*ksi**(-1.)+c(3)*ksi**(-2.)-c(4)*ksi**(-3.)+
     &          c(5)*ksi**(-4.)-c(6)*ksi**(-5.))

         Else 

            ksi=-(2./3.)*DABS(x)**(3./2.)
         
            dAiz=-CDEXP(-(0.,1.)*pi/3)*DABS(x)**(1./4.)*
     &            CDEXP((0.,1.)*pi/6.)
     &        *(d(1)-d(2)*ksi**(-1.)+d(3)*ksi**(-2.)-d(4)*ksi**(-3.)+
     &          d(5)*ksi**(-4.)-d(6)*ksi**(-5.))
            Aiz=DABS(x)**(-1./4.)*CDEXP(-(0.,1.)*pi/6.)*
     &         (c(1)-c(2)*ksi**(-1.)+c(3)*ksi**(-2.)-c(4)*ksi**(-3.)+
     &          c(5)*ksi**(-4.)-c(6)*ksi**(-5.))

         Endif
 
      Endif
   


      Return
      End



      Subroutine eval_airy_2(x,Ai,dAi)

      implicit none

      Real *8 x,Ai,dAi,Bi,dBi,x_moins,c(6),d(6)
      Real *8 un,pi,ksi
      Complex *16 Aiz,dAiz

      un=1.
      pi=4.*DATAN(un)

      c(1)=1.
      c(2)=0.069444444444444
      c(3)=0.037133487654321
      c(4)=0.037993059127801
      c(5)=0.057649190412670
      c(6)=0.116099064025515


      d(1)=1.
      d(2)=-0.097222222222222
      d(3)=-0.043885030864198
      d(4)=-0.042462830789895
      d(5)=-0.062662163492032
      d(6)=-0.124105896027275


      If (DABS(x).LT.50.) Then

         Call AIRYA(x,Ai,Bi,dAi,dBi)
	  
      Else 

c valable pour |x|>>1

         ksi=(2./3.)*DABS(x)**(3./2.)

         If (x.GT.0.) Then

            dAi=-x**(1./4.)*(d(1)-d(2)*ksi**(-1.)+d(3)*ksi**(-2.)-
     &           d(4)*ksi**(-3.)+d(5)*ksi**(-4.)-d(6)*ksi**(-5.))
            Ai=x**(-1./4.)*(c(1)-c(2)*ksi**(-1.)+c(3)*ksi**(-2.)-
     &           c(4)*ksi**(-3.)+c(5)*ksi**(-4.)-c(6)*ksi**(-5.))

         Else

            x_moins=-x

            dAi=-x_moins**(1./4.)*(DCOS(ksi+pi/4.)*(d(1)-
     &           d(3)*ksi**(-2.)+d(5)*ksi**(-4.))+DSIN(ksi+pi/4.)*
     &          (d(2)*ksi**(-1.)-d(4)*ksi**(-3.)+d(6)*ksi**(-5.)))
            Ai=x_moins**(-1./4.)*(DSIN(ksi+pi/4.)*(c(1)-
     &           c(3)*ksi**(-2.)+c(5)*ksi**(-4.))-DCOS(ksi+pi/4.)*
     &          (c(2)*ksi**(-1.)-c(4)*ksi**(-3.)+c(6)*ksi**(-5.)))

         Endif

 
      Endif
   


      Return
      End





c ********************************************
c calcule la phase d'un nombre complexe

      Subroutine phase(z,phi)

      implicit none

      Real *8 real_z,imag_z,phi,un,pi
      Complex *16 z,z_conj

      un=1.
      pi=4.*DATAN(un)
      
      z_conj=CONJG(z)
      real_z=(z+z_conj)/2.
      imag_z=-(0.,1.)*(z-z_conj)/2.

      if (real_z.GT.0.) Then
          phi=ATAN(imag_z/real_z)
      Else
          phi=ATAN(imag_z/real_z)+pi
      Endif


      Return
      End


























