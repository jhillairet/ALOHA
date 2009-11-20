


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
      Common/admittance_surface/gs_int,gf_int,ys_yy_stok_sp,
     &ys_yz_stok_sp,ys_zy_stok_sp,ys_zz_stok_sp
      Common/tab_admittance/tab_rempli
      Common/quad/nz_quad,error,pertes,quad_key


      Integer *4 m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Integer *4 tab_rempli,ind_i,quad_key
      Real *8 pi,cl,Eps0,k0,Y0,f,un
      Real *8 a1,b1,y1,z1,a2,b2,y2,z2
      Real *8 S,D,dD,P,dP
      Real *8 min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz,ny_tampon
      Real *8 nz_quad,error,pertes
      Complex *16 K
      Complex *16 eyt_ny_nz(30000),ezt_ny_nz(30000)
      Complex *16 hyt_ny_nz(30000),hzt_ny_nz(30000)
      Complex *16 gs_int(361501),gf_int(361501)
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
    

      quad_key=6
c MODIF JH 19-05-2009. Was: 1e-4
      error=1.E-4
c      error=1.E-5
      
      pertes=-10.E-2 
  
      If (tab_rempli.EQ.0) Then
         Call map_gS_gF
      Endif

      Call integral_discr(K)   

      Call spect_discr(eyt_ny_nz,ezt_ny_nz,hyt_ny_nz,hzt_ny_nz)


      Return
      End


c **************************************************************************************************
c **************************************************************************************************
c expression de la fonction


c ******************************

      
      subroutine chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                      hy_ny_nz,hz_ny_nz,int_sp,num)      

      implicit none

      Common/dim/a1,b1,y1,z1,a2,b2,y2,z2
      Common/mode/m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/admittance_surface/gs_int,gf_int,ys_yy_stok_sp,
     &ys_yz_stok_sp,ys_zy_stok_sp,ys_zz_stok_sp
      Common/tab_admittance/tab_rempli
      Common/quad/nz_quad,error,pertes,quad_key

           
      Integer *4 m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Integer *4 ifail,n_Ai
      Integer *4 tab_rempli,quad_key
      Integer *4 int_sp,num,case_plasma
      Real *8 k0,Y0
      Real *8 S,D,dD,P,dP
      Real *8 ny,nz,a1,b1,y1,z1,a2,b2,y2,z2
      Real *8 coeff_mode_1y,coeff_mode_1z
      Real *8 coeff_mode_2y,coeff_mode_2z
      Real *8 un,deux,pi
      Real *8 sign_y,sign_z
      Real *8 alpha_r,a_r,neta_r
      Real *8 U,dU,Ai,dAi
      Real *8 nz_quad,error,pertes
      Complex *16 TF_cos_m1_cj,TF_sin_n1_cj,TF_sin_m1_cj,TF_cos_n1_cj
      Complex *16 TF_cos_m2,TF_sin_n2,TF_sin_m2,TF_cos_n2
      Complex *16 Hy1_cj,Hz1_cj,Ey2,Ez2
      Complex *16 f_ny_nz,ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz  
      Complex *16 ys_yy,ys_yz,ys_zy,ys_zz          
      Complex *16 neta,alpha,Aiz,dAiz
      Complex *16 E_conj,dE_conj
      Complex *16 g_F,g_S,nxF,nxS
      Complex *16 gs_int(361501),gf_int(361501)
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)
 


      un=1.
      deux=2.
      pi=4.*DATAN(un)


 
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
    

!      Call eval_fem_gS_gF(ny,nz,g_S,g_F)
      Call eval_fem_gS(ny,nz,g_S)
      Call eval_fem_gF(ny,nz,g_F) 

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

      Call estim_gS_gF(ny,nz,g_S,g_F)

      ys_yy= -(0.,1.)*ny*nz*g_F/
     &               ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes)
     &       +(0.,1.)*nz*D/
     &               ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes)

      ys_yz=-(0.,1.)*((1.,0.)*S+(0.,1.)*pertes)*g_S/
     &                      ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)
     &      -(0.,1.)*ny**2.*nz**2.*g_F/
     &              (((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)*
     &                    ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes))
     &      +(0.,1.)*ny*nz**2.*D/
     &              (((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)*
     &                    ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes))
               
      ys_zy=+(0.,1.)*((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)*g_F/
     &                      ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)
     &      -(0.,1.)*ny*D/
     &                      ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)
     &                         
 
      ys_zz=+(0.,1.)*ny*nz*g_F/
     &                      ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)
     &      -(0.,1.)*ny**2.*nz*D/
     &              (((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)*
     &                    ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes))


 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      EndIf 
    


      If (int_sp.EQ.2) Then

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

         
 8    Return
      End
      
      




 
 
c **************************************************************************************************
c **************************************************************************************************
c int�gration


      Real *8  function real_f(ny)

      implicit none

      Common/quad/nz_quad,error,pertes,quad_key

      Integer *4 int_sp,num,quad_key    
      Real *8 ny,nz,nz_quad,error,pertes
      Complex *16 f_ny_nz,f_ny_nz_conj
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz 

      int_sp=1
      num=1
      nz=nz_quad
      Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                      hy_ny_nz,hz_ny_nz,int_sp,num)    
      f_ny_nz_conj=CONJG(f_ny_nz)
      real_f=(f_ny_nz+f_ny_nz_conj)/2.

      Return
      End


      Real *8 function imag_f(ny)

      implicit none

      Common/quad/nz_quad,error,pertes,quad_key

      Integer *4 int_sp,num,quad_key    
      Real *8 ny,nz,nz_quad,error,pertes
      Complex *16 f_ny_nz,f_ny_nz_conj
      Complex *16 ey_ny_nz,ez_ny_nz,hy_ny_nz,hz_ny_nz 

      int_sp=1
      num=1
      nz=nz_quad
      Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                      hy_ny_nz,hz_ny_nz,int_sp,num)    
      f_ny_nz_conj=CONJG(f_ny_nz)
      imag_f=-(0.,1.)*(f_ny_nz-f_ny_nz_conj)/2.

      Return
      End




      Real *8 function real_F_nz(nz)

      implicit none

      Common/plasma/S,D,dD,P,dP    
      Common/quad/nz_quad,error,pertes,quad_key
 
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_ny,b_ny,nz,nz_quad,real_f,epsabs,epsrel,points(2)
      Real *8 abserr,work(2000),S,D,dD,P,dP,tamp,epsabsp,epsrelp 
      Real *8 borne_max,demi_intervalle,error,pertes

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
      Common/quad/nz_quad,error,pertes,quad_key
 
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_ny,b_ny,nz,nz_quad,imag_f,epsabs,epsrel,points(2)
      Real *8 abserr,work(2000),S,D,dD,P,dP,tamp,epsabsp,epsrelp
      Real *8 borne_max,demi_intervalle,error,pertes

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
      Common/quad/nz_quad,error,pertes,quad_key
    
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      integer *4 quad_key
      Real *8 a_nz,b_nz,epsabs,epsrel,epsabsp,epsrelp,real_F_nz
      Real *8 Kr,points(2),abserr,work(2000),S,D,dD,P,dP,K_tamp 
      Real *8 borne_max,demi_intervalle,error,nz_quad,pertes

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
      Common/quad/nz_quad,error,pertes,quad_key
    
      Integer *4 npts2,key,neval,ier,limit,leniw,lenw,last,iwork(802)
      Integer *4 quad_key
      Real *8 a_nz,b_nz,epsabs,epsrel,epsabsp,epsrelp,imag_F_nz
      Real *8 Ki,points(2),abserr,work(2000),S,D,dD,P,dP,K_tamp 
      Real *8 borne_max,demi_intervalle,error,nz_quad,pertes

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

  
      Real *8 Kr,Ki
      Complex *16 K

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

      Integer *4 nbre,n_1,n_2,incr,int_sp   
      Real *8 min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Real *8 pas_ny,pas_nz,ny,nz
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
            Call chps_spect(ny,nz,f_ny_nz,ey_ny_nz,ez_ny_nz,
     &                      hy_ny_nz,hz_ny_nz,int_sp,incr) 

            eyt_ny_nz(incr)=ey_ny_nz
            ezt_ny_nz(incr)=ez_ny_nz
            hyt_ny_nz(incr)=hy_ny_nz
            hzt_ny_nz(incr)=hz_ny_nz

            incr=incr+1    
   
 50      Continue  
 
 40   Continue  
        
      Return
      End









c **************************************************************************************************
c **************************************************************************************************       
c calcul numerique de gS et gF



      Subroutine eval_fem_gS(ny,nz,g_S)


      implicit none

      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/quad/nz_quad,error,pertes,quad_key
       
      Integer *4 nbre_intervalles,info,ind_piv(2000)
      Integer *4 incr,incr_l,incr_c,quad_key
      Real *8 k0,Y0
      Real *8 S,D,dD,P,dP,nz_quad,error,pertes
      Real *8 P_tab(1001),M1(2001,2001),M2(2001,2001),L(2001,2001)
      Real *8 M1_ref_l1(3,3),M1_ref_l2(3,3),M2_ref(3,3),L_ref(3,3)
      Real *8 ny,nz,Long,Long_0,Long_1,Long_norm,h,h_norm
      Real *8 Eps_eq,dEps_eq,un,pi,nbre_per,nbre_eva
      Complex *16 A(7,2000),B(2000,1)
      Complex *16 g_S,K

      un=1.
      pi=4.*DATAN(un)

      nbre_per=25
      nbre_eva=10

      Eps_eq = (P/S)*(S-nz**2.)-ny**2.


      If ((dP/S)*(S-nz**2.).GT.0.) then
      
         If (Eps_eq.GT.0.) Then

            Long=nbre_per*2*pi/(k0*DSQRT(Eps_eq))
            Eps_eq = (P/S)*(S-nz**2.)-ny**2.+
     &                   (dP/S)*(S-nz**2.)*k0*Long
 10         If (Long.GT.1.5*nbre_per*2*pi/(k0*DSQRT(Eps_eq))) Then
               Long=Long/1.5
               Eps_eq = (P/S)*(S-nz**2.)-ny**2.+
     &                   (dP/S)*(S-nz**2.)*k0*Long 
               GOTO 10
            EndIf
            nbre_intervalles=500
            dEps_eq = (dP/S)*(S-nz**2.)
            K=(1.,0.)*dEps_eq/Eps_eq+(0.,1.)*DSQRT(Eps_eq)

         Else

            Long_0=(ny**2.-(P/S)*(S-nz**2.))/
     &                        ((dP/S)*(S-nz**2.)*k0)
            If (DSQRT(-Eps_eq).GT.nbre_eva/(k0*Long_0)) Then
               Long=nbre_eva/(k0*DSQRT(-Eps_eq))
               nbre_intervalles=100
               K = (0.,0.)
            Else
               Long_1=((nbre_per*2.*pi/DSQRT((dP/S)*(S-nz**2.)))
     &                                                   **(2./3.))/k0
               Long=Long_0+Long_1
               Eps_eq = (P/S)*(S-nz**2.)-ny**2.+
     &                          (dP/S)*(S-nz**2.)*k0*Long 
               dEps_eq = (dP/S)*(S-nz**2.)
               K=(1.,0.)*dEps_eq/Eps_eq+(0.,1.)*DSQRT(Eps_eq)
               nbre_intervalles=MIN(1000,
     &                              MAX(600,100*(1+INT(Long_1/Long_0))))

            EndIf

         Endif


      Else
     
         If (Eps_eq.LT.0.) Then

            Long=nbre_eva/(k0*DSQRT(-Eps_eq))
            Eps_eq = (P/S)*(S-nz**2.)-ny**2.+
     &                   (dP/S)*(S-nz**2.)*k0*Long
 20         If (Long.GT.1.5*nbre_eva/(k0*DSQRT(-Eps_eq))) Then
               Long=Long/1.5
               Eps_eq = (P/S)*(S-nz**2.)-ny**2.+
     &                   (dP/S)*(S-nz**2.)*k0*Long 
               GOTO 20
            EndIf
            nbre_intervalles=100
            K=(0.,0.)

         Else

            Long_0=(ny**2.-(P/S)*(S-nz**2.))/
     &                        ((dP/S)*(S-nz**2.)*k0)

            Long_1=((nbre_eva/DSQRT(-(dP/S)*(S-nz**2.)))**(2./3.))/k0
            Long=Long_0+Long_1
            K=(0.,0.)
            nbre_intervalles=MIN(1000,
     &                           MAX(600,100*(1+INT(Long_0/Long_1))))

          EndIf

  
      Endif


c$$$      Long=0.1
c$$$      nbre_intervalles=1000
c$$$         
c$$$      If ((dP/S)*(S-nz**2.).GT.0.) Then
c$$$        Eps_eq=(P/S)*(S-nz**2.)-ny**2.+
c$$$     &                   (dP/S)*(S-nz**2.)*k0*Long
c$$$        dEps_eq=(dP/S)*(S-nz**2.)
c$$$        If (Eps_eq.GT.0.) Then
c$$$          K=(1.,0.)*dEps_eq/Eps_eq+(0.,1.)*DSQRT(Eps_eq)
c$$$        Else
c$$$          K=(0.,0.)
c$$$        Endif
c$$$      Else
c$$$          K=(0.,0.)
c$$$      Endif    


      h=Long/nbre_intervalles
      Long_norm=k0*Long
      h_norm=k0*h
      
      Do 60 incr=1,nbre_intervalles+1
         P_tab(incr)=P+dP*h_norm*(incr-1)
60    Continue


      Do 70 incr_c=1,2*nbre_intervalles+1
         Do 80 incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles+1,incr_c+2)
            M1(incr_l,incr_c)=0.
            M2(incr_l,incr_c)=0.
            L(incr_l,incr_c)=0.
80       Continue
70    Continue


      M1_ref_l1(1,1)=(h_norm/60.)*7.
      M1_ref_l1(1,2)=(h_norm/60.)*4.    
      M1_ref_l1(1,3)=(h_norm/60.)*(-1)
      M1_ref_l1(2,1)=(h_norm/60.)*4.
      M1_ref_l1(2,2)=(h_norm/60.)*16.
      M1_ref_l1(2,3)=(h_norm/60.)*0.
      M1_ref_l1(3,1)=(h_norm/60.)*(-1.)
      M1_ref_l1(3,2)=(h_norm/60.)*0.
      M1_ref_l1(3,3)=(h_norm/60.)*1.

      M1_ref_l2(1,1)=(h_norm/60.)*1.
      M1_ref_l2(1,2)=(h_norm/60.)*0.    
      M1_ref_l2(1,3)=(h_norm/60.)*(-1.)
      M1_ref_l2(2,1)=(h_norm/60.)*0.
      M1_ref_l2(2,2)=(h_norm/60.)*16.
      M1_ref_l2(2,3)=(h_norm/60.)*4.
      M1_ref_l2(3,1)=(h_norm/60.)*(-1.)
      M1_ref_l2(3,2)=(h_norm/60.)*4.
      M1_ref_l2(3,3)=(h_norm/60.)*7.

      M2_ref(1,1)=(h_norm/15.)*2.
      M2_ref(1,2)=(h_norm/15.)*1.   
      M2_ref(1,3)=(h_norm/15.)*(-1./2.)
      M2_ref(2,1)=(h_norm/15.)*1.
      M2_ref(2,2)=(h_norm/15.)*8.
      M2_ref(2,3)=(h_norm/15.)*1.
      M2_ref(3,1)=(h_norm/15.)*(-1./2.)
      M2_ref(3,2)=(h_norm/15.)*1.
      M2_ref(3,3)=(h_norm/15.)*2.
 
      L_ref(1,1)=(1./(3.*h_norm))*7.
      L_ref(1,2)=(1./(3.*h_norm))*(-8.)   
      L_ref(1,3)=(1./(3.*h_norm))*1.
      L_ref(2,1)=(1./(3.*h_norm))*(-8.)
      L_ref(2,2)=(1./(3.*h_norm))*16.
      L_ref(2,3)=(1./(3.*h_norm))*(-8.)
      L_ref(3,1)=(1./(3.*h_norm))*1.
      L_ref(3,2)=(1./(3.*h_norm))*(-8.)
      L_ref(3,3)=(1./(3.*h_norm))*7.


      Do 90 incr=1,nbre_intervalles

         M1(2*(incr-1)+1,2*(incr-1)+1)=M1(2*(incr-1)+1,2*(incr-1)+1)
     &      +M1_ref_l1(1,1)*P_tab(incr)+M1_ref_l2(1,1)*P_tab(incr+1)
         M1(2*(incr-1)+1,2*(incr-1)+2)=M1(2*(incr-1)+1,2*(incr-1)+2)
     &      +M1_ref_l1(1,2)*P_tab(incr)+M1_ref_l2(1,2)*P_tab(incr+1)
         M1(2*(incr-1)+1,2*(incr-1)+3)=M1(2*(incr-1)+1,2*(incr-1)+3)
     &      +M1_ref_l1(1,3)*P_tab(incr)+M1_ref_l2(1,3)*P_tab(incr+1)
         M1(2*(incr-1)+2,2*(incr-1)+1)=M1(2*(incr-1)+2,2*(incr-1)+1)
     &      +M1_ref_l1(2,1)*P_tab(incr)+M1_ref_l2(2,1)*P_tab(incr+1)
         M1(2*(incr-1)+2,2*(incr-1)+2)=M1(2*(incr-1)+2,2*(incr-1)+2)
     &      +M1_ref_l1(2,2)*P_tab(incr)+M1_ref_l2(2,2)*P_tab(incr+1)
         M1(2*(incr-1)+2,2*(incr-1)+3)=M1(2*(incr-1)+2,2*(incr-1)+3)
     &      +M1_ref_l1(2,3)*P_tab(incr)+M1_ref_l2(2,3)*P_tab(incr+1)   
         M1(2*(incr-1)+3,2*(incr-1)+1)=M1(2*(incr-1)+3,2*(incr-1)+1)
     &      +M1_ref_l1(3,1)*P_tab(incr)+M1_ref_l2(3,1)*P_tab(incr+1)
         M1(2*(incr-1)+3,2*(incr-1)+2)=M1(2*(incr-1)+3,2*(incr-1)+2)
     &      +M1_ref_l1(3,2)*P_tab(incr)+M1_ref_l2(3,2)*P_tab(incr+1)      
         M1(2*(incr-1)+3,2*(incr-1)+3)=M1(2*(incr-1)+3,2*(incr-1)+3)
     &      +M1_ref_l1(3,3)*P_tab(incr)+M1_ref_l2(3,3)*P_tab(incr+1)


         M2(2*(incr-1)+1,2*(incr-1)+1)=M2(2*(incr-1)+1,2*(incr-1)+1)
     &                                +M2_ref(1,1)
         M2(2*(incr-1)+1,2*(incr-1)+2)=M2(2*(incr-1)+1,2*(incr-1)+2)
     &                                +M2_ref(1,2)
         M2(2*(incr-1)+1,2*(incr-1)+3)=M2(2*(incr-1)+1,2*(incr-1)+3)
     &                                +M2_ref(1,3)
         M2(2*(incr-1)+2,2*(incr-1)+1)=M2(2*(incr-1)+2,2*(incr-1)+1)
     &                                +M2_ref(2,1)
         M2(2*(incr-1)+2,2*(incr-1)+2)=M2(2*(incr-1)+2,2*(incr-1)+2)
     &                                +M2_ref(2,2)
         M2(2*(incr-1)+2,2*(incr-1)+3)=M2(2*(incr-1)+2,2*(incr-1)+3)
     &                                +M2_ref(2,3)  
         M2(2*(incr-1)+3,2*(incr-1)+1)=M2(2*(incr-1)+3,2*(incr-1)+1)
     &                                +M2_ref(3,1)
         M2(2*(incr-1)+3,2*(incr-1)+2)=M2(2*(incr-1)+3,2*(incr-1)+2)
     &                                +M2_ref(3,2)    
         M2(2*(incr-1)+3,2*(incr-1)+3)=M2(2*(incr-1)+3,2*(incr-1)+3)
     &                                +M2_ref(3,3)

         L(2*(incr-1)+1,2*(incr-1)+1)=L(2*(incr-1)+1,2*(incr-1)+1)
     &                                +L_ref(1,1)
         L(2*(incr-1)+1,2*(incr-1)+2)=L(2*(incr-1)+1,2*(incr-1)+2)
     &                                +L_ref(1,2)
         L(2*(incr-1)+1,2*(incr-1)+3)=L(2*(incr-1)+1,2*(incr-1)+3)
     &                                +L_ref(1,3)
         L(2*(incr-1)+2,2*(incr-1)+1)=L(2*(incr-1)+2,2*(incr-1)+1)
     &                                +L_ref(2,1)
         L(2*(incr-1)+2,2*(incr-1)+2)=L(2*(incr-1)+2,2*(incr-1)+2)
     &                                +L_ref(2,2)
         L(2*(incr-1)+2,2*(incr-1)+3)=L(2*(incr-1)+2,2*(incr-1)+3)
     &                                +L_ref(2,3)  
         L(2*(incr-1)+3,2*(incr-1)+1)=L(2*(incr-1)+3,2*(incr-1)+1)
     &                                +L_ref(3,1)
         L(2*(incr-1)+3,2*(incr-1)+2)=L(2*(incr-1)+3,2*(incr-1)+2)
     &                                +L_ref(3,2)    
         L(2*(incr-1)+3,2*(incr-1)+3)=L(2*(incr-1)+3,2*(incr-1)+3)
     &                                +L_ref(3,3)



90    Continue



      Do 100 incr_c=1,2*nbre_intervalles

         Do 110 incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles,incr_c+2)

            A(5+incr_l-incr_c,incr_c) = -(1.,0.)*L(incr_l+1,incr_c+1) 
     &             + ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/
     &                 ((1.,0.)*S+(0.,1.)*pertes)*
     &                      ((1.,0.)*M1(incr_l+1,incr_c+1)+
     &                          (0.,1.)*pertes*M2(incr_l+1,incr_c+1))
     &                    - (1.,0.)*(ny**2.)*M2(incr_l+1,incr_c+1)

110      Continue                                      

         B(incr_c,1)=(0.,0.)

100    Continue


       B(1,1)=(1.,0.)*L(2,1)
     &              -((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/
     &                   ((1.,0.)*S+(0.,1.)*pertes)*
     &                      ((1.,0.)*M1(2,1)+(0.,1.)*pertes*M2(2,1))
     &              + (1.,0.)*(ny**2.)*M2(2,1)
       B(2,1)=(1.,0.)*L(3,1)
     &               -((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/
     &                    ((1.,0.)*S+(0.,1.)*pertes)*
     &                      ((1.,0.)*M1(3,1)+(0.,1.)*pertes*M2(3,1)) 
     &              + (1.,0.)*(ny**2.)*M2(3,1)
  
  
      A(5,2*nbre_intervalles)=A(5,2*nbre_intervalles)+K

      Call ZGBSV(2*nbre_intervalles,2,2,1,A,7,ind_piv,B,
     &                                  2*nbre_intervalles,info)

      g_S=(1/h_norm)*((-3.,0)+4.*B(1,1)-B(2,1))


      Return
      End



      Subroutine eval_fem_gF(ny,nz,g_F)


      implicit none

      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/quad/nz_quad,error,pertes,quad_key
       
      Integer *4 nbre_intervalles,info,ind_piv(2000)
      Integer *4 incr,incr_l,incr_c,quad_key
      Real *8 k0,Y0
      Real *8 S,D,dD,P,dP,nz_quad,error,pertes
      Real *8 D_carre_tab(1001)
      Real *8 M1_bis(2001,2001),M2(2001,2001),L(2001,2001)
      Real *8 M1_ref_l1(3,3),M1_ref_l2(3,3),M2_ref(3,3),L_ref(3,3)
      Real *8 ny,nz,Long,Long_0,Long_1,Long_norm,h,h_norm
      Real *8 Eps_eq,dEps_eq,un,pi,nbre_per,nbre_eva
      Complex *16 A_bis(7,2000),B_bis(2000,1)
      Complex *16 g_F,K

      nbre_per=25
      nbre_eva=10

      un=1.
      pi=4.*DATAN(un)

      If ((D/dD).LT.0.) then
         write(*,*) 'warning D/dD<0 non traite dans g_f'
      Endif


      Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((D/dD)**2.)


      If ((nz**2.).GT.S) Then
      
         If (Eps_eq.GT.0.) Then

            Long=nbre_per*2*pi/(k0*DSQRT(Eps_eq))
            Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)   
 10         If (Long.GT.1.5*nbre_per*2*pi/(k0*DSQRT(Eps_eq))) Then
               Long=Long/1.5
               Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)  
               GOTO 10
            EndIf
            nbre_intervalles=500
            dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
            K=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*DSQRT(Eps_eq)

         Else

            Long_0=(DSQRT((ny*dD/(S-nz**2.)-(S-ny**2.-nz**2.))/ 
     &            (dD**2./(nz**2.-S)))-D/dD)/k0
            If (DSQRT(-Eps_eq).GT.nbre_eva/(k0*Long_0)) Then
               Long=nbre_eva/(k0*DSQRT(-Eps_eq))
               nbre_intervalles=100
               K = (0.,0.)
            Else
               Long_1=((nbre_per*2*pi/
     &                     DSQRT((2*(dD**2.)/(nz**2.-S))*
     &                             (k0*Long_0+D/dD)))**(2./3.))/k0
               Long=Long_0+Long_1
               Eps_eq=S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.) 
 15            If (Long_1.GT.1.5*nbre_per*2*pi/(k0*DSQRT(Eps_eq))) Then
                  Long_1=Long_1/1.5
                  Long=Long_0+Long_1
                  Eps_eq=S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &                     (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)  
                  GOTO 15
               EndIf
               dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
               K=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*DSQRT(Eps_eq)
               nbre_intervalles=MIN(1000,
     &                              MAX(600,100*(1+INT(Long_1/Long_0))))
            EndIf

         Endif


      Else
     
         If (Eps_eq.LT.0.) Then

            Long=nbre_eva/(k0*DSQRT(-Eps_eq))
            Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.) 
 20         If (Long.GT.1.5*nbre_eva/(k0*DSQRT(-Eps_eq))) Then
               Long=Long/1.5
               Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.) 
               GOTO 20
            EndIf
            nbre_intervalles=100
            K=(0.,0.)

         Else

            Long_0=(DSQRT((ny*dD/(S-nz**2.)-(S-ny**2.-nz**2.))/ 
     &            (dD**2./(nz**2.-S)))-D/dD)/k0
            Long_1=((nbre_eva/
     &              DSQRT((-2*(dD**2.)/(nz**2.-S))*
     &                       (k0*Long_0+D/dD)))**(2./3.))/k0
            Long=Long_0+Long_1
            Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.) 
 25         If (Long_1.GT.1.5*nbre_eva/(k0*DSQRT(-Eps_eq))) Then
               Long_1=Long_1/1.5
               Long=Long_0+Long_1
               Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
     &                      (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)  
               GOTO 25
            EndIf
            K=(0.,0.)
            nbre_intervalles=MIN(1000,
     &                           MAX(600,100*(1+INT(Long_0/Long_1))))

          EndIf
  
      Endif


c$$$      Long=0.2
c$$$      nbre_intervalles=1000
c$$$         
c$$$      If ((nz**2.).GT.S) Then
c$$$         Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ 
c$$$     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
c$$$         dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
c$$$         If (Eps_eq.GT.0.) Then
c$$$             K=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*DSQRT(Eps_eq)
c$$$         Else
c$$$             K=(0.,0.)
c$$$         Endif
c$$$      Else
c$$$         K = (0.,0.)
c$$$      Endif

      h=Long/nbre_intervalles
      Long_norm=k0*Long
      h_norm=k0*h
      
      Do 60 incr=1,nbre_intervalles+1
         D_carre_tab(incr)=(D+dD*h_norm*(incr-1))**2.
60    Continue


      Do 70 incr_c=1,2*nbre_intervalles+1
         Do 80 incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles+1,incr_c+2)
            M1_bis(incr_l,incr_c)=0.
            M2(incr_l,incr_c)=0.
            L(incr_l,incr_c)=0.
80       Continue
70    Continue


      M1_ref_l1(1,1)=(h_norm/60.)*7.
      M1_ref_l1(1,2)=(h_norm/60.)*4.    
      M1_ref_l1(1,3)=(h_norm/60.)*(-1)
      M1_ref_l1(2,1)=(h_norm/60.)*4.
      M1_ref_l1(2,2)=(h_norm/60.)*16.
      M1_ref_l1(2,3)=(h_norm/60.)*0.
      M1_ref_l1(3,1)=(h_norm/60.)*(-1.)
      M1_ref_l1(3,2)=(h_norm/60.)*0.
      M1_ref_l1(3,3)=(h_norm/60.)*1.

      M1_ref_l2(1,1)=(h_norm/60.)*1.
      M1_ref_l2(1,2)=(h_norm/60.)*0.    
      M1_ref_l2(1,3)=(h_norm/60.)*(-1.)
      M1_ref_l2(2,1)=(h_norm/60.)*0.
      M1_ref_l2(2,2)=(h_norm/60.)*16.
      M1_ref_l2(2,3)=(h_norm/60.)*4.
      M1_ref_l2(3,1)=(h_norm/60.)*(-1.)
      M1_ref_l2(3,2)=(h_norm/60.)*4.
      M1_ref_l2(3,3)=(h_norm/60.)*7.

      M2_ref(1,1)=(h_norm/15.)*2.
      M2_ref(1,2)=(h_norm/15.)*1.   
      M2_ref(1,3)=(h_norm/15.)*(-1./2.)
      M2_ref(2,1)=(h_norm/15.)*1.
      M2_ref(2,2)=(h_norm/15.)*8.
      M2_ref(2,3)=(h_norm/15.)*1.
      M2_ref(3,1)=(h_norm/15.)*(-1./2.)
      M2_ref(3,2)=(h_norm/15.)*1.
      M2_ref(3,3)=(h_norm/15.)*2.
 
      L_ref(1,1)=(1./(3.*h_norm))*7.
      L_ref(1,2)=(1./(3.*h_norm))*(-8.)   
      L_ref(1,3)=(1./(3.*h_norm))*1.
      L_ref(2,1)=(1./(3.*h_norm))*(-8.)
      L_ref(2,2)=(1./(3.*h_norm))*16.
      L_ref(2,3)=(1./(3.*h_norm))*(-8.)
      L_ref(3,1)=(1./(3.*h_norm))*1.
      L_ref(3,2)=(1./(3.*h_norm))*(-8.)
      L_ref(3,3)=(1./(3.*h_norm))*7.


      Do 90 incr=1,nbre_intervalles

         M1_bis(2*(incr-1)+1,2*(incr-1)+1)=
     &               M1_bis(2*(incr-1)+1,2*(incr-1)+1)
     &                     +M1_ref_l1(1,1)*D_carre_tab(incr)
     &                             +M1_ref_l2(1,1)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+1,2*(incr-1)+2)=
     &               M1_bis(2*(incr-1)+1,2*(incr-1)+2)
     &                      +M1_ref_l1(1,2)*D_carre_tab(incr)
     &                             +M1_ref_l2(1,2)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+1,2*(incr-1)+3)=
     &               M1_bis(2*(incr-1)+1,2*(incr-1)+3)
     &                      +M1_ref_l1(1,3)*D_carre_tab(incr)
     &                             +M1_ref_l2(1,3)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+2,2*(incr-1)+1)=
     &               M1_bis(2*(incr-1)+2,2*(incr-1)+1)
     &                      +M1_ref_l1(2,1)*D_carre_tab(incr)
     &                             +M1_ref_l2(2,1)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+2,2*(incr-1)+2)=
     &               M1_bis(2*(incr-1)+2,2*(incr-1)+2)
     &                      +M1_ref_l1(2,2)*D_carre_tab(incr)
     &                            +M1_ref_l2(2,2)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+2,2*(incr-1)+3)=
     &               M1_bis(2*(incr-1)+2,2*(incr-1)+3)
     &                      +M1_ref_l1(2,3)*D_carre_tab(incr)
     &                            +M1_ref_l2(2,3)*D_carre_tab(incr+1)   
         M1_bis(2*(incr-1)+3,2*(incr-1)+1)=
     &               M1_bis(2*(incr-1)+3,2*(incr-1)+1)
     &                      +M1_ref_l1(3,1)*D_carre_tab(incr)
     &                            +M1_ref_l2(3,1)*D_carre_tab(incr+1)
         M1_bis(2*(incr-1)+3,2*(incr-1)+2)=
     &               M1_bis(2*(incr-1)+3,2*(incr-1)+2)
     &                      +M1_ref_l1(3,2)*D_carre_tab(incr)
     &                            +M1_ref_l2(3,2)*D_carre_tab(incr+1)      
         M1_bis(2*(incr-1)+3,2*(incr-1)+3)=
     &               M1_bis(2*(incr-1)+3,2*(incr-1)+3)
     &                      +M1_ref_l1(3,3)*D_carre_tab(incr)
     &                            +M1_ref_l2(3,3)*D_carre_tab(incr+1)

         M2(2*(incr-1)+1,2*(incr-1)+1)=M2(2*(incr-1)+1,2*(incr-1)+1)
     &                                +M2_ref(1,1)
         M2(2*(incr-1)+1,2*(incr-1)+2)=M2(2*(incr-1)+1,2*(incr-1)+2)
     &                                +M2_ref(1,2)
         M2(2*(incr-1)+1,2*(incr-1)+3)=M2(2*(incr-1)+1,2*(incr-1)+3)
     &                                +M2_ref(1,3)
         M2(2*(incr-1)+2,2*(incr-1)+1)=M2(2*(incr-1)+2,2*(incr-1)+1)
     &                                +M2_ref(2,1)
         M2(2*(incr-1)+2,2*(incr-1)+2)=M2(2*(incr-1)+2,2*(incr-1)+2)
     &                                +M2_ref(2,2)
         M2(2*(incr-1)+2,2*(incr-1)+3)=M2(2*(incr-1)+2,2*(incr-1)+3)
     &                                +M2_ref(2,3)  
         M2(2*(incr-1)+3,2*(incr-1)+1)=M2(2*(incr-1)+3,2*(incr-1)+1)
     &                                +M2_ref(3,1)
         M2(2*(incr-1)+3,2*(incr-1)+2)=M2(2*(incr-1)+3,2*(incr-1)+2)
     &                                +M2_ref(3,2)    
         M2(2*(incr-1)+3,2*(incr-1)+3)=M2(2*(incr-1)+3,2*(incr-1)+3)
     &                                +M2_ref(3,3)

         L(2*(incr-1)+1,2*(incr-1)+1)=L(2*(incr-1)+1,2*(incr-1)+1)
     &                                +L_ref(1,1)
         L(2*(incr-1)+1,2*(incr-1)+2)=L(2*(incr-1)+1,2*(incr-1)+2)
     &                                +L_ref(1,2)
         L(2*(incr-1)+1,2*(incr-1)+3)=L(2*(incr-1)+1,2*(incr-1)+3)
     &                                +L_ref(1,3)
         L(2*(incr-1)+2,2*(incr-1)+1)=L(2*(incr-1)+2,2*(incr-1)+1)
     &                                +L_ref(2,1)
         L(2*(incr-1)+2,2*(incr-1)+2)=L(2*(incr-1)+2,2*(incr-1)+2)
     &                                +L_ref(2,2)
         L(2*(incr-1)+2,2*(incr-1)+3)=L(2*(incr-1)+2,2*(incr-1)+3)
     &                                +L_ref(2,3)  
         L(2*(incr-1)+3,2*(incr-1)+1)=L(2*(incr-1)+3,2*(incr-1)+1)
     &                                +L_ref(3,1)
         L(2*(incr-1)+3,2*(incr-1)+2)=L(2*(incr-1)+3,2*(incr-1)+2)
     &                                +L_ref(3,2)    
         L(2*(incr-1)+3,2*(incr-1)+3)=L(2*(incr-1)+3,2*(incr-1)+3)
     &                                +L_ref(3,3)



90    Continue



      Do 100 incr_c=1,2*nbre_intervalles

         Do 110 incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles,incr_c+2)

           A_bis(5+incr_l-incr_c,incr_c) = -(1.,0.)*L(incr_l+1,incr_c+1) 
     &           + ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes))
     &                                    *M1_bis(incr_l+1,incr_c+1) 
     &           + ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes-
     &                (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes))
     &                                    *M2(incr_l+1,incr_c+1)

110      Continue                                      

         B_bis(incr_c,1)=(0.,0.)

100    Continue

       B_bis(1,1)=(1.,0.)*L(2,1)
     &     - ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes))*M1_bis(2,1)
     &       - ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes-
     &       (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes))*M2(2,1)

       B_bis(2,1)=(1.,0.)*L(3,1)
     &     - ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes))*M1_bis(3,1)
     &       - ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes-
     &      (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes))*M2(3,1)


 
       
       A_bis(5,2*nbre_intervalles)=A_bis(5,2*nbre_intervalles)+K


      Call ZGBSV(2*nbre_intervalles,2,2,1,A_bis,7,ind_piv,B_bis,
     &                                  2*nbre_intervalles,info)

      g_F=(1/h_norm)*((-3.,0)+4.*B_bis(1,1)-B_bis(2,1))

      Return
      End



  

  


      Subroutine map_gS_gF

      implicit none

      Common/admittance_surface/gs_int,gf_int,ys_yy_stok_sp,
     &ys_yz_stok_sp,ys_zy_stok_sp,ys_zz_stok_sp


      Integer *4 n1,n2,incr
      Real *8 ny,nz
      Complex *16 g_S,g_F,gs_int(361501),gf_int(361501)
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)

      incr=1
      
      Do 120 n1=0,1200 

         nz=0.0+n1*0.02
	 
	 Do 130 n2=0,300
	 
            ny=-3.0+n2*0.02

            Call eval_fem_gS(ny,nz,g_S)
            Call eval_fem_gF(ny,nz,g_F)
            gs_int(incr)=g_S
            gf_int(incr)=g_F

            incr=incr+1	 
 
 130      Continue  

c          Write(*,*) incr

 120   Continue  

      Return
      End



      Subroutine estim_gS_gF(ny,nz,g_S,g_F)

      implicit none

      Common/admittance_surface/gs_int,gf_int,ys_yy_stok_sp,
     &ys_yz_stok_sp,ys_zy_stok_sp,ys_zz_stok_sp


      Integer *4 n1,n2
      Real *8 ny,nz,abs_nz
      Real *8 d1,d2,d3,d4
      Complex *16 f1,f2,f3,f4
      Complex *16 g_S,g_F,gs_int(361501),gf_int(361501)
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)

      abs_nz=DABS(nz)

      n1=INT(abs_nz/0.02)+1
      n2=INT((ny+3.)/0.02)+1

      f1=gs_int(n2+(n1-1)*301)
      f2=gs_int(n2+1+(n1-1)*301)
      f3=gs_int(n2+1+n1*301)
      f4=gs_int(n2+n1*301)

      d1=DSQRT((abs_nz-(n1-1)*0.02)**2.+(ny-(n2-1)*0.02+3.)**2.)
      d2=DSQRT((abs_nz-(n1-1)*0.02)**2.+(ny-n2*0.02+3.)**2.)
      d3=DSQRT((abs_nz-n1*0.02)**2.+(ny-n2*0.02+3.)**2.)      
      d4=DSQRT((abs_nz-n1*0.02)**2.+(ny-(n2-1)*0.02+3.)**2.)


      g_S=(d2*d3*d4*f1 + d1*d3*d4*f2 + d1*d2*d4*f3 + d1*d2*d3*f4)/
     &        (d2*d3*d4 + d1*d3*d4 + d1*d2*d4 + d1*d2*d3)  

      f1=gf_int(n2+(n1-1)*301)
      f2=gf_int(n2+1+(n1-1)*301)
      f3=gf_int(n2+1+n1*301)
      f4=gf_int(n2+n1*301)


      g_F=(d2*d3*d4*f1 + d1*d3*d4*f2 + d1*d2*d4*f3 + d1*d2*d3*f4)/
     &        (d2*d3*d4 + d1*d3*d4 + d1*d2*d4 + d1*d2*d3)  



      Return
      End
