
      program couplage_2D

c **************************************************************************************************
c **************************************************************************************************
c Declaration de variable

      implicit none
      
c      200 guides au maximum
c	4 modes au maximum


      Common/dim/a1,b1,y1,z1,a2,b2,y2,z2
      Common/mode/m1,n1,TE_TM_1,m2,n2,TE_TM_2
      Common/plasma/S,D,dD,P,dP
      Common/onde/k0,Y0
      Common/pi/pi
      Common/integral/min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Common/admittance_surface/ys_yy_stok_int,ys_yz_stok_int,
     &ys_zy_stok_int,ys_zz_stok_int,ys_yy_stok_sp,ys_yz_stok_sp,
     &ys_zy_stok_sp,ys_zz_stok_sp
      Common/tab_admittance/ tab_rempli
      
      Integer *4 m1,n1,TE_TM_1,m2,n2,TE_TM_2,i,j,incr,ind_ref,ind
      Integer *4 nbre_guides,nbre_modes,T_grill,class_modes
      Integer *4 tab_ind_m(200*2*4*4),tab_ind_n(200*2*4*4)
      Integer *4 tab_TE_TM(200*2*4*4)
      Integer *4 tab_class(2*4*4),tab_TE_TM_tot(200*4)
      Integer *4 tab_ind_m_tot(200*4),tab_ind_n_tot(200*4)
      Integer *4 abs_ligne,abs_col,ind_g_col,ind_g_ligne
      Integer *4 ind_ligne,ind_col,ind_m,ind_n,ind_guide
      Integer *4 quot_i,quot_j,quot,rest_i,rest_j
      Integer *4 abs_tab
      Integer *4 tab_ligne(200*4*200*4),tab_col(200*4*200*4)
      Integer *4 ind_tab,K_known
      Integer *4 tab_TE_TM_1(200*4*200*4)
      Integer *4 tab_TE_TM_2(200*4*200*4)
      Integer *4 tab_rempli
c GNU
      Real *8 nbre_guides_tamp,nbre_modes_tamp
c
      Real *8 a(200),b(200),y(200),z(200),kc_carre(200*2*4*4)
      Real *8 tab_m1(200*4*200*4),tab_m2(200*4*200*4)
      Real *8 tab_n1(200*4*200*4),tab_n2(200*4*200*4)
      Real *8 tab_a1(200*4*200*4),tab_a2(200*4*200*4)
      Real *8 tab_b1(200*4*200*4),tab_b2(200*4*200*4)
      Real *8 tab_dy(200*4*200*4),tab_dz(200*4*200*4)
      Real *8 pi,cl,qe,me,ne,dne,wpe,dwpe,wce,B0,Eps0,mu0,k0,Y0,f
      Real *8 mD,nD,dnD,wpD,dwpD,wcD,mT,nT,dnT,wpT,dwpT,wcT
      Real *8 mH,nH,dnH,wpH,dwpH,wcH,mHe,nHe,dnHe,wpHe,dwpHe,wcHe
      Real *8 a1,b1,y1,z1,a2,b2,y2,z2
      Real *8 S,D,dD,P,dP
      Real *8 min_ny,min_nz,max_ny,max_nz,nbre_ny,nbre_nz
      Complex *16 K,K_cpl(200*4,200*4),Zhe(200*4)
      Complex *16 eyt_ny_nz(30000),ezt_ny_nz(30000)
      Complex *16 hyt_ny_nz(30000),hzt_ny_nz(30000)
      Complex *16 eytt_ny_nz(200*4,30000)
      Complex *16 eztt_ny_nz(200*4,30000)
      Complex *16 hytt_ny_nz(200*4,30000)
      Complex *16 hztt_ny_nz(200*4,30000)
      Complex *16 ys_yy_stok_int(160000),ys_yz_stok_int(160000)
      Complex *16 ys_zy_stok_int(160000),ys_zz_stok_int(160000)
      Complex *16 ys_yy_stok_sp(30000),ys_yz_stok_sp(30000)
      Complex *16 ys_zy_stok_sp(30000),ys_zz_stok_sp(30000)






c **************************************************************************************************
c **************************************************************************************************
c constantes universelles

      pi=4.*ATAN(1.)
      cl=3.E8
      Y0=1/(120.*pi)
      Eps0=(1/(36*pi))*1.E-9
      mu0=4.*pi*1E-7
      qe=1.6021E-19
      me=9.1091E-31 
      mD=3.345E-27
      mT=5.018E-27
      mH=1.673E-27
      mHe=5.018E-27   

 
 
c **************************************************************************************************
c **************************************************************************************************      
c variables du programme

      Open(10,File='par_grill_2D.dat',Status='unknown',
     &Form='formatted')

      Read(10,*) nbre_guides_tamp
      nbre_guides=nbre_guides_tamp
      Read(10,*) nbre_modes_tamp
      nbre_modes=nbre_modes_tamp
      Read(10,*) B0
      Read(10,*) f
      Read(10,*) ne
      Read(10,*) dne
      Read(10,*) (a(ind),ind=1,nbre_guides)
      Read(10,*) (b(ind),ind=1,nbre_guides)
      Read(10,*) (y(ind),ind=1,nbre_guides)
      Read(10,*) (z(ind),ind=1,nbre_guides)
      Read(10,*) min_ny
      Read(10,*) min_nz
      Read(10,*) max_ny
      Read(10,*) max_nz
      Read(10,*) nbre_ny
      Read(10,*) nbre_nz
      
      Rewind 10
      Close(10)


c class_modes = 0 TE10 + TM1n
c class_modes = 1 modes par ordre croissant de kc

      class_modes=0
      
      

      nD=0.9*ne
      dnD=0.9*dne
      nT=0.
      dnT=0.
      nH=0.1*ne
      dnH=0.1*dne
      nHe=0.
      dnHe=0.


      k0=2.*pi*f/cl
      

      wce=qe*B0/me   
      wcD=qe*B0/mD
      wcT=qe*B0/mT
      wcH=qe*B0/mH
      wcHe=2*qe*B0/mHe

      wpe=sqrt(ne*qe**2/(Eps0*me))      
      wpD=sqrt(nD*qe**2/(Eps0*mD))
      wpT=sqrt(nT*qe**2/(Eps0*mT))
      wpH=sqrt(nH*qe**2/(Eps0*mH))
      wpHe=sqrt(nHe*(2*qe)**2/(Eps0*mHe))  
      
      dwpe=sqrt(dne*qe**2/(Eps0*me))      
      dwpD=sqrt(dnD*qe**2/(Eps0*mD))
      dwpT=sqrt(dnT*qe**2/(Eps0*mT))
      dwpH=sqrt(dnH*qe**2/(Eps0*mH))
      dwpHe=sqrt(dnHe*(2*qe)**2/(Eps0*mHe))  

      S=1-( wpe**2/((2*pi*f)**2-wce**2) +
     &      wpD**2/((2*pi*f)**2-wcD**2) +
     &      wpT**2/((2*pi*f)**2-wcT**2) +
     &      wpH**2/((2*pi*f)**2-wcH**2) +
     &      wpHe**2/((2*pi*f)**2-wcHe**2) )      
      D=-wpe**2*wce/(2*pi*f*((2*pi*f)**2 - wce**2)) +
     &   wpD**2*wcD/(2*pi*f*((2*pi*f)**2 - wcD**2)) +
     &   wpT**2*wcT/(2*pi*f*((2*pi*f)**2 - wcT**2)) +
     &   wpH**2*wcH/(2*pi*f*((2*pi*f)**2 - wcH**2)) +
     &   wpHe**2*wcHe/(2*pi*f*((2*pi*f)**2 - wcHe**2))            
      dD=(-dwpe**2*wce/(2*pi*f*((2*pi*f)**2 - wce**2)) +
     &   dwpD**2*wcD/(2*pi*f*((2*pi*f)**2 - wcD**2)) +
     &   dwpT**2*wcT/(2*pi*f*((2*pi*f)**2 - wcT**2)) +
     &   dwpH**2*wcH/(2*pi*f*((2*pi*f)**2 - wcH**2)) +
     &   dwpHe**2*wcHe/(2*pi*f*((2*pi*f)**2 - wcHe**2)))
     &   /k0  
      P=1-( (wpe/(2*pi*f))**2 +
     &      (wpD/(2*pi*f))**2 +
     &      (wpT/(2*pi*f))**2 +
     &      (wpH/(2*pi*f))**2 +
     &      (wpHe/(2*pi*f))**2 )
     
      dP=-((dwpe/(2*pi*f))**2 +
     &   (dwpD/(2*pi*f))**2 +
     &   (dwpT/(2*pi*f))**2 +
     &   (dwpH/(2*pi*f))**2 +
     &   (dwpHe/(2*pi*f))**2)/k0

      

       
c **************************************************************************************************
c **************************************************************************************************
c programme
   
          
      Do 1 ind_guide=1,nbre_guides
      
         abs_ligne=(ind_guide-1)*nbre_modes+1
         tab_ind_m_tot(abs_ligne)=1
	 tab_ind_n_tot(abs_ligne)=0
         tab_TE_TM_tot(abs_ligne)=1
	 Zhe(abs_ligne)=(0.,1.)*2*pi*f*mu0/
     &	          CDSQRT((1.,0.)*((pi/a(ind_guide))**2-k0**2))
     
   
      
         Do 2 ind=1,nbre_modes-1
	 
	    abs_ligne=(ind_guide-1)*nbre_modes+ind+1
            tab_ind_m_tot(abs_ligne)=1
	    tab_ind_n_tot(abs_ligne)=ind
            tab_TE_TM_tot(abs_ligne)=2
	    
	    Zhe(abs_ligne)=CDSQRT((1.,0.)*((pi/a(ind_guide))**2+ 
     &                    (ind*pi/b(ind_guide))**2-k0**2))
     &                                       /((0.,1.)*2*pi*f*Eps0) 
	    
2	 Continue
      
1     Continue     
           
      If (class_modes.EQ.1) Then 
      
      Do 10 ind_guide=1,nbre_guides
       
         ind=1
       	       
         Do 20 ind_n=1,nbre_modes
	 
	    tab_ind_m(ind)=0
	    tab_ind_n(ind)=ind_n
	    kc_carre(ind)=(ind_n*pi/b(ind_guide))**2
            tab_TE_TM(ind)=1
	    
	    ind=ind+1	    

 20      Continue	 
	 
	 Do 30 ind_m=1,nbre_modes
	 
	    Do 40 ind_n=0,nbre_modes-1
	    
	       tab_ind_m(ind)=ind_m
	       tab_ind_n(ind)=ind_n
	       kc_carre(ind)=(ind_m*pi/a(ind_guide))**2+ 
     &                            (ind_n*pi/b(ind_guide))**2
               tab_TE_TM(ind)=1
	    
	    ind=ind+1	 
	  
 40	    Continue 
	  
 30	 Continue
 
 	 
	 
	 Do 50 ind_m=1,nbre_modes
	 
	    Do 60 ind_n=1,nbre_modes
	    
	       tab_ind_m(ind)=ind_m
	       tab_ind_n(ind)=ind_n
	       kc_carre(ind)=(ind_m*pi/a(ind_guide))**2+ 
     &                            (ind_n*pi/b(ind_guide))**2
               tab_TE_TM(ind)=2
	    
	    ind=ind+1	 
	  
 60	    Continue 
	  
 50	 Continue
 
         tab_class(1)=1 
 
         Do 70 ind=2,2*nbre_modes*nbre_modes
	 
	    tab_class(ind)=1
	 
	    Do 80 ind_ref=1,ind-1
       
	    
	       If (kc_carre(ind).LT.kc_carre(ind_ref)) Then
	       
	           tab_class(ind_ref)=tab_class(ind_ref)+1
	
	       Else   
	       
                   tab_class(ind)=tab_class(ind)+1
		   
	       Endif

 80	    Continue 
	  
 70	 Continue
 

 	 Do 90 ind=1,2*nbre_modes*nbre_modes
	 	   	 
	    If (tab_class(ind).LE.nbre_modes) Then
	    
	        abs_ligne=tab_class(ind)+nbre_modes*(ind_guide-1)	    
	        tab_ind_m_tot(abs_ligne)=tab_ind_m(ind)
		tab_ind_n_tot(abs_ligne)=tab_ind_n(ind)
		tab_TE_TM_tot(abs_ligne)=tab_TE_TM(ind)
	
		
		If (tab_TE_TM(ind).EQ.1) Then
		
		
		   Zhe(abs_ligne)=(0.,1.)*2*pi*f*mu0/
     &	                         CDSQRT((1.,0.)*(kc_carre(ind)-k0**2))
     
                Else 	   
	
		   Zhe(abs_ligne)=CDSQRT((1.,0.)*(kc_carre(ind)-k0**2))/
     &	                                         ((0.,1.)*2*pi*f*Eps0)               	
		
                Endif
	     
	    Endif 

 90      Continue  
	 
 10   Continue 
 
      Endif     
         
      ind=1	 
      tab_rempli=0
 
      Do 100 ind_g_ligne=1,nbre_guides
      
         a1=a(ind_g_ligne)
	 b1=b(ind_g_ligne)
	 y1=y(ind_g_ligne)
	 z1=z(ind_g_ligne)
	       
			    
	 Do 110 ind_ligne=1,nbre_modes
	     
	    abs_ligne=ind_ligne+nbre_modes*(ind_g_ligne-1)
	     
	    m1=tab_ind_m_tot(abs_ligne)
            n1=tab_ind_n_tot(abs_ligne)
            TE_TM_1=tab_TE_TM_tot(abs_ligne)
	    
	    
	    Do 120 ind_g_col=1,nbre_guides
	    
	        a2=a(ind_g_col)
	        b2=b(ind_g_col)
	        y2=y(ind_g_col)
	        z2=z(ind_g_col)

	        Do 130 ind_col=1,nbre_modes

                   abs_col=ind_col+nbre_modes*(ind_g_col-1)
	 
		   		   
c r�ciprocit�
c		   If (abs_ligne.GT.abs_col) Then
c	    
c	                K_cpl(abs_ligne,abs_col)=K_cpl(abs_col,abs_ligne)
c			 
c		   Else
		   
		   m2=tab_ind_m_tot(abs_col)
                   n2=tab_ind_n_tot(abs_col)
                   TE_TM_2=tab_TE_TM_tot(abs_col)
		   
		   abs_tab=(abs_ligne-1)*nbre_guides*nbre_modes+abs_col
		   tab_ligne(abs_tab)=abs_ligne
		   tab_col(abs_tab)=abs_col
		   tab_m1(abs_tab)=m1
		   tab_m2(abs_tab)=m2
		   tab_n1(abs_tab)=n1
		   tab_n2(abs_tab)=n2
		   tab_TE_TM_1(abs_tab)=TE_TM_1
		   tab_TE_TM_2(abs_tab)=TE_TM_2
		   tab_a1(abs_tab)=a1
		   tab_a2(abs_tab)=a2
		   tab_b1(abs_tab)=b1
		   tab_b2(abs_tab)=b2
		   tab_dy(abs_tab)=y1-y2
		   tab_dz(abs_tab)=z1-z2
		   
		   			   
		   ind_tab=1
		   K_known=0
		   
		   Do While ((K_known.EQ.0).AND.(ind_tab.LT.abs_tab)) 
		   
		     If ((tab_m1(ind_tab).EQ.m1).AND.
     &		         (tab_m2(ind_tab).EQ.m2).AND.
     &		         (tab_n1(ind_tab).EQ.n1).AND.
     &		         (tab_n2(ind_tab).EQ.n2).AND.
     &		         (tab_TE_TM_1(ind_tab).EQ.TE_TM_1).AND.
     &		         (tab_TE_TM_2(ind_tab).EQ.TE_TM_2).AND.
     &          	 (tab_a1(ind_tab).EQ.a1).AND.
     &		         (tab_a2(ind_tab).EQ.a2).AND.
     &		         (tab_b1(ind_tab).EQ.b1).AND.
     &		         (tab_b2(ind_tab).EQ.b2).AND.
     &		         (DABS(tab_dy(ind_tab)-(y1-y2)).LT.1.E-6).AND.
     &		         (DABS(tab_dz(ind_tab)-(z1-z2)).LT.1.E-6)) Then


				     
			 K_known=1
 			 
		     Else
		     
		         ind_tab=ind_tab+1
			 			 
		     Endif
		     
		   End Do


		     
		   If (K_known.EQ.0) Then
		   	   	
                       Call Calcul_K(K,eyt_ny_nz,ezt_ny_nz,
     &                               hyt_ny_nz,hzt_ny_nz)	
                       K_cpl(abs_ligne,abs_col)=K
                       tab_rempli=1


                       if (abs_ligne.EQ.1) Then

      
                            Do 140 incr=1,(nbre_ny+1)*(nbre_nz+1)
                          
                               eytt_ny_nz(abs_col,incr)=eyt_ny_nz(incr)                              
                               eztt_ny_nz(abs_col,incr)=ezt_ny_nz(incr)
                               hytt_ny_nz(abs_col,incr)=hyt_ny_nz(incr)
                               hztt_ny_nz(abs_col,incr)=hzt_ny_nz(incr)

 140                        Continue


                       Endif
                          
                       Write(*,*) ind
		       ind=ind+1


		       
		   Else		   

                       K_cpl(abs_ligne,abs_col)=
     &  	            K_cpl(tab_ligne(ind_tab),tab_col(ind_tab))


                 
		       
		   Endif
		   
c		 Endif
	   
		 
			           
 130             Continue
 
 120         Continue
 
 110     Continue
 
 100  Continue


 
      Open(10,File='K_plasma.dat',Status='unknown'
     &                                 ,Form='unformatted')
  
      Write(10) ((K_cpl(i,j),i=1,nbre_modes*nbre_guides),
     &j=1,nbre_modes*nbre_guides)

     
      Write(10) (Zhe(i),i=1,nbre_modes*nbre_guides)
     

      Rewind 10
      Close(10)

      Open(20,File='Spect_plasma.dat',Status='unknown'
     &                                   ,Form='unformatted')
  
      Write(20) ((eytt_ny_nz(i,j),i=1,nbre_modes*nbre_guides),
     &j=1,(nbre_ny+1)*(nbre_nz+1))
      Write(20) ((eztt_ny_nz(i,j),i=1,nbre_modes*nbre_guides),
     &j=1,(nbre_ny+1)*(nbre_nz+1)) 
      Write(20) ((hytt_ny_nz(i,j),i=1,nbre_modes*nbre_guides),
     &j=1,(nbre_ny+1)*(nbre_nz+1))
      Write(20) ((hztt_ny_nz(i,j),i=1,nbre_modes*nbre_guides),
     &j=1,(nbre_ny+1)*(nbre_nz+1))

      Rewind 20
      Close(20)


      Open(30,File='Modes.dat',Status='unknown'
     &                                   ,Form='unformatted')
  
      Write(30) (a(i),i=1,nbre_guides)
      Write(30) (b(i),i=1,nbre_guides)
      Write(30) (y(i),i=1,nbre_guides)
      Write(30) (z(i),i=1,nbre_guides)
      Write(30) (tab_ind_m_tot(i),i=1,nbre_guides*nbre_modes)
      Write(30) (tab_ind_n_tot(i),i=1,nbre_guides*nbre_modes)
      Write(30) (tab_TE_TM_tot(i),i=1,nbre_guides*nbre_modes)

      Rewind 30
      Close(30)

c ecriture des resultats dans un fichier texte
      Open(20,File='K_plasma2.dat', Status='replace', Form='formatted')
c header
      write(20,*) '==================== ALOHA 2D result file ========'
      write(20,*) '\t\t K_cpl '
c ecriture des resultats en colonnes
      do i=1,nbre_modes*nbre_guides
        do j=1,nbre_modes*nbre_guides
            write(20,*) K_cpl(i,j)
        end do
      end do
      Close(20)

      Open(30,File='Zhe2.dat', Status='replace', Form='formatted')
c header
      write(30,*) '==================== ALOHA 2D result file ========'
      write(30,*) '\t\t Zhe '  
      do j=1,nbre_modes*nbre_guides
        write(30,*) Zhe(j)
      end do
      Close(30)

c   ecriture du spectre dans un fichier ascii
      Open(20, File='Spect_plasma2.dat', Status='replace', 
     &Form='formatted')
c header
      write(20,*) '==================== ALOHA 2D result file ========'
      write(20,*) '\t\t    nbre_modes, nbre_guides'
      Write(20,*) nbre_modes, nbre_guides
      Write(20,*) '\t\t    nbre_ny, nbre_nz'
      Write(20,*) nbre_ny, nbre_nz 
      Write(20,*) '\t\t    Ey, Ez, Hy, Hz'
      do i=1,nbre_modes*nbre_guides
          do j=1,(nbre_ny+1)*(nbre_nz+1)
              Write(20,*) eytt_ny_nz(i,j), eztt_ny_nz(i,j), 
     &hytt_ny_nz(i,j), hztt_ny_nz(i,j)
          end do
      end do

      Close(20)

      Open(30,File='Modes2.dat', Status='replace', 
     &Form='formatted')
      Write(30,*) nbre_guides, nbre_modes
      Write(30,*) (a(i),i=1,nbre_guides)
      Write(30,*) (b(i),i=1,nbre_guides)
      Write(30,*) (y(i),i=1,nbre_guides)
      Write(30,*) (z(i),i=1,nbre_guides)
      Write(30,*) (tab_ind_m_tot(i),i=1,nbre_guides*nbre_modes)
      Write(30,*) (tab_ind_n_tot(i),i=1,nbre_guides*nbre_modes)
      Write(30,*) (tab_TE_TM_tot(i),i=1,nbre_guides*nbre_modes)

      Rewind 30
      Close(30)

       
      End

      
      
