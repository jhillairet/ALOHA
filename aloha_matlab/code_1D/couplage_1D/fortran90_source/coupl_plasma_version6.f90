! ********************************************************************* 
! Programme de calcul de matrice de scattering
! Couplage Grill Plasma avec differents cas poloidaux
! ********************************************************************* 
program coupl_plasma_version6
      ! Physics and mathematics constants
      Use aloha_constants
      ! Simulation configuration variables
      ! TODO : this module should not exist. It exists in replacement of the F77 "Common", but it makes
      ! the code quite complicated to follow...
      !       Common /err/ epsabs
      !       Common /nlim/ nlimit
      !       Common /dims/ b,z
      !       Common /idem/ bcte
      !       Common /onde/ k0
      !       Common /plasma/ X0,D0
      !       Common /integration/ max_nz
      Use aloha_config, Only : epsabs, nlimit, b, z, bcte, k0, X0, D0, max_nz, Gmax, Nmhm, Nmem, X1, D1, d_couche, pertes, d_vide
      ! aloha main routines 
      Use aloha_sgrill

      Implicit None

      Integer       :: ios,i,j,intg,ni

      ! JH 03/10/08 : Pour lire les donnees ascii, on a besoin d'utiliser des variables reelles. 
      ! Lorsqu'il s'agit en réalité d'eniter, on devra les caster comme tels une fois les
      ! valeurs reelles lues.
      Real(kind=wp) :: T_grill_r, D_guide_max_r
      Integer       :: T_grill, D_guide_max
      Real(kind=wp) :: Pg_r,Nmh_r,Nme_r
      Integer       :: Pg, Nmh, Nme
      Integer :: ok

      Real(kind=wp) :: f,epsrel 
      Real(kind=wp) :: a!,b(60),z(60)!,bc!,b0!,e(60)!,ei!,ac
      Real(kind=wp) :: n0, nc!,D0,X0,max_nz
      Real(kind=wp) :: n1, dn1, dn0 ! V6

      Integer         , allocatable, dimension(:)   :: prh, pre
      Complex(kind=wp), allocatable, dimension(:,:) :: S, Spr
      Complex(kind=wp), allocatable, dimension(:,:) :: K_cpl
      Complex(kind=wp), allocatable, dimension(:,:) :: RZheg, RZhegpr
      
      Integer, parameter :: VERSION = 6


!       Open(10,File='par_grill.dat',Status='unknown',Form='formatted')
      Open(unit=10, file='par_grill.dat', status='old', form='formatted', action='read', position='rewind', iostat=ios) ! JH 17/04/2008 unknow->old
      if (ios /= 0) then
        write (*,*) "[ERROR] : erreur lecture fichier par_grill.dat"
        stop
      end if

      Read(10,*) Nmh_r
      Read(10,*) Nme_r
      Nmh = int(Nmh_r)
      Nme = int(Nme_r)
      Read(10,*) f
      Read(10,*) n0
      Read(10,*) dn0 !V6
      Read(10,*) d_couche !V6
      Read(10,*) dn1 ! V6
      Read(10,*) Pg_r ! Nombre total de guides sur une ligne poloidale
      Pg = int(Pg_r)
      Read(10,*) a
      Read(10,*) (b(i),i=1,Pg)
      Read(10,*) (z(i),i=1,Pg)
      Read(10,*) T_grill_r
      T_grill = int(T_grill_r)
      Read(10,*) D_guide_max_r
      D_guide_max = int(D_guide_max_r)
      Read(10,*) epsrel
      Read(10,*) pertes !V6
      Read(10,*) max_nz
      Read(10,*) d_vide !V6 couche de vide
      
      Rewind 10
      Close(10)

      ! copy read values to shared memory through aloha_config module variables :
      k0=2.*pi*f/cl  
      nc=(2.*pi*f)**2*me*Eps0/qe**2
      X0=n0/nc
      Nmhm=Nmh
      Nmem=Nme

      ! V6
      D0=k0*n0/dn0
      n1=n0+dn0*d_couche
      X1=n1/nc
      D1=k0*n1/dn1

        
      
      ! TODO : intg is not use anymore 
      intg=4 ! Mode d'integration numerique
      nlimit=0
      ni = 32 
      epsabs=0.

      ! Initialisation des matrices 
      allocate(    S((Nme+Nmh)*Pg,(Nme+Nmh)*Pg), stat=ok)
      allocate(K_cpl((Nme+Nmh)*Pg,(Nme+Nmh)*Pg), stat=ok)
      allocate(RZheg((Nme+Nmh)*Pg,(Nme+Nmh)*Pg), stat=ok)

      allocate(    Spr(Pg,Pg), stat=ok)
      allocate(RZhegpr(Pg,Pg), stat=ok)

      allocate(pre(Pg), prh(Pg), stat=ok)


      Call SMGrill(VERSION, intg,ni,epsrel,T_grill,D_guide_max, &
                    Nmh, Nme, Pg, a, prh, pre, & 
                    S, RZheg, Spr, RZhegpr, K_cpl)

! ------------------------------------------------------
      ! Sauvegarde des resultats pour matlab
      write(*,*) "Sauvegarde des resultats..."
      
      Open(10,File='S_plasma.dat', Status='replace', Form='unformatted')
 
      Write(10) ((real(S(i,j))     ,i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
      Write(10) ((aimag(S(i,j))    ,i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
      Write(10) ((real(RZheg(i,j)) ,i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
      Write(10) ((aimag(RZheg(i,j)),i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
      Write(10) ((real(K_cpl(i,j)) ,i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
      Write(10) ((aimag(K_cpl(i,j)),i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)

      Close(10)

      ! ecriture des resultats dans un fichier texte
      Open(20,File='S_plasma2.dat', Status='replace', Form='formatted')
      ! header
      write(20,*) '==================== ALOHA result file ===================='
      write(20,*) '\t\t S \t\t RZheg \t\t K_cpl '
      ! ecriture des resultats en colonnes
      do i=1,size(S,1)
        do j=1,size(S,2)
            write(20,*) S(i,j), RZheg(i,j), K_cpl(i,j)
        end do
      end do

! i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
!       Write(20,*) ((RZheg(i,j),i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)
!       Write(20,*) ((K_cpl(i,j),i=1,(Nme+Nmh)*Pg), j=1,(Nme+Nmh)*Pg)

      Close(20)
    
      Stop
end program coupl_plasma_version6
! Fin du programme principal