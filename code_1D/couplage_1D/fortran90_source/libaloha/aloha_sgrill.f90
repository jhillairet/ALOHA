module aloha_sgrill

    Use aloha_constants
    Use aloha_function
    Use aloha_scattering_matrix_utils

    ! NAG
    !Use nag_airy_fun, Only : nag_airy_ai 

    implicit none

! ********************************************************************* 
! Programme de calcul de matrice de scattering
! DIScontinuite GRILl PLasma tire du prg bersNSbR2.f
! valide (en spectre) %SWAN le 17 octdbre 1995
! ********************************************************************* 
! LAST UPDATES :
! V3 MAJ F90 JH 2008
! V6 MAJ F90 JH 07/2008

contains
    ! 
    ! 
    ! SMGrill
    ! 
    ! 
    Subroutine SMGrill(version, intg,ni,epsrel,T_grill,D_guide_max,Nmh,Nme,Pg,a,prh,pre,S,RZheg,Spr,RZhegpr,K_cpl)

    Use aloha_config, Only : Nmhm,Nmem,app,epsabs,nlimit,Nmax,Gmax,b,z,bcte,k0,X0,D0,max_nz,d_vide

    ! NAG
    !Use nag_gen_lin_sys, Only : nag_gen_lin_sol

    ! LAPACK
    Use LA_PRECISION
    Use F95_LAPACK, Only : LA_GESV

    implicit none

    !     Common/nbmod/Nmhm,Nmem
    !     Common/Airy/app
    !     Common/err/epsabs
    !     Common/nlim/nlimit
    !     Common/Nax/Nmax
    !     Common/dims/b,z
    !     Common/idem/bcte
    !     Common/onde/k0
    !     Common/plasma/X0,D0
    !     Common/integration/max_nz

    Integer,          intent(in)  :: version, intg, ni, T_grill, D_guide_max, Nmh, Nme, Pg
    Real(kind=wp),    intent(in)  :: a
    Integer,          intent(out) :: prh(Pg), pre(Pg)
    Complex(kind=wp), intent(out) :: S((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)
    Complex(kind=wp), intent(out) :: Spr(Pg,Pg), RZhegpr(Pg,Pg) 
    Complex(kind=wp), intent(out) :: K_cpl((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)

    ! ----------------------------------------------- 

    Integer :: i,j,k,l!,Nmheg,dNmaxg,ifail
    Integer :: mhi(2,10*10),mei(2,10*10)
    Integer :: prhi,prei,nprhi,nprei
    Integer :: rest_i,rest_l,quot,quot_i,quot_l,guide_j
    Integer :: nprh(Pg),npre(Pg),prhe
    Integer :: Statsym(Pg,Pg)
    Integer :: Statsymgg((Nmh+Nme)*Pg,(Nmh+Nme)*Pg)

    Real(kind=wp)    :: epsrel

    Complex(kind=wp) :: zohe(Nmh,Nme), zoeh(Nme,Nmh)
    Complex(kind=wp) :: RZh(Nmh,Nmh), RYh(Nmh,Nmh)
    Complex(kind=wp) :: RZe(Nme,Nme), RYe(Nme,Nme)
    Complex(kind=wp) :: RZhe(Nme+Nmh, Nme+Nmh), RZheg((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)
    Complex(kind=wp) :: Y(Nmh+Nme,(Nme+Nmh)*Pg)
    Complex(kind=wp) :: H((Nme+Nmh)*Pg,(Nme+Nmh)*Pg) 
    Complex(kind=wp) :: RZYRZg((Nme+Nmh),(Nme+Nmh)*Pg)
    Complex(kind=wp) :: Idg((Nme+Nmh)*Pg,(Nme+Nmh)*Pg), IdgpH((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)
    Complex(kind=wp) :: IdgpHI((Nme+Nmh)*Pg,(Nme+Nmh)*Pg), IdgmH((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)

    ! initialisation des variables
    Spr = 0.0
    RZhegpr = 0.0

    write(*,*) 'Verification parametres en entree de SMgrill :'
    write(*,*) 'intg=', intg
    write(*,*) 'ni=', ni
    write(*,*) 'epsrel=', epsrel      
    write(*,*) 'T_grill=', T_grill
    write(*,*) 'D_guide_max=', D_guide_max
    write(*,*) 'Nmh=', Nmh
    write(*,*) 'Nme=', Nme
    write(*,*) 'Pg=', Pg
    write(*,*) 'a=', a
!     write(*,*) 'prh=', prh
!     write(*,*) 'pre=', pre

    write(*,*) 'Variables globales'
!     write(*,*) 'Nmhm=', Nmhm
!     write(*,*) 'Nmem=', Nmem
    write(*,*) 'app=', app
    write(*,*) 'epsabs=', epsabs
    write(*,*) 'nlimit=', nlimit
    write(*,*) 'Nmax=', Nmax ! maximum |N_parallel|
    write(*,*) 'b(1:5)=', b(1:5)
    write(*,*) 'z(1:5)=', z(1:5)
    write(*,*) 'bcte=', bcte
    write(*,*) 'k0=', k0
    write(*,*) 'X0=', X0
    write(*,*) 'D0=', D0
    write(*,*) 'max_nz=',max_nz    
    write(*,*) 'd_vide=',d_vide

    Nmhm=Nmh
    Nmem=Nme

    ! Matrice identite JH
    Idg = eye((Nme+Nmh)*Pg)

    ! Calcul de la matrice d'impedance globale
    Call CalMatig(Pg,a,b,Nmh,Nme,RZheg) 

!     do i=1,(Nme+Nmh)*Pg
!         write(*,*) 'i=', i,' RZheg(i,i)=', RZheg(i,i)
!     end do
! ! OK

    ! Boucle sur les differents guides
    Do i=1,Pg
        ! calcule les impedances caracteristiques des guides selon les modes
        Call CalRZh(i,a,b,Nmh,RZh,RYh)
        Call CalRZe(i,a,b,Nme,RZe,RYe)

!         Call Zero(Nmh,Nme,zohe)
!         Call Zero(Nme,Nmh,zoeh)
!         Call Mpm(Nmh,Nmh,Nme,Nme,RZh,zohe,zoeh,RZe,RZhe)

        ! Concatenation des matrices RZ sous la forme d'une matrice bloc : 
        !        [ RZh   0  ]
        ! RZhe = [          ]
        !        [  0   RZe ]
        RZhe(:,:) = 0.0 ! on s'assure de l'initialisation 
        RZhe(1:Nmh, 1:Nmh) = RZh ! bloc 11
        RZhe(Nmh+1:Nme+Nmh, Nmh+1:Nme+Nmh) = RZe ! bloc 22
        
        !  Calcul du couplage guides-plasma     
        !  Retourne la matrice admittance plasma Y

        write(*,*) 'Calcul de la matrice admittance'
        Call CalcouPgi(version, intg,i,a,epsrel,T_grill,D_guide_max,Nmh,Nme,Pg,Y)         
        
! ! Calcul de RYhe^-1*(Y*RZheg)=RZhe*(Y*RZheg)=RZYRZg en fortran 77
!          Call F06ZAF(TransA,TransB,Nmh+Nme,(Nmh+Nme)*Pg,(Nmh+Nme)*Pg,
!      &Alpha,Y,2*Nmax,RZheg,2*Nmax*Gmax,Beta,YRZheg,2*Nmax)
! 
!          Call F06ZAF(TransA,TransB,Nmh+Nme,(Nmh+Nme)*Pg,Nmh+Nme,
!      &Alpha,RZhe,2*Nmax,YRZheg,2*Nmax,Beta,RZYRZg,2*Nmax)
        
        ! Calcul de RYhe^-1*(Y*RZheg)=RZhe*(Y*RZheg)=RZYRZg en fortran 90
        RZYRZg = matmul(RZhe, matmul(Y, RZheg))

        ! Constrution de la matrice H
        do k=1,Nmh+Nme
            do l=1,(Nmh+Nme)*Pg
                if (MOD(l,(Nmh+Nme)) == 0) then
                    guide_j=l/(Nmh+Nme)
                else
                    guide_j=1+(l-MOD(l,(Nmh+Nme)))/(Nmh+Nme)
                endif

                If (i > guide_j) then
                    H((Nmh+Nme)*(i-1)+k,l)=H(l,(Nmh+Nme)*(i-1)+k)
                    K_cpl((Nmh+Nme)*(i-1)+k,l)=K_cpl(l,(Nmh+Nme)*(i-1)+k)
                ElseIf ((i > T_grill).AND.(l > T_grill*(Nmh+Nme))) then
                    quot_i=(i-MOD(i,T_grill))/T_grill
                    quot_l=(l-MOD(l,T_grill*(Nmh+Nme)))/(T_grill*(Nmh+Nme))
                    quot=MIN(quot_i,quot_l)
                    rest_i=i-T_grill*quot
                    rest_l=l-T_grill*quot*(Nme+Nmh)
                    If ((rest_i == 0).OR.(rest_l == 0)) Then
                        rest_i=rest_i+T_grill
                        rest_l=rest_l+T_grill*(Nmh+Nme)
                    Endif
                    H((Nmh+Nme)*(i-1)+k,l)=H((Nmh+Nme)*(rest_i-1)+k,rest_l)	 
                    K_cpl((Nmh+Nme)*(i-1)+k,l)=K_cpl((Nmh+Nme)*(rest_i-1)+k,rest_l)	              
                Else
                    H((Nmh+Nme)*(i-1)+k,l)=RZYRZg(k,l) 
                    K_cpl((Nmh+Nme)*(i-1)+k,l)=Y(k,l)  
                endif
           end do ! l=1,(Nmh+Nme)*Pg
       end do ! k=1,Nmh+Nme
    end do ! i / Pg

    write(*,*) 'Calcul de la matrice de scattering...'
    ! Calcul de la matrice de Scattering S=(I+H)^(-1).(I-H) en fortran 90
    IdgpH = Idg + H
    IdgmH = Idg - H
    IdgpHI= eye((Nme+Nmh)*Pg)

! !     inversion de la matrice IdgpH (NAG)
! !     on inverse en realité le systeme lineaire IdgpH.X=IdgpHI, 
! !     ou IdgpHI est la matrice identite.
! !     L'argument IdgpHI est ecrase par la sortie IdgpH^(-1)
! !     
!     CALL nag_gen_lin_sol(IdgpH, IdgpHI) 

    ! inversion de la matrice IdgpH (LAPACK)
    ! on inverse en realité le systeme lineaire IdgpH.X=IdgpHI, 
    ! ou IdgpHI est la matrice identite.
    ! L'argument IdgpHI est ecrase par la sortie IdgpH^(-1)
    ! 
    call LA_GESV(IdgpH, IdgpHI)

    S = matmul(IdgpHI, IdgmH)

    ! Calcul du nombre de modes propageants - construction de la matrice
    ! le nombre de mode propageant est propre a chque guide
    ! Verification de l'unitarite et de la symetrie
    Do i=1,Pg
        Call CalNbmpr(i,Nmh,Nme,a,b,mhi,mei,prhi,prei,nprhi,nprei)
        prh(i)=prhi
        pre(i)=prei
        nprh(i)=nprhi
        npre(i)=nprei
    end do

    prhe=0
    Call ConstMatpr(Pg,prh,pre,Nmh,Nme,S,Spr)
    Call ConstMatpr(Pg,prh,pre,Nmh,Nme,RZheg,RZhegpr)

    Do i=1,Pg
        prhe=prhe+(prh(i)+pre(i))
    end do

    
    ! Calcul de la matrice amplitude et phase
    
    ! Integration 
    If (intg == 1) Then
        Write (*,*) 'AKF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
    ElseIf (intg == 2) Then
        Write (*,*) 'ARF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
    ElseIf (intg == 3) Then
        Write (*,*) 'BDF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
    ElseIf (intg == 4) Then
        Write (*,*) 'AHF/epsrel,Nmh,Nme=',epsrel,Nmh,Nme
    Else
        Write (*,*) 'BAF/ni,Nmh,Nme=',ni,Nmh,Nme
    EndIf

    Call ContSymgg((Nmh+Nme)*Pg,S,Statsymgg)
    Call ContSym1(prhe,Spr,Statsym)


    write(*,*) '---------Sortie de SMGrill-------'   
    
    write(*,*) 'sum(RZh)=', sum(RZh)
    write(*,*) 'sum(RZe)=', sum(RZe)
    write(*,*) 'sum(RZhe)=', sum(RZhe)

    write(*,*) 'sum(Y)=', sum(Y)
    write(*,*) 'sum(RZYRZg)=', sum(RZYRZg)
    write(*,*) 'sum(H)=',sum(H)

    write(*,*) 'sum(Spr)=',sum(Spr)
    write(*,*) 'sum(RZhegpr)=',sum(RZhegpr)
    write(*,*) 'sum(S)='    ,sum(S)
    write(*,*) 'sum(RZheg)=',sum(RZheg)
    write(*,*) 'sum(K_cpl)=',sum(K_cpl)

!      write (*,*) 'Statsym OK'
!      Call ContUni1(prhe,Spr,Statuni)
!      write (*,*) 'Statuni OK'

!      Stop
    
    End subroutine SMGrill


    ! CARACTERISTIQUE DU GRILL

    ! ***
    ! Subroutine CLASSant les MODEs TE par ordre croissant
    ! 
    ! INPUT
    !  - Nmh : number of TE modes
    !  - p : indice du guide
    !  - a : longueur grand cote guide
    !  - b : hauteur petit cote guide (1x?)
    ! OUPUT
    !  - mh  : mode indice (2,Nmh)
    !  - kch : cut off wavenumbers for the mode TEm0
    ! 
    ! ***
    Subroutine ClassModeTE(Nmh,p,a,b,mh,kch)
    
        use aloha_config, only : Gmax
    
    !     Common/nbmod/Nmhm,Nmem
    
        Integer      , Intent(in) :: Nmh,p
        Real(kind=wp), Intent(in) :: a
        Real(kind=wp), Intent(in), Dimension(:) :: b

        Integer      , Intent(out), Dimension(2,Nmh*Nmh) :: mh
        Real(kind=wp), Intent(out), Dimension(Nmh*Nmh)   :: kch


        Integer :: i,m,n

        i=1
    
    !      Do 100 m=0,Nmax-1
    !         Do 110 n=0,Nmax-1
    !            If (m.EQ.0.AND.n.EQ.0) Then
    !               Goto 110
    !            EndIf
    
    ! on force la polarisation en TEm0
        n=0
        Do m=1,Nmh
            kch(i)=sqrt(m**2*b(p)/a+n**2*a/b(p))*PI/sqrt(a*b(p)) ! JH 17/04/2008 sqrt 
            mh(1,i)=m
            mh(2,i)=n
            i=i+1
        end do 
    !      ri=i-1
    
        ! classement des modes se referant au rang de kch
        ! la troncature des Nmh premier terme ce fera ds les subroutines
        !      IFAIL=0
        !      Call M01DAF(kch,1,ri,'Ascending',IRANK,IFAIL)
        !      Call M01EAF(mh,1,ri,IRANK,IFAIL)
        ! classement des coupures par ordre croissant
        !      Call M01CAF(kch,1,ri,'Ascending',IFAIL) 
    
    End subroutine ClassModeTE


    ! ***
    ! Subroutine CLASSant les MODEs TM par ordre croissant
    ! 
    ! INPUT
    !  - Nmax : number of TM modes
    !  - p : indice du guide
    !  - a : longueur grand cote guide
    !  - b : hauteur petit cote guide (1x?)
    ! OUPUT
    !  - me  : mode indice (2,Nmh)
    !  - kce : cut off wavenumbers for the mode TEm0
    ! 
    ! ***
    Subroutine ClassModeTM(Nmax,p,a,b,me,kce)
    
        use aloha_config, only : Nmhm, Nmem
   
    !     Common/nbmod/Nmhm,Nmem
    
        Integer,       Intent(in) :: Nmax, p 
        Real(kind=wp), Intent(in) , dimension(:)            :: b
        Integer,       Intent(out), dimension(2,Nmax*Nmax)  :: me
        Real(kind=wp), Intent(out), dimension(Nmax*Nmax)        :: kce

        ! ----------------

        Integer       :: i, m, n!,Nmhm,Nmem!,ri
        Real(kind=wp) :: a, mn!,PI


        i=1
        mn=Nmem/Nmhm
    
    !      Do 100 m=0,Nmax-1
    !         Do 110 n=0,Nmax-1
    !            If (m.EQ.0.OR.n.EQ.0) Then
    !               Goto 110
    !            EndIf
    
        ! on ne prend que les TM avec m identique au m du TE
        Do m=1,Nmhm
    !         Do 110 n=m,(m+IDNINT(mn)-1)
            ! En multimode on peut avoir n < m car "a" augmente Ez/Ey#a/b.n/m
            Do n=1,Nint(mn)
                kce(i)=sqrt(m**2*b(p)/a+n**2*a/b(p))*PI/sqrt(a*b(p))
                me(1,i)=m
                me(2,i)=n
                i=i+1
            end do
        end do
    !      ri=i-1
    
    ! classement des modes se referant au rang de kce
    ! la troncature des Nme premier terme ce fera ds les subroutines
    !      IFAIL=0
    !      Call M01DAF(kce,1,ri,'Ascending',IRANK,IFAIL)
    !      Call M01EAF(me,1,ri,IRANK,IFAIL)
    ! classement des coupures par ordre croissant
    !      Call M01CAF(kce,1,ri,'Ascending',IFAIL) 
        
    
    End subroutine ClassModeTM




    Subroutine CalRZh(p,a,b,Nmh, RZhi,RYhi)
    
        use aloha_config, only : Nmax,Gmax, k0    
    
    !     Common/Nax/Nmax
    !     Common/onde/k0
        Integer          :: Nmh, Nme, i, p
        Integer          :: mh(2,Nmh*Nmh)
        Real(kind=wp)    :: a,kch(Nmh*Nmh)
        Complex(kind=wp) :: RZh(Nmh),RYh(Nmh),tampon,tamponbis
        Complex(kind=wp) :: RZhi(Nmh,Nmh),RYhi(Nmh,Nmh)
        Real(kind=wp), dimension(:) :: b
    
        Call ClassModeTE(Nmh,p,a,b,mh,kch)
   
        Do i=1,Nmh
    !         RZh(i)=csqrt(Z0/csqrt((1.,0.)*(1-(kch(i)/k0)**2)))
            tampon=(0.,1.)/((cl**2*Eps0)*sqrt((1.,0.)*(kch(i)**2-k0**2)))
            tamponbis=k0*cl*tampon
            RZh(i)=sqrt(tamponbis)
            RYh(i)=1/RZh(i)
        end do

        RZhi = diag(RZh)
        RYhi = diag(1/RZh)
        
        
    End subroutine CalRZh



    ! ***
    Subroutine CalRZe(p,a,b,Nme,RZei,RYei)
    
        use aloha_config, only : Nmax, Gmax, k0
    
    !     Common/Nax/Nmax
    !     Common/onde/k0
        Integer          :: Nme,i,p
        Integer          :: me(2,Nme*Nme)
        Real(kind=wp)    :: a,kce(Nme*Nme)!
        Complex(kind=wp) :: RZe(Nme),RYe(Nme),tampon
        Complex(kind=wp) :: RZei(Nme,Nme),RYei(Nme,Nme)
        Real(kind=wp), dimension(:) :: b
        
        Call ClassModeTM(Nmax,p,a,b,me,kce)
    
        Do i=1,Nme
    !         RZe(i)=csqrt(Z0*csqrt((1.,0.)*(1-(kce(i)/k0)**2)))
            tampon=sqrt((1.,0.)*(kce(i)**2-k0**2))/((0.,1.)*k0*cl*Eps0) 
            RZe(i)=sqrt(tampon)
            RYe(i)=1/RZe(i)
        end do

        RZei = diag(RZe)
        RYei = diag(1/RZe)
    
        
    End subroutine CalRZe

    ! ***
    ! CALcul de la MATrice d'Impedance Global
    ! ***
    Subroutine CalMatig(Pg,a,b,Nmh,Nme,RZheg)   
        use aloha_config, only : Nmax, Gmax
    
        Integer          :: Nmh,Nme,p,d,Pg,i,j
        Real(kind=wp)    :: a
        Complex(kind=wp) :: RZh(Nmh,Nmh),RYh(Nmh,Nmh) 
        Complex(kind=wp) :: RZe(Nme,Nme),RYe(Nme,Nme) 
        Complex(kind=wp) :: RZheg((Nme+Nmh)*Pg,(Nme+Nmh)*Pg)
        Real(kind=wp), dimension(:) :: b
    
        ! initialisation
        RZheg(:,:)=0.0
        
        d=0
        Do p=1,Pg
            Call CalRZh(p,a,b,Nmh,RZh,RYh)
            Call CalRZe(p,a,b,Nme,RZe,RYe)
    
            forall (i=1:Nmh, j=1:Nmh, i==j) RZheg(i+d,j+d)=RZh(i,j)
            forall (i=1:Nme, j=1:Nme, i==j) RZheg(i+d+Nmh,j+d+Nmh)=RZe(i,j)
    
            d=(Nmh+Nme)*p
        end do ! Pg
    
        forall(i=1:Pg*(Nmh+Nme), j=1:Pg*(Nmh+Nme), i==j) RZheg(i,j)=RZheg(i,j)
    
    
        
    End subroutine


! ***
! Sudroutine donnant la normalisation des fct. de forme sin et cos
! ***

    !
    ! INPUT
    !  - p : indice du guide
    !  - a : longueur grand cote du guide
    !  - b : largeur petit cote du guide
    !  - m : indice de mode m
    !  - n : indice de mode n
    ! 
    ! OUTPUT
    !  - E : 
    !  - sg: 
    !  
    ! coeff. de normalisation en modE TE 
    Subroutine CE(p,a,b,m,n,E,sg)
        use aloha_config, only : Nmax, Gmax
    
        Integer          :: m, n, p
        Real(kind=wp)    :: a, Epsm, Epsn, E
        Character(len=3) :: sg
        Real(kind=wp), dimension(:) :: b
    
        Call Eps(m,Epsm)
        Call Eps(n,Epsn)
    
        If ((m == 0) .AND. (n == 0)) Then
            Write (*,*) 'Pb m=0 et n=0 en Subroutine CE'
            Stop
        EndIf
        E=(sqrt(Epsm*Epsn)/sqrt(m**2*b(p)/a+n**2*a/b(p)))*m*m/a
    
    !      If (sg.EQ.'sgp') Then
    !         E=E
    !      Else
        If (sg == 'sgm') Then
            E=-1.*E
        EndIf

    End subroutine CE

    !
    ! INPUT
    !  - p : indice du guide
    !  - a : longueur grand cote du guide
    !  - b : largeur petit cote du guide
    !  - m : indice de mode m
    !  - n : indice de mode n
    ! 
    ! OUTPUT
    !  - E : 
    !  - sg: 
    !  
    ! coeff. de normalisation en moDe TM
    Subroutine CM(p,a,b,m,n,D,sg)   
        Integer          :: m, n, p
        Real(kind=wp)    :: a, D
        Character(len=3) :: sg
        Real(kind=wp), dimension(:) :: b
    
        If ((m == 0).OR.(n == 0)) Then
            Write (*,*) 'Pb m=0 ou n=0 en Subroutine CM'
            Stop
        EndIf
        D=(2./sqrt(m**2*b(p)/a+n**2*a/b(p)))*m*n/b(p)
    
    !      If (sg.EQ.'sgp') Then
    !         D=D
    !      Else
        If (sg == 'sgm') Then
            D=-1.*D
        EndIf

    End subroutine CM

    ! 
    !  Neumann factor :
    !  
    !  eps(m) = 1 if m = 0
    !         = 2 if m > 0
    ! 
    ! INPUT
    !  - m : indice (0 ou <>0)
    ! 
    ! OUPUT
    !  - Epsm : Neumann factor
    ! 
    Subroutine Eps(m, Epsm)
        Integer      , intent(in) :: m
        Real(kind=wp), intent(out):: Epsm

        If (m == 0) Then
            Epsm=1.0
        Else 
            Epsm=2.0
        EndIf
    End subroutine Eps

! ***
! Couplage au plasma
! ***

    ! ***
    ! Subroutine de CALcul de la matrice de COUplage Plasma/Guide 
    ! du guide d'indice p
    ! 
    ! Input arguments :
    !  - version [int] : version du calcul (3 ou 6)
    !  - intg [int] : mode d'integration 
    !  - p [int] : indice du guide
    !  - a [real] : largeur guide (dir. pol)
    !  - epsrel [real] : 
    !  - T_grill [int] : indice max du guide couple avec le guide p
    !  - D_guide_max [int] :
    !  - Nmh [int]: nombre de mode TE
    !  - Nme [int]: nombre de mode TM
    !  - Pg [int] : nombre total de guides
    ! 
    ! Output argument 
    !  - Y [Complex(Nmh+Nme,(Nme+Nmh)*Pg)] : admittance matrix
    ! 
    ! ***
    Subroutine CalcouPgi(version,intg,p,a,epsrel,T_grill,D_guide_max,Nmh,Nme,Pg,Y)

        use aloha_config, only : knout, Nmax, b, z, m, n, i, j, k0, X0, D0, X1, D1, d_couche, pertes, d_vide
        use aloha_function, only : Calresiduy
        use aloha_integration

        implicit none

        Integer,       intent(in) :: version, intg, p, T_grill, D_guide_max, Nmh, Nme, Pg
        Real(kind=wp), intent(in) :: a, epsrel
        Complex(kind=wp), Dimension(Nmh+Nme,(Nme+Nmh)*Pg), intent(out) :: Y
        

        Integer :: ct!knout,
        Integer :: k,l!Nmax,i,j
        Integer :: mi,mj!,ma,m,n
        Integer :: mhi(2,Nmh),mhj(2,Nmh)
        Integer :: mei(2,Nmh*Nme),mej(2,Nmh*Nme)
        !,b(60),z(60),k0,X0,D0,Z0,Y0
        Real(kind=wp) :: kch(Nmh),kce(Nmh*Nme)
        Real(kind=wp) :: absrr(Nmh*Nme*2 + Nme*Nme + Nmh*Nmh + 1), absrm
        Real(kind=wp) :: resultrhh,resultihh,resultrhe,resultihe
        Real(kind=wp) :: resultreh,resultieh,resultree,resultiee
        Real(kind=wp) :: cte,cofm,Di,Dj,Ei,Ej
        Character(len=3) :: sgp,sgm
    
        cte=-Y0*(X0/D0)**(1./3.)*(PI/a)**2./k0**4./((2*PI)**2.)
    
        sgp='sgp'
        sgm='sgm'
    
        j=p
        ! Pour tous les guides d'une ligne poloidale
        ! TODO : paralleliser ici !
        Do i=1,Pg
            !write(*,*) 'DEBUG: WG#',i,'/',Pg
            If (((i > T_grill).AND.(j > T_grill)).OR.(j > i).OR.(abs(i-j) > D_guide_max)) Then
                forall(k=1:Nmh, l=1:Nmh) Y(k,l+(Nmh+Nme)*(i-1))=(0.,0.)
                forall(k=1:Nmh, l=1:Nme) Y(k,l+Nmh+(Nmh+Nme)*(i-1))=(0.,0.)
                forall(l=1:Nme, k=1:Nmh) Y(Nmh+l,k+(Nmh+Nme)*(i-1))=(0.,0.)    
                forall(l=1:Nme, k=1:Nme) Y(Nmh+l,k+Nmh+(Nmh+Nme)*(i-1))=(0.,0.)
            Else
                ct=1
                Call ClassModeTE(Nmh,i,a,b,mhi,kch)
                Call ClassModeTE(Nmh,j,a,b,mhj,kch)
                Call ClassModeTM(Nmax,i,a,b,mei,kce)
                Call ClassModeTM(Nmax,j,a,b,mej,kce)

                Do k=1,Nmh
                    Do l=1,Nmh
                        Call CE(j,a,b,mhj(1,k),mhj(2,k),Ej,sgp)
                        Call CE(i,a,b,mhi(1,l),mhi(2,l),Ei,sgm)

                        mj=mhj(1,k)
                        mi=mhi(1,l)
                        If (mj == mi) Then
                            Call Calresiduy(a,mi,mj,cofm)
                        Else
                            cofm=0.0
                        EndIf

                        m=mhj(2,k)
                        n=mhi(2,l)
                        
                        !if (i.EQ.63) then
                        ! write(*,*) 'k=',k,' l=',l,'i=',i,' j=',j
                        !end if

                        if (cofm /= 0.0) then
                            if (version == 3) then
                                Call Int01AHF(epsrel,i,j,resultrhh,resultihh,absrm)

                            elseif (version == 6) then
                                Call Int01AJF(epsrel,i,j,resultrhh,resultihh,absrm)
                            endif
                        endif

                        !write(*,*) 'E-E :', resultrhh, resultihh
                        Y(k,l+(Nmh+Nme)*(i-1))=Ej*Ei*cofm*cte*(resultrhh+(0.,1.)*resultihh)
                        absrr(ct)=absrm
                        ct=ct+1
                    end do
                end do

                !if (i.EQ.63) then
                ! write(*,*) 'premier quadrant: sum(Y)=', sum(Y)
                !end if

                Do k=1,Nmh
                    Do l=1,Nme
                        Call CE(j,a,b,mhj(1,k),mhj(2,k),Ej,sgp)
                        Call CM(i,a,b,mei(1,l),mei(2,l),Di,sgm)
                        mj=mhj(1,k)
                        mi=mei(1,l)
                        If (mj == mi) Then
                            Call Calresiduy(a,mi,mj,cofm)
                        Else
                            cofm=0.0
                        EndIf
                        m=mhj(2,k)
                        n=mei(2,l)

                        if (cofm /= 0.0) then
                            if (version == 3) then
                                Call Int01AHF(epsrel,i,j,resultrhe,resultihe,absrm)
                            elseif (version == 6) then
                                Call Int01AJF(epsrel,i,j,resultrhe,resultihe,absrm)
                            endif
                        endif
                        !write(*,*) 'E-M :', resultrhe, resultihe
                        !write(*,*) Ej, Di, cofm, cte
                        
                        Y(k,l+Nmh+(Nmh+Nme)*(i-1))=Ej*Di*cofm*cte*(resultrhe+(0.,1.)*resultihe)
                            absrr(ct)=absrm
                            ct=ct+1
                    end do
                end do
                !if (i.EQ.63) then
                ! write(*,*) 'deuxieme quadrant : sum(Y)=',sum(Y)
                !end if
                Do l=1,Nme
                    Do k=1,Nmh
                        Call CM(j,a,b,mej(1,l),mej(2,l),Dj,sgp)
                        Call CE(i,a,b,mhi(1,k),mhi(2,k),Ei,sgm)
                        mj=mej(1,l)
                        mi=mhi(1,k)
                        If (mj == mi) Then
                            Call Calresiduy(a,mi,mj,cofm)
                        Else
                            cofm=0.0
                        EndIf
                        m=mej(2,l)
                        n=mhi(2,k)

                        if (cofm /= 0.0) then
                            if (version == 3) then
                                Call Int01AHF(epsrel,i,j,resultreh,resultieh,absrm)
                            elseif (version == 6) then
                                Call Int01AJF(epsrel,i,j,resultreh,resultieh,absrm)
                            endif
                        endif
                        !write(*,*) 'M-E :', resultreh, resultieh
                        Y(Nmh+l,k+(Nmh+Nme)*(i-1))=Dj*Ei*cofm*cte*(resultreh+(0.,1.)*resultieh)
                        absrr(ct)=absrm
                        ct=ct+1
                    end do
                end do
                !if (i.EQ.63) then
                ! write(*,*) 'troisieme quadrant : sum(Y)=',sum(Y)
                !end if
                Do l=1,Nme
                    Do k=1,Nme
                        Call CM(j,a,b,mej(1,l),mej(2,l),Dj,sgp)
                        Call CM(i,a,b,mei(1,k),mei(2,k),Di,sgm)
                        mj=mej(1,l)
                        mi=mei(1,k)
                        If (mj == mi) Then
                            Call Calresiduy(a,mi,mj,cofm)
                        Else
                            cofm=0.0
                        EndIf
                        m=mej(2,l)
                        n=mei(2,k)

                        if (cofm /= 0.0) then
                            if (version == 3) then
                                Call Int01AHF(epsrel,i,j,resultree,resultiee,absrm)
                            elseif (version == 6) then
                                Call Int01AJF(epsrel,i,j,resultree,resultiee,absrm)
                            endif
                        endif
                        !write(*,*) 'M-M :', resultree, resultiee
                        Y(Nmh+l,k+Nmh+(Nmh+Nme)*(i-1))=Dj*Di*cofm*cte*(resultree+(0.,1.)*resultiee)
                        
                        absrr(ct)=absrm
                        ct=ct+1
                    end do
                end do
                !      ct=ct-1
                Call tri(ct,absrr,absrm)
            Endif
        end do ! Pg

        !  DEBUG      
        !if (p<=T_grill) then
        write(*,*) 'DEBUG : final sum(Y)=',sum(Y)
        !endif
    End subroutine CalcouPgi


    !
    !+ tri
    !
    Subroutine tri(ct,absrr,absrm)
        implicit none

        Integer       :: i,ct
        Real(kind=wp) :: absrr(:), absrm
    
        absrm=minval(absrr)
!         Do i=1,ct
!             If (absrm.LT.absrr(i+1)) Then
!                 absrm=absrr(i+1)
!             EndIf
!         end do
!     Return
    End subroutine tri



    ! ***
    ! Subroutine creant une matrice DIAGOnale a partir d un vect. colonne
    ! ***
    ! JH 24/04/2008 : Obsolete : renmplace par la fonction diag(v) du module aloha_function
    Subroutine Diago(Mc,Nm,Md)   
        Integer :: Nm,i,j
        Complex(kind=wp) :: Mc(10),Md(10,10)
    
        Do j=1,Nm
            Do i=1,Nm
                If (i.EQ.j) Then
                Md(i,j)=Mc(i)
                Else
                Md(i,j)=(0.,0.)
                EndIf
            end do
        end do
    Return
    End subroutine


    ! ***
    ! Subroutine donnant une matrice de ZEROs 
    ! OBSOLETE : fortran 90 permet d'initialiser les tableaux directement par tab=0
    ! ***
    Subroutine Zero(Mm,Nm,zox)

    Integer :: i,j,Mm,Nm
    Complex(kind=wp) :: zox(10,10)

    Do j=1,Nm
        Do i=1,Mm
            zox(i,j)=(0.,0.)
        end do
    end do

    Return
    End subroutine


    ! ***
    ! Subroutine de construction d une Matrice Par Morceau 
    ! OBSOLETE : Fortran 90 permet de manipuler les tableaux par indice facilement
    ! ***
    Subroutine Mpm(Nm11,Nm12,Nm21,Nm22,Mor11,Mor12,Mor21,Mor22,Mat)
        Integer :: i,j,Nm11,Nm12,Nm21,Nm22
        Complex(kind=wp) :: Mor11(10,10),Mor12(10,10)
        Complex(kind=wp) :: Mor21(10,10),Mor22(10,10)
        Complex(kind=wp) :: Mat(2*10,2*10)

        Do j=1,Nm11
            Do i=1,Nm12
                Mat(i,j)=Mor11(i,j)
            end do
        end do
    
        Do j=1,Nm22
            do i=1,Nm12 
                Mat(i,Nm11+j)=Mor12(i,j)
            end do
        end do
    
        Do  j=1,Nm11
            Do i=1,Nm21
                Mat(Nm12+i,j)=Mor21(i,j)
            end do
        end do
    
        Do j=1,Nm22
            Do i=1,Nm21
                Mat(Nm12+i,Nm11+j)=Mor22(i,j)
            end do
        end do

    
    End subroutine Mpm

    ! ***
    ! Subroutine de CALcul  du NomBre de Mode PRopageant
    ! ***
    Subroutine CalNbmpr(i,Nmh,Nme,a,b,mhi,mei,prhi,prei,nprhi,nprei)
        use aloha_config, only : Nmax, k0
         !     Common/Nax/Nmax
        !     Common/onde/k0
       
        Integer :: i,l,Nmh,Nme,prhi,prei,nprhi,nprei!,Nmax
        Integer :: mhi(2,10*10),mei(2,10*10)
        Real(kind=wp) :: a,b(60)
        Real(kind=wp) :: kch(10*10),kce(10*10)!k0,
    
        !&&&      Integer :: i,l,Nmax,Nmh,Nme,prhi,prei,nprhi,nprei
        !&&&      Integer :: mhi(2,3*3),mei(2,3*3)
        !&&&      Real(kind=wp) :: a,b(15)
        !&&&      Real(kind=wp) :: k0,kch(3*3),kce(3*3)
    
        Call ClassModeTE(Nmh,i,a,b,mhi,kch)
        prhi=0
        Do l=1,Nmh
            If (k0.GT.kch(l)) Then
                prhi=prhi+1
            EndIf
        end do
        nprhi=0
        Do l=1,Nmax
            If (k0.GT.kch(l)) Then
                nprhi=nprhi+1
            EndIf
        end do
    
        Call ClassModeTM(Nmax,i,a,b,mei,kce)
        prei=0
        Do l=1,Nme
            If (k0.GT.kce(l)) Then
                prei=prei+1
            EndIf
        end do
        nprei=0
        Do l=1,Nmax
            If (k0.GT.kce(l)) Then
                nprei=nprei+1
            EndIf
        end do
    
    End subroutine CalNbmpr

    ! ***
    ! Subroutine de CONSTruction de la MATrice des modes PRopageants
    ! ***
    Subroutine ConstMatpr(Pg,prh,pre,Nmh,Nme,S,Spr)

        implicit none
    ! ccc
    ! modif altair.partenaires
    
        Integer :: p,q,i,j,Pg,deci,Ndeci,decj,Ndecj
        Integer :: Nmh,Nme,prh(Nmh),pre(Nme)
        Complex(kind=wp) :: S((Nme+Nmh)*Pg,(Nme+Nmh)*Pg),Spr(Pg,Pg)
    
    !&&&      Integer :: p,q,i,j,Pg,deci,Ndeci,decj,Ndecj
    !&&&      Integer :: Nmh,Nme,prh(15),pre(15)
    !&&&      Complex(kind=wp) :: S(2*3*15,2*3*15),Spr(2*3/2*15,2*3/2*15)
    
    ! Construction de Sprhh
        Spr(:,:) = 0! initalisation

        deci=0
        Ndeci=0
        Do p=1,Pg
            decj=0
            Ndecj=0
            Do q=1,Pg
                Do i=1,prh(p)
                    Do j=1,prh(q)
                        Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
                    end do
                end do
                decj=decj+prh(q)+pre(q)
                Ndecj=Ndecj+Nmh+Nme
            end do
            deci=deci+prh(p)+pre(p)
            Ndeci=Ndeci+Nmh+Nme
        end do
    
    ! Construction de Sprhe
        deci=0
        Ndeci=0
        Do p=1,Pg
            decj=prh(1)
            Ndecj=Nmh
            Do q=1,Pg
                Do i=1,prh(p)
                    Do j=1,pre(q)
                        Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
                    end do
                end do
                decj=decj+pre(q)+prh(q+1)
                Ndecj=Ndecj+Nme+Nmh
            end do
            deci=deci+prh(p)+pre(p)
            Ndeci=Ndeci+Nmh+Nme
        end do
    
    ! Construction de Spreh
        deci=prh(1)
        Ndeci=Nmh
        Do p=1,Pg
            decj=0
            Ndecj=0
            Do q=1,Pg
                Do i=1,pre(p)
                    Do j=1,prh(q)
                        Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
                    end do
                end do
                decj=decj+prh(q)+pre(q)
                Ndecj=Ndecj+Nmh+Nme
            end do
            deci=deci+pre(p)+prh(p+1)
            Ndeci=Ndeci+Nme+Nmh
        end do
    
    ! Construction de Spree
        deci=prh(1)
        Ndeci=Nmh
        Do p=1,Pg
            decj=prh(1)
            Ndecj=Nmh
            Do q=1,Pg
                Do i=1,pre(p)
                    Do j=1,pre(q)
                        Spr(i+deci,j+decj)=S(i+Ndeci,j+Ndecj)
                    end do
                end do
                decj=decj+pre(q)+prh(q+1)
                Ndecj=Ndecj+Nme+Nmh
            end do
            deci=deci+pre(p)+prh(p+1)
            Ndeci=Ndeci+Nme+Nmh
        end do
            
    End subroutine ConstMatpr

end module aloha_sgrill
