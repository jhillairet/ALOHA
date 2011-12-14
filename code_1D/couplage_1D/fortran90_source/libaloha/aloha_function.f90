! ***
! Couplage au plasma
! ***
MODULE aloha_function

      Use aloha_constants
      Use airy_functions_complex_double, Only: airy_ai
      ! TODO : la librairie utilisee ne permet que de calculer Airy Bi en reel !
      ! Il faudra trouver de quoi calculer Bi pour des valeurs complexes.
      Use airy_functions_real_double, Only: airy_bi

!       ! Equivalent NAG F90
!       Use nag_airy_fun, Only: nag_airy_ai

      implicit none

CONTAINS

!     !
!     !+ Classe les modes TM
!     !
!     function kc_e(Nmax,p,a,b,me)
! 
!         use aloha_config, only : Nmhm, Nmem
! 
!         ! input arguments
!         integer, intent(in)       :: Nmax, p, me(:,:)
!         real(kind=wp), intent(in) :: a, b(:)
!         real(kind=wp)             :: kc_e(Nmax,Nmax)
!         ! working variables
!         integer                   :: i,m,n
!         real(kind=wp)             :: mn
! 
!         i=1
!         mn=Nmem/Nmhm
!     
!         ! on ne prend que les TM avec m identique au m du TE
!             Do m=1,Nmhm
!         !         Do 110 n=m,(m+IDNINT(mn)-1)
!                 ! En multimode le n peut etre < a m car "a" augmente Ez/Ey#a/b.n/m
!                 Do n=1,Nint(mn)
!                     kc_e(i)=sqrt(m**2*b(p)/a+n**2*a/b(p))*PI/sqrt(a*b(p))
!                     me(1,i)=m
!                     me(2,i)=n
!                     i=i+1
!                 end do
!             end do
!         !      ri=i-1
!     
!     end function


!     ! ***
!     ! Fonction calculant la racine de l'impedance carcteristique du mode h 
!     ! ***
!     function sqrt_Zc_h(p,a,b,Nmh)
!         use aloha_config, only : k0,Nmax
!         
!         implicit none
!         
!         ! input arguments
!         integer, intent(in)       :: p, Nmh
!         real(kind=wp), intent(in) :: a, b(:)
!         ! result
!         complex(kind=wp)          :: sqrt_Zc_h(Nmh,Nmh)
!         
!         ! working variables 
!         integer                   :: i, mh(2,Nmh*Nmh)
!         real(kind=wp)             :: kch(Nmh*Nmh)
!         complex(kind=wp)          :: RZh(Nmh),RYh(Nmh),tampon, tamponbis
!         complex(kind=wp)          :: RZhi(Nmh,Nmh),RYhi(Nmh,Nmh)        
!                 
!         Call ClassModeTE(Nmax,p,a,b,mh,kch)
!             
!         Do i=1,Nmh
! !         RZh(i) = csqrt(Z0/csqrt((1.,0.)*(1-(kch(i)/k0)**2)))
!             tampon = (0.,1.)/((cl**2*Eps0)*sqrt((1.,0.)*(kch(i)**2-k0**2)))
!             tamponbis = k0*cl*tampon
!             RZh(i) = sqrt(tamponbis)
! !          RYh(i) = 1/RZh(i)
!         end do
! 
!         sqrt_Zc_h = diag(RZh)    
! 
!     end function 


    !
    !+ Generate an identity matrix of size dim x dim
    ! 
    function eye(dim)
        ! input var.
        Integer, intent(in) :: dim
        ! output var.
        Real(kind=wp)       :: eye(dim,dim)
        ! working var.
        Integer             :: i,j
        
        ! on s'assure de l'initialisation 
        ! (meme si elle est parfois intrinseque selon le compilateur...)
        eye(:,:) = 0.0
        ! 1 sur la diagonale
        forall (i=1:dim, j=1:dim, i == j) eye(i,j) = 1.0
    end function eye

    
    !
    !+ Generate a diagonal matrix from a vector v
    !
    ! The matrix diag as the size (N,N), where N is the length of the vector v
    ! and diag(i,i)=v(i)
    !
    function diag(v)
        ! input var.
        complex(kind=wp), intent(in) :: v(:)
        ! output var.
        complex(kind=wp)             :: diag(size(v),size(v))
        
        diag = eye(size(v))*spread(v, 1, size(v))        
    end function diag
        
        
    ! 
    !+ Real part of f011 
    ! 
    function rf011(t)
        implicit none
        real(kind=wp), intent(in) :: t(:)
        real(kind=wp)             :: rf011(size(t))
    
        rf011=Real(f011(t))
    end function rf011

    !
    !+ Imaginary part of if011
    !
    function if011(t)
        real(kind=wp), intent(in) :: t(:)
        real(kind=wp)             :: if011(size(t))
    
        if011=aimag(f011(t))
    end function if011
    
    !
    !+ Real part of f012
    !
    function rf012(t)
        real(kind=wp), intent(in) :: t(:)
        real(kind=wp)             :: rf012(size(t))
    
        rf012=Real(f012(t))
    end function rf012

    !
    !+ Imaginary part of f012
    !
    function if012(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: if012(size(t))
    
        if012=aimag(f012(t))
    end function if012

      function rf0131(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: rf0131(size(t))

        rf0131=Real(f0131(t))
      end function rf0131

      function if0131(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: if0131(size(t))

        if0131=aimag(f0131(t))
      end function if0131


      function rf0132(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: rf0132(size(t))

        rf0132=Real(f0132(t))
      end function rf0132

      function if0132(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: if0132(size(t))

        if0132=aimag(f0132(t))
      end function if0132

      function rf0133(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: rf0133(size(t))

        rf0133=Real(f0133(t))
      end function rf0133

      function if0133(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: if0133(size(t))

        if0133=aimag(f0133(t))
      end function if0133

      function rf0134(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: rf0134(size(t))

        rf0134=Real(f0134(t))
      end function rf0134

      function if0134(t)
        Real(kind=wp), intent(in) :: t(:)
        Real(kind=wp)             :: if0134(size(t))

        if0134=aimag(f0134(t))
      end function if0134

      !
      !+ f011
      !
      function f011(t)  
        use aloha_config, only : b,z,m,n,i,j,bcte,k0,X0,D0

        intrinsic epsilon, size
        ! input var.   
        Real(kind=wp), intent(in) :: t(:)
        ! output var.
        Complex(kind=wp)          :: f011(size(t))
        ! working var.
        integer          :: p
        Real(kind=wp)    :: v(size(t))
        Real(kind=wp)    :: pc = epsilon(8.0)
        Complex(kind=wp) :: y(size(t)), ye(size(t))
        Complex(kind=wp) :: Aie(size(t)), Aied(size(t)),Ai(size(t)), Aid(size(t)), Airy(size(t))
        Complex(kind=wp) :: mode(size(t)), coupl(size(t)), f11(size(t))
        ! airy function error flag
        Integer          :: ierc

        v=(1.-t)/t

        y = (-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
    !      y=(0.,0.)
        ye= y*exp((0.,1.)*2.*pi/3.)

        
!         write(*,*) 't=',t, ', v=', v
!         write(*,*) 'y=',y
!         write(*,*) 'ye=',ye
        
        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc)
            call airy_ai(y(p), Ai(p), Aid(p), ierc)

!             Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!              Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!              Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!               Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)   
        end do
!         write(*,*) 'pc=', pc
!         write(*,*) 'Ai=', Ai
!         write(*,*) 'Aie=', Aie
!         write(*,*) ''
!         pause
    
        Airy=(Aid/Ai)-exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) &
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        coupl=(2.,0.)
    
        f011=(-1.,0.)*Airy*mode*coupl/(t*t)
        ! MAJ JH 30/04/2009
        ! Le code suivant est une reprise du code fortran 77 :
        ! Dans le code originel, pc=1e-300.
        ! lorsque le resultat des fontions d'Airy est inferieur a cette valeur
        ! alors l'intégrande est nulle.
        ! Cela correspond a la limite lorsque y->+oo, et on a Q->0 soit f011->0.
        ! Dans la version fortran90, j'avais limite pc au plus petit reel non negligeable.
        ! Toutefois, du a l'existence du OR, cette valeur pouvait etre trop petite, et ainsi
        ! fausser le resultat final en le mettant a 0 au lieu d'une petite valeur.
        ! Le code suivant est donc laisse en commentaire, tout comme dansles fonctions suivantes.
!         do p=1,size(t)
!             if (abs(Ai(p)) < pc.OR.abs(Aie(p)) < pc) then
!                 f011(p)=(0.,0.)
!             end if
!         end do
    end function

    !
    !+ f012
    !
    function f012(t)
        use aloha_config, only : b,z,m,n,i,j,bcte,k0,X0,D0
        intrinsic epsilon, size

        ! input var.
        real(kind=wp), intent(in) :: t(:)
        ! output var.
        complex(kind=wp)          :: f012(size(t))

        ! working var.
        integer             :: p
        real(kind=wp)       :: v(size(t)), pc = epsilon(8.0)!k0,X0,D0,z(60),b(60),pi,pc!,mn
        complex(kind=wp)    :: y(size(t)), ye(size(t))
        complex(kind=wp)    :: Aie(size(t)), Aied(size(t)), Ai(size(t)), Aid(size(t)), Airy(size(t))
        complex(kind=wp)    :: mode(size(t)), coupl(size(t)), f12(size(t))
        ! airy function error flag
        Integer             :: ierc

        v=(1.-t)/t
    
        y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
        !      y=(0.,0.)
        ye=y*exp((0.,1.)*2.*pi/3.)


        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc )
            call airy_ai(y(p), Ai(p), Aid(p), ierc )
!             Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!              Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!              Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!               Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)
        end do
      
        Airy=(Aid/Ai)-exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) & 
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        !      coupl=(2.,0.)*(-(-1.)**m*exp((0.,1.)*k0*(-1+(0.,1.)*v)*b(i)))
        coupl=(2.,0.)*(-(-1.)**n*exp((0.,1.)*k0*(-1+(0.,1.)*v)*b(i)))
    
        f012=(-1.,0.)*Airy*mode*coupl/(t*t)

!          do p=1,size(t)
!             if (abs(Ai(p)) < pc.OR.abs(Aie(p)) < pc) Then
!                 f012(p)=(0.,0.)
!             endif
!         end do
    end function

    !
    !+ f0131
    !
    function f0131(t)
        use aloha_config, only : b,z,m,n,i,j,k0,X0,D0 
        intrinsic epsilon, size
    
        ! input var.
        real(kind=wp), intent(in) :: t(:)
        ! output var.
        complex(kind=wp)          :: f0131(size(t))
        ! working var.
        integer          :: p
        real(kind=wp)    :: v(size(t)), pc = epsilon(8.0)
        complex(kind=wp) :: y(size(t)), ye(size(t))
        complex(kind=wp) :: Aie(size(t)), Aied(size(t)), Ai(size(t)), Aid(size(t)), Airy(size(t))
        complex(kind=wp) :: mode(size(t)), coupl(size(t))
        ! airy function error flag
        Integer          :: ierc
    
        v=(1.-t)/t
    
        y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
    !      y=(0.,0.)
        ye=y*exp((0.,1.)*2.*pi/3.)
        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc )
            call airy_ai(y(p), Ai(p), Aid(p), ierc )
!              Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!               Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!               Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!                Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)
        end do
    
        Airy=(Aid/Ai)-exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) & 
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        coupl=exp((0.,1.)*k0*(-1+(0.,1.)*v)*abs(z(j)-z(i)))
    
        f0131=(-1.,0.)*Airy*mode*coupl/(t*t)

!         do p=1,size(t)
!             if (abs(Ai(p)) < pc.OR.abs(Aie(p)) < pc) Then
!                 f0131(p)=(0.,0.)
!             end if
!         end do
    end function

    !
    !+ f0132
    !
    function f0132(t)
        use aloha_config, only : b,z,m,n,i,j,k0,X0,D0 
        intrinsic epsilon, size
    
        ! input var.
        real(kind=wp), intent(in) :: t(:)
        ! output var.
        complex(kind=wp)          :: f0132(size(t))
        ! working var.
        integer          :: p
        real(kind=wp)    :: v(size(t)), pc = epsilon(8.0)
        complex(kind=wp) :: y(size(t)), ye(size(t))
        complex(kind=wp) :: Aie(size(t)), Aied(size(t)), Ai(size(t)), Aid(size(t)), Airy(size(t))
        complex(kind=wp) :: mode(size(t)), coupl(size(t))
        ! airy function error flag
        Integer          :: ierc

        v=(1.-t)/t
    
        y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
        !      y=(0.,0.)
        ye=y*exp((0.,1.)*2.*pi/3.)
        
        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc )
            call airy_ai(y(p), Ai(p), Aid(p), ierc )
!                 Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!                 Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!                 Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!                 Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)
        end do
        
        Airy=(Aid/Ai)-exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) & 
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        coupl=exp((0.,1.)*k0*(-1+(0.,1.)*v)*abs(z(j)-z(i)+b(j)-b(i)))
    
        f0132=(-1.,0.)*Airy*mode*coupl/(t*t)
    
!         do p=1,size(t)
!             If (abs(Ai(p)) < pc.OR.abs(Aie(p)) < pc) Then
!                 f0132(p)=(0.,0.)
!             end if
!         end do
    End function


    !
    !+ f0133
    !
    function f0133(t)
        use aloha_config, only : b,z,m,n,i,j,k0,X0,D0 
        intrinsic epsilon, size
    
        ! input var.
        real(kind=wp), intent(in) :: t(:)
        ! output var.
        complex(kind=wp)          :: f0133(size(t))
        ! working var.
        integer          :: p
        real(kind=wp)    :: v(size(t)), pc = epsilon(8.0)
        complex(kind=wp) :: y(size(t)), ye(size(t))
        complex(kind=wp) :: Aie(size(t)), Aied(size(t)), Ai(size(t)), Aid(size(t)), Airy(size(t))
        complex(kind=wp) :: mode(size(t)), coupl(size(t))
        ! airy function error flag
        Integer          :: ierc

        v=(1.-t)/t
    
        y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
    !      y=(0.,0.)
        ye=y*exp((0.,1.)*2.*pi/3.)
        
        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc )
            call airy_ai(y(p), Ai(p), Aid(p), ierc )
!               Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!                Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!                Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!                 Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)
        end do

        Airy=(Aid/Ai)-exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) & 
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        coupl=exp((0.,1.)*k0*(-1+(0.,1.)*v)*abs(z(j)-z(i)+b(j)))
    
        f0133=(-1.,0.)*Airy*mode*coupl/(t*t)
  
!         do p=1,size(t)   
!             If (abs(Ai(p)) < pc.OR.abs(Aie(p)) < pc) Then
!                 f0133(p)=(0.,0.)
!             endif
!         end do

    End function


    !
    !+ f0134
    !
    function f0134(t)
        use aloha_config, only : b,z,m,n,i,j,k0,X0,D0 
        intrinsic epsilon, size

        ! input var.
        real(kind=wp), intent(in) :: t(:)
        ! output var.
        complex(kind=wp)          :: f0134(size(t))
        ! working var.
        integer          :: p
        real(kind=wp)    :: v(size(t)), pc = epsilon(8.0)
        complex(kind=wp) :: y(size(t)), ye(size(t))
        complex(kind=wp) :: Aie(size(t)), Aied(size(t)), Ai(size(t)), Aid(size(t)), Airy(size(t))
        complex(kind=wp) :: mode(size(t)), coupl(size(t))
        ! airy function error flag
        Integer          :: ierc
!         write(*,*) ' in f0134...'
        v=(1.-t)/t 
    
        y=(-v*(v+2.*(0.,1.)))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)
    !      y=(0.,0.)
        ye=y*exp((0.,1.)*2.*pi/3.)

        do p=1,size(t)
            call airy_ai(ye(p), Aie(p), Aied(p), ierc )
            call airy_ai(y(p), Ai(p), Aid(p), ierc )
!               Aied(p)=nag_airy_ai(ye(p), deriv=.true. , scale=.false.) 
!                Aie(p)=nag_airy_ai(ye(p), deriv=.false., scale=.false.) 
!                Aid(p)=nag_airy_ai(y(p) , deriv=.true. , scale=.false.)   
!                 Ai(p)=nag_airy_ai(y(p) , deriv=.false., scale=.false.)
        end do
        
        Airy=(Aid/Ai) - exp((0.,-1.)*4.*pi/3.)*(Aied/Aie)
    
        mode=(-1.+(0.,1.)*v)**2/(-v*(v+2.*(0.,1.)))**(2./3.) & 
        /(((-1.+(0.,1.)*v)**2-(n*pi/(k0*b(i)))**2)*((-1.+(0.,1.)*v)**2-(m*pi/(k0*b(j)))**2))*exp((0.,1.)*pi/6.)*(0.,1.)
    
        coupl=exp((0.,1.)*k0*(-1+(0.,1.)*v)*abs(z(j)-z(i)-b(i)))
    
        f0134=(-1.,0.)*Airy*mode*coupl/(t*t)
        

!         do p=1,size(t)
!             If (abs(Ai(p)) < pc .OR. abs(Aie(p)) < pc) Then
!                 f0134(p)=(0.,0.) 
!             end if
!         end do


    End function


    ! ***
    ! Suroutine de CALcul de RESIDUs
    ! 
    ! INPUT
    ! 
    ! OUTPUT
    ! 
    Subroutine Calresidu(residu)

    use aloha_config, only : b, z, m, n, i, j, k0, X0, D0

    implicit none

!     Common/dims/b,z
!     Common/modmn/m,n
!     Common/guid/i,j
!     Common/onde/k0
!     Common/plasma/X0,D0

    Complex(kind=wp), intent(out) :: residu

    !Integer :: m,n!,ifail,nz,i,j
    !Real(kind=wp) :: PI,z(60),b(60),k0,X0,D0!,pc
    Complex(kind=wp) :: y, Ai, Aid
    ! airy function error flag
    Integer          :: ierc

    If ((m == n).AND.(i == j)) Then
!        write(*,*) 'k0=',k0,'n=',n,'b(i)=',b(i),'D0=',D0,'X0=',X0,'pi=',pi 
       ! TODO ? JH 10/2008: Il y'a un bug dans gfortran 4.2. Lorsqu'il evalue
       ! la soustraction d'un complexe, le signe de la partie reelle influe anormalement
       ! sur le resultat. Ainsi, on est oblige de tester si la partie reelle est nulle.
       ! Si c'est le cas, on doit alors "forcer le signe du 0", car selon le mode de 
       ! compilation (optimisation ou non), le compilateur prendra la partie reelle de 0
       ! negative ou positive...
       ! Pour obtenir le bon resultat avec gfortran, il ne faut donc pas compiler avec le 
       ! flag 'optimisation' -O[1-3]
       y=((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp((0.,-1.)*pi/3.)

!        write(*,*) '(n*pi/(k0*b(i)))**2=',(n*pi/(k0*b(i)))**2
!        write(*,*) '(-(1.,0.))**(1./3.)=',(-(1.,0.))**(1./3.)
!        write(*,*) '((0.,0.)-(1.,0.))**(1./3.)=',((0.,0.)-(1.,0.))**(1./3.)
!        write(*,*) '(+(0.,0.)-(1.,0.))**(1./3.)=',(+(0.,0.)-(1.,0.))**(1./3.) 
!        write(*,*) '(-(0.,0.)-(1.,0.))**(1./3.)=',(-(0.,0.)-(1.,0.))**(1./3.)
!        write(*,*) ''
!        write(*,*) '((n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)=',((n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)
!        write(*,*) '((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)=',((1.,0.)*(n*pi/(k0*b(i)))**2-(1.,0.))**(1./3.)


!        write(*,*) 'y=', y
       call airy_ai(y, Ai, Aid, ierc)
!        write(*,*) 'Ai=', Ai, 'Aid=', Aid

! !      NAG
!        Aid = nag_airy_ai(y, deriv=.true., scale=.false. ) ! JH 17/04/2008 F90
!        Ai  = nag_airy_ai(y, deriv=.false., scale=.false. ) ! JH 17/04/2008 F90

        If (n == 0) Then
            residu=-2.*PI*k0*b(i)*Aid/Ai/((1.,0.)*(n*PI/(k0*b(i)))**2-(1.,0.))**(2./3.)*exp((0.,1.)*PI/6.)
        Else
            residu=-PI*k0*b(i)*Aid/Ai/((1.,0.)*(n*PI/(k0*b(i)))**2-(1.,0.))**(2./3.)*exp((0.,1.)*PI/6.)
        EndIf  

    Else
        residu=(0.,0.)
    EndIf

    End subroutine Calresidu


    ! ***
    ! Suroutine de CALcul de RESIDUs en dimension y (integrale I_y)
    ! 
    ! Les poles apparaissent uniquement lorsque les indices mi et mj
    ! sont egaux. Dans un tel cas, la contribution des residus a 
    ! l'integrale I_y est donnée dans [Berio, p146-148, eq. A.III.5]
    ! 
    ! Input arguments : 
    !  - a  : largeur du guide (dir. poloidale)
    !  - mi : premier indice (dir. poloidal) m du mode (m,n) du guide p
    !  - mj : premier indice (dir. poloidal) m du mode (m,n) du guide q
    ! 
    ! Output argument :
    ! - residuy 
    ! 
    Subroutine Calresiduy(a, mi, mj, residuy)
        use aloha_config, only : k0

        Integer,       intent(in) :: mi, mj
        Real(kind=wp), intent(in) :: a
        Real(kind=wp), intent(out):: residuy

        If (mi == mj) Then
            If (mi /= 0) Then
                residuy=pi*k0*a/((mi*pi/(k0*a))**2.)
            Else
                write (*,*) 'pb mi=mj=0'
                stop
                ! residuy=PI*(k0*a)**3./3.
            EndIf
        Else ! mi /= mj
            residuy=(0.,0.)
        EndIf  
    End subroutine


      !*********
      ! Couplage au plasma (V6)
      ! 
      !call airy_ai(ye(p), Aie(p), Aied(p), ierc)
      !
      !TODO : les fonctions d'airy bi ne sont calculees que pour des valeurs reelles.
      !       il faudrait trouver une routine permettant le calcul complexe.
      !
      !*********
      function f_nz(nz)

        use aloha_config, only : b, z, m, n, i, j, k0, X0, D0, d_couche, X1, D1, pertes, d_vide

        implicit none

        Real(kind=wp)   , intent(in)  :: nz(:)
        Complex(kind=wp)              :: f_nz(size(nz))

        Real(kind=wp), parameter :: PC_limit=epsilon(8.0)
        Complex(kind=wp) :: alpha
        Complex(kind=wp), Dimension(size(nz)) :: neta
        Complex(kind=wp), Dimension(size(nz)) :: Ai, Aid, Ai_0, Aid_0, Ai_d, Aid_d
        Real(kind=wp)   , Dimension(size(nz)) :: Bi_0, Bid_0, Bi_d, Bid_d 
        Complex(kind=wp), Dimension(size(nz)) :: ya, y_L, yb, num, denom, ys, ys_sur_factor, gamma, tanh_d_vide, y0
        Integer  :: ierc, p
        Real(kind=wp) :: z_r, z_i
        Complex(kind=wp) :: u

        alpha=(1.,0.)-(0.,1.)*pertes

        neta=((1.,0.)*nz**2-alpha)**(1./3.)*(D1/X1)**(2./3.)*(X1-1.)*exp(-(0.,1.)*pi/3.)
        
        do p=1,size(nz)
            call airy_ai(neta(p), Ai(p), Aid(p), ierc)
    
            If (abs(Ai(p)) < PC_limit) Then
                f_nz = (0.,0.)
                return
            Else
                y_L(p) = -alpha*(Aid(p)/Ai(p))*(X1/D1)**(1./3.)/((1.,0.)*nz(p)**2-alpha)**(2./3.)*exp((0.,1.)*pi/6.)
            EndIf

            neta(p)=((1.,0.)*nz(p)**2-alpha)**(1./3.)*(D0/X0)**(2./3.)*(X0-1.)*exp(-(0.,1.)*pi/3.)

            call airy_ai(neta(p), Ai_0(p), Aid_0(p), ierc)


            If (abs(Ai_0(p)) < PC_limit) Then
                f_nz = (0.,0.)
                return
            Else
                ya(p) = -alpha*(Aid_0(p)/Ai_0(p))*(X0/D0)**(1./3.)/((1.,0.)*nz(p)**2-alpha)**(2./3.)*exp((0.,1.)*pi/6.)
            EndIf
    
            call airy_bi(real(neta(p)), Bi_0(p), Bid_0(p), ierc, .FALSE.)
    
            If (abs(Bi_0(p)) < PC_limit) Then
                f_nz = (0.,0.) 
                return
            Else
                yb(p) = -alpha*(Bid_0(p)/Bi_0(p))*(X0/D0)**(1./3.)/((1.,0.)*nz(p)**2-alpha)**(2./3.)*exp((0.,1.)*pi/6.)
            EndIf
    
            neta(p) = ((1.,0.)*nz(p)**2-alpha)**(1./3.)*(X0/D0)**(1./3.)*(k0*d_couche+(D0/X0)*(X0-1.))*exp(-(0.,1.)*pi/3.)
    
            call airy_ai(neta(p), Ai_d(p), Aid_d(p), ierc)
            call airy_bi(real(neta(p)), Bi_d(p), Bid_d(p), ierc, .FALSE.)
    
            If ((abs(Bid_0(p)) < PC_limit).OR.(abs(Aid_0(p)) < PC_limit)) Then
                f_nz = (0.,0.)
                return    
            Else 
                num(p) = y_L(p)*(Ai_d(p)/Ai_0(p)*yb(p)-Bi_d(p)/Bi_0(p)*ya(p))+ya(p)*yb(p)*(Bid_d(p)/Bid_0(p)-Aid_d(p)/Aid_0(p))
                denom(p) = y_L(p)*(Ai_d(p)/Ai_0(p)-Bi_d(p)/Bi_0(p))-Aid_d(p)/Aid_0(p)*ya(p)+Bid_d(p)/Bid_0(p)*yb(p)
            Endif
    
            ys(p)=num(p)/denom(p)
    
!             ! Special cases to avoid : nz=+/-1
!             If (abs(nz(p)) == 1.0)then
!                 write(*,*) 'nz=1, ys(nz)=',ys(p)
!                 ys(p) = 0.0
!             EndIf


            ! Admittance propagation inside the vacuum layer of length d_vide
            gamma(p) = k0*(nz(p)**2-alpha)**(1./2.)

            ! Hyperbolic tangent calculation : tanh(gamma*v_dide)
            ! This calculation may generate NaN, which ultimatly lead to crappy NaN integration results... 
            ! 
            ! Expression 1 (definition) --> produces NaN
            !tanh_d_vide(p) = (exp(2*gamma(p)*d_vide)-1)/(exp(2*gamma(p)*d_vide)+1)
            ! Expression 2 (built in function)  --> Complex tanh only in Fortran 2008 !
            !tanh_d_vide(p) = tanh(gamma(p)*d_vide) 
            ! Expression 3 (alternate x+iy form) --> produces NaN
            !z_r = real(gamma(p)*d_vide)
            !z_i = aimag(gamma(p)*d_vide)
            !tanh_d_vide(p) = (sinh(2*z_r) + (0.,1.)*sin(2*z_i))/(cosh(2*z_r) + cos(2*z_i))
            ! Expression 4 : limit --> works !
            !tanh_d_vide(p) = gamma(p)*d_vide
            ! Expression 5 : Taylor expansion of the tanh function (source: wikipedia) --> works also !
            u = gamma(p)*d_vide
            tanh_d_vide(p) = u - (u**3.)/3. + 2.*(u**5.)/15. - 17.*(u**7.)/315.


            y0(p) = (0.,1.)*k0/gamma(p)

            ys(p) = y0(p)*(ys(p)+y0(p)*tanh_d_vide(p))/(ys(p)*tanh_d_vide(p)+y0(p))

            ys_sur_factor(p)=ys(p)/((X0/D0)**(1./3.))

            If (isnan(real(ys(p)))) then
                write(*,*) 'NaN detected ! '
                write(*,*) 'nz=',nz(p)
                write(*,*) 'gamma=', gamma(p)
                write(*,*) 'tanh_d_vide', tanh_d_vide(p)
                write(*,*) 'y0(p)', y0(p)
                write(*,*) 'ys(p)', ys(p)
            EndIf

        endDo
       

        f_nz = nz**2.*ys_sur_factor*(1.-(-1.)**n*exp((0.,1.)*k0*nz*b(i)))* &
         1/(nz**2-(n*pi/(k0*b(i)))**2)* &
         (1.-(-1.)**m*exp((0.,-1.)*k0*nz*b(j)))* &
         1/(nz**2-(m*pi/(k0*b(j)))**2)* &
         exp((0.,1.)*k0*nz*(z(i)-z(j)))
        

      End function f_nz

      ! 
      ! rf_nz(nz)
      ! 
      ! Real part of f_nz
      ! 
      function rf_nz(nz)
        implicit none
        Real(kind=wp), intent(in) :: nz(:)
        Real(kind=wp)             :: rf_nz(size(nz))

        rf_nz = real(f_nz(nz))

      End function rf_nz

      ! 
      ! if_nz(nz)
      ! 
      ! imaginary part of f_nz
      ! 
      function if_nz(nz)
        implicit none
        Real(kind=wp), intent(in) :: nz(:)
        Real(kind=wp)             :: if_nz(size(nz))

        if_nz = aimag(f_nz(nz))

      End function if_nz

END MODULE aloha_function
