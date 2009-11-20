module aloha_integration

Use aloha_constants

contains
    ! INTegration de 0 a 1 routine NAG D01AHF
    ! 
    ! INPUT
    !  - epsrel
    !  - i
    !  - j
    ! 
    ! OUTPUT
    !  - resultr
    !  - resulti
    !  - absrm
    ! 
    Subroutine Int01AHF(epsrel,i,j,resultr,resulti,absrm)
    
        Use aloha_config, Only : nlimit, b, z, m, n, max_nz
        Use aloha_function

        ! Numerical integration quadpack
        ! http://orion.math.iastate.edu/burkardt/f_src/quadpack/quadpack.f90
        Use quadpack, Only : qag

!         ! Numerical integration NAG
!         ! http://www.nag.co.uk/numeric/fn/manual/pdf/c11/c11m01_quad_1d_fn03.pdf
!         Use nag_quad_1d,  Only : nag_quad_1d_gen    
!         Use nag_error_handling

        implicit none

        intrinsic epsilon

        Real(kind=wp), intent(in)  :: epsrel
        Integer      , intent(in)  :: i, j
        Real(kind=wp), intent(out) :: resultr, resulti, absrm


        Integer       :: ct, ifail, ni, k
        Real          :: mn
        Real(kind=wp) :: d, f, relr
        Real(kind=wp), dimension(12) :: absr = 0.0
        Real(kind=wp) :: resultr0112, resultr013
        Real(kind=wp) :: resulti0112, resulti013
        Real(kind=wp) :: resultr011, resultr012
        Real(kind=wp) :: resulti011, resulti012
        Real(kind=wp) :: resultr0131, resultr0132, resultr0133, resultr0134
        Real(kind=wp) :: resulti0131, resulti0132, resulti0133, resulti0134
        Complex(kind=wp) :: residu

        ! quadpack, qag
        ! JH 10/2008 : les constantes qui suivent influent evidemment sur la precision
        ! des integrations, et forcement aussi enormement sur le temps de calcul.
        ! Baisser la precision absolue par rapport a la valeur par defaut NAG amene un gain
        ! en temps x7 environ sur le total de toutes les integrations.
        real(kind=wp) :: eps_rel=1.0e-4 ! relative accuracy required. ex: 1.0e-4 (defaut value for NAG)
        real(kind=wp) :: eps_abs=1.0e-3 ! absolute accuracy required. ex: sqrt(epsilon(1.0_wp)) (defaut value for NAG) 
        integer :: neval(12), ier(12)
        integer :: key=1

!         ! NAG
!         type(nag_error) :: error
!         ! ! On remonte le niveau d'arret sur erreur de 2 (defaut) à 3
!         ! ! Si on gardait 2, l'integration par la lib nag provoquerait une erreur
!         ! ! due à la non-convergence des intégrales.
!         ! ! On supprime également l'affichage des warnings si  print_level=3
!         Call nag_set_error(error, halt_level=3, print_level=3)
    
        d=1./(max_nz+1.)
        f=1.0 !- epsilon(1.0) ! JH/DV : epsilon pour eviter f=1 qui provoque une singularite 
    !      nlimit=0
        ifail=-1

        if (i == j) then
            mn=(m+n)/2.
            if (mn == int(mn)) then
                ! Integration quadpack
                call qag(rf011, d, f, eps_abs, eps_rel, key, resultr011, absr(1), neval(1), ier(1))
                call qag(if011, d, f, eps_abs, eps_rel, key, resulti011, absr(2), neval(2), ier(2))
                call qag(rf012, d, f, eps_abs, eps_rel, key, resultr012, absr(3), neval(3), ier(3))
                call qag(if012, d, f, eps_abs, eps_rel, key, resulti012, absr(4), neval(4), ier(4))

!                 ! Integration NAG d01alf_9
!                 call nag_quad_1d_gen(rf011, d, f, resultr011, abs_err=absr(1))
!                 call nag_quad_1d_gen(if011, d, f, resulti011, abs_err=absr(2))
!                 call nag_quad_1d_gen(rf012, d, f, resultr012, abs_err=absr(3))
!                 call nag_quad_1d_gen(if012, d, f, resulti012, abs_err=absr(4))

                resultr0112=resultr011+resultr012
                resulti0112=resulti011+resulti012
            Else
                resultr0112=(0.,0.)
                resulti0112=(0.,0.)
            EndIf
            resultr013=(0.,0.)
            resulti013=(0.,0.)
            absr(5:12)=0.
        Else ! i/= j
            resultr0112=(0.,0.)
            resulti0112=(0.,0.)
            absr(1:4)=0.

                ! Integration quadpack
            call qag(rf0131, d, f, eps_abs, eps_rel, key, resultr0131, absr(5), neval(5), ier(5))
            call qag(if0131, d, f, eps_abs, eps_rel, key, resulti0131, absr(6), neval(6), ier(6))
            call qag(rf0132, d, f, eps_abs, eps_rel, key, resultr0132, absr(7), neval(7), ier(7))
            call qag(if0132, d, f, eps_abs, eps_rel, key, resulti0132, absr(8), neval(8), ier(8))
            call qag(rf0133, d, f, eps_abs, eps_rel, key, resultr0133, absr(9), neval(9), ier(9))
            call qag(if0133, d, f, eps_abs, eps_rel, key, resulti0133, absr(10), neval(10), ier(10))
            call qag(rf0134, d, f, eps_abs, eps_rel, key, resultr0134, absr(11), neval(11), ier(11))
            call qag(if0134, d, f, eps_abs, eps_rel, key, resulti0134, absr(12), neval(12), ier(12))

!             ! Integration NAG d01alf_9d
!             call nag_quad_1d_gen(rf0131, d, f, resultr0131, abs_err=absr(5), error=error)
!             call nag_quad_1d_gen(if0131, d, f, resulti0131, abs_err=absr(6), error=error)
!             call nag_quad_1d_gen(rf0132, d, f, resultr0132, abs_err=absr(7), error=error)
!             call nag_quad_1d_gen(if0132, d, f, resulti0132, abs_err=absr(8), error=error)
!             call nag_quad_1d_gen(rf0133, d, f, resultr0133, abs_err=absr(9), error=error)
!             call nag_quad_1d_gen(if0133, d, f, resulti0133, abs_err=absr(10), error=error)
!             call nag_quad_1d_gen(rf0134, d, f, resultr0134, abs_err=absr(11), error=error)
!             call nag_quad_1d_gen(if0134, d, f, resulti0134, abs_err=absr(12), error=error)

            resultr013=resultr0131+(-1.)**(m+n)*resultr0132-(-1.)**m*resultr0133-(-1.)**n*resultr0134
            resulti013=resulti0131+(-1.)**(m+n)*resulti0132-(-1.)**m*resulti0133-(-1.)**n*resulti0134

    
        EndIf
    
        

!         write(*,*) 'resultr0112=',resultr0112, 'resulti0112=', resulti0112
!         write(*,*) 'resultr013=', resultr013, 'resulti013=',resulti013

        If ((m == n).AND.(i == j)) Then
            Call Calresidu(residu)
        Else
            residu=(0.,0.)
        EndIf
    
!         write(*,*) 'residu=',residu
!         write(*,*) 'resultr0112=',resultr0112, 'resulti0112=',resulti0112
!         write(*,*) 'resultr013=',resultr013, 'resulti013=',resulti013

        If (i == j) Then
            resultr=real(residu)-resultr0112
            resulti=aimag(residu)-resulti0112
        Else
            resultr=real(residu)-resultr013
            resulti=aimag(residu)-resulti013
        EndIf
    
!         write(*,*) 'resultr=',resultr,'resulti=',resulti

        absrm=minval(absr)
!         do ct=2,12
!             If (absrm < absr(ct)) Then
!                 absrm=absr(ct)
!             EndIf
!         end do
    End subroutine Int01AHF


! ******************************************************
!  subroutine d'integration de -1000 a 1000
    Subroutine Int01AJF(epsrel,i,j,resultr,resulti,absrm)

        Use aloha_config, Only : nlimit, b, z, m, n, max_nz
        Use aloha_function, Only : rf_nz, if_nz 

        Use quadpack, Only : qag

        implicit none
        intrinsic epsilon

!     Common/err/epsabs
!     Common/integration/max_nz      
!     Common/dim/b,z
!     Common/mod/m,n
!     Common/guid/i,j


        Integer :: LW, LIW, i, j
        Integer, Dimension(5000) :: IW

        Real(kind=wp), Dimension(20000) :: W
        Real(kind=wp) :: epsrel
        Real(kind=wp) :: min_, max_, resultr, resulti
        Real(kind=wp), Dimension(2) :: absr = 0
        Real(kind=wp) :: absrm

        ! quadpack, qag
        real(kind=wp) :: eps_rel=1.0e-4 
        real(kind=wp) :: eps_abs=sqrt(epsilon(1.0_wp)) 
        integer, Dimension(2) :: neval, ier
        integer :: key=1  

    If ((b(i)-b(j) == 0).AND.((z(i)-z(j)) == 0.).AND.(MOD((m+n),2) == 1)) Then
        resultr = 0.0
        resulti = 0.0
    Else 
        LW = 20000
        LIW = LW/4
        min_ = -max_nz
        max_ = max_nz

        call qag(rf_nz, min_, max_, eps_abs, eps_rel, key, resultr, absr(1), neval(1), ier(1))
        call qag(if_nz, min_, max_, eps_abs, eps_rel, key, resulti, absr(2), neval(2), ier(2))
    Endif
    End subroutine Int01AJF

    ! Riemann sum (from the left)
    ! Simplest approximation of the integral of f(t) for t in [a,b].
    ! The integration is performed on N points.
    !
    ! Input :
    !  - f : name of the function to perform
    !  - a : lower integration point
    !  - b : upper integration point
    !  - N : number of integration points
    !
    ! Output : s
    !
    subroutine riemann_sum(func, a, b, N, s)
        interface
            function func(x)
                use aloha_constants, only : wp
                intrinsic size
                real(kind=wp), intent(in) :: x(:)
                real(kind=wp)             :: func(size(x))
            end function func
        end interface

        integer, intent(in)         :: N
        real(kind=wp), intent(in)   :: a, b
        real(kind=wp), intent(out)  :: s

        integer                     :: i
        real(kind=wp)               :: delta, x(N)

        delta = (b-a)/N;

        do i=0,N-1
            x(i) = a + i*delta
        end do

        s = delta*sum(func(x))
    end subroutine riemann_sum

    ! polynom interpolation
    SUBROUTINE polint(xa,ya,n,x,y,dy)
        integer       :: n
        real(kind=wp) :: dy,x,y,xa(n),ya(n)
        !         Largest anticipated value of n.            
        integer, parameter :: NMAX=10 
        !         Given arrays xa and ya, each of length n, and given a value x, this routine returns a
        !         value y, and an error estimate dy. If P(x) is the polynomial of degree N ? 1 such that
        !         P(xai) = yai, i = 1, . . . , n, then the returned value y = P(x).
        INTEGER       :: i,m,ns
        real(kind=wp) :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
        ns=1
        dif=abs(x-xa(1))
        !         Here we find the index ns of the closest table entry,
        do i=1,n 
            dift=abs(x-xa(i))
            if (dift.lt.dif) then
                ns=i
                dif=dift
            endif
        !         and initialize the tableau of c?s and d?s.
            c(i)=ya(i) 
            d(i)=ya(i)
        enddo
        !         This is the initial approximation to y.
        y=ya(ns) 
        ns=ns-1
        !         For each column of the tableau,
        do m=1,n-1 
        !         we loop over the current c?s and d?s and update them.
            do i=1,n-m 
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                if(den.eq.0.) write (*,*) 'failure in polint'
        !         This error can occur only if two input xa?s are (to within roundoff) identical.
                den=w/den
        !         Here the c?s and d?s are updated.
                d(i)=hp*den 
                c(i)=ho*den
            enddo
            if (2*ns.lt.n-m) then 
            !         After each column in the tableau is completed, we decide
            !         which correction, ! or d, we want to add to our accumulating
            !         value of y, i.e., which path to take through
            !         the tableau?forking up or down. We do this in such a
            !         way as to take the most ?straight line? route through the
            !         tableau to its apex, updating ns accordingly to keep track
            !         of where we are. This route keeps the partial approximations
            !         centered (insofar as possible) on the target x. T he
            !         last dy added is thus the error indication.
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            endif
            y=y+dy
        enddo
        return
    end subroutine polint


    ! integration par la methode de Romberg        
    Subroutine qromb(func,a,b,ss)
        interface
            function func(x)
                use aloha_constants, only : wp
                intrinsic size
                real(kind=wp), intent(in) :: x(:)
                real(kind=wp)             :: func(size(x))
            end function func
        end interface

        real(kind=wp) :: a,b,ss

        integer      , parameter :: JMAX=20
        integer      , parameter :: JMAXP=JMAX+1
        integer      , parameter :: K=5
        integer      , parameter :: KM=K-1
        real(kind=wp), parameter :: EPS=1.e-6
        ! USES polint,trapzd
        ! Returns as ss the integral of the function fun! from a to b. Integration is performed by
        ! Romberg's method of order 2K, where, e.g., K=2 is Simpson?s rule.
        ! Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation
        ! error estimate; JMAX limits the total number of steps; K is the number of points used in
        ! the extrapolation.
        integer       :: j
        ! These store the successive trapezoidal approximations
        real(kind=wp) ::  dss,h(JMAXP),s(JMAXP) 
        ! and their relative stepsizes.
        h(1)=1. 
        do j=1,JMAX
            call trapzd(func,a,b,s(j),j)
            if (j.ge.K) then
                call polint(h(j-KM),s(j-KM),K,0.0_wp,ss,dss)
                if (abs(dss).le.EPS*abs(ss)) return
            endif
            s(j+1)=s(j)

        ! This is a key step: The factor is 0.25 even though
        ! the stepsize is decreased by only 0.5. This makes
        ! the extrapolation a polynomial in h2 as allowed
        ! by equation (4.2.1), not just a polynomial in h.

            h(j+1)=0.25*h(j) 
        enddo
        !  write (*,*) 'too many steps in qromb'
    end subroutine qromb

    ! subroutine d'integration classique (trapeze)
    ! This routine computes the nth stage of refinement of an extended trapezoidal rule. fun is
    ! input as the name of the function to be integrated between limits a and b, also input. When
    ! called with n=1, the routine returns as s the crudest estimate of integral from a to 0 of f(x)dx. 
    ! Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2^(n-2)
    ! additional interior points. s should not be modified between sequential calls
    Subroutine trapzd(func,a,b,s,n)
        interface
            function func(x)
                use aloha_constants, only : wp
                intrinsic size
                real(kind=wp), intent(in) :: x(:)
                real(kind=wp)             :: func(size(x))
            end function func
        end interface
        
        integer       :: n
        real(kind=wp) :: a, b, s(1)
    !     EXTERNAL func
    
        integer       :: it, j
        real(kind=wp) :: del, som(1), tnm, x
        
        if (n == 1) then
            s=0.5*(b-a)*(func((/a/))+func((/b/)))
        else
            it=2**(n-2)
            tnm=it
            ! This is the spacing of the points to be added
            del=(b-a)/tnm 
            x=a+0.5*del
            som=0.
            do j=1,it
                som=som+func((/x/))
                x=x+del
            enddo
            ! This replaces s by its refined value.
            s=0.5*(s+(b-a)*som/tnm) 
        endif
        return
    end subroutine trapzd

    ! *********************************************************
    ! Returns as s the integral of the function fun! from a to b. 
    !  Integration is performed by the trapezoidal rule.
    Subroutine qtrap(func,a,b,s)
        interface
            function func(x)
                use aloha_constants, only : wp
                intrinsic size
                real(kind=wp), intent(in) :: x(:)
                real(kind=wp)             :: func(size(x))
            end function func
        end interface
        intrinsic epsilon 
        ! On utilise un tableau s_tab(1) pour être copatible avec les fonctions
        ! func qui sont generiques (ie : admette en entree un tableau et non un scalaire)
        ! [ car en fortran, un scalaire ne correspond pas à un tableau(1) ! ]
        real(kind=wp) :: a, b, s_tab(1), s
        ! The parameters EPS can be set to the desired fractional accuracy 
        ! and JMAX so that 2 to the power JMAX-1 is the maximum allowed number of steps.
        ! MODIF JH : JMAX=12 au lieu de 20 (pour accelerer les calculs au detriment de la convergence)
        integer      , parameter :: JMAX=12
        real(kind=wp), parameter :: EPS=epsilon(8.0)
    
        integer       :: j
        real(kind=wp) :: olds(1)
        ! Any number that is unlikely to be the average of the function
        ! at its endpoints will do here.
        olds=-1.e30 
        do j=1,JMAX 
            call trapzd(func,a,b,s_tab,j)
            ! Avoid spurious early convergence.
            if (j > 5) then 
                if (abs(s_tab(1)-olds(1)) < (EPS*abs(olds(1))).or.((s_tab(1) == 0).and.(olds(1) == 0))) return
            endif
            olds=s_tab(1)
        enddo
        s=s_tab(1)
        ! pause 'too many steps in qtrap'
        ! write (*,*) '[S_grill_version6] (WARNING) : too many steps in qtrap'
    end subroutine qtrap


end module aloha_integration


