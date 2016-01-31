module aloha2d_coupling
  use aloha2D_constants, only : wp

  implicit none

    ! Ports attributes
    real, private :: a1, b1, y1, z1
    real, private :: a2, b2, y2, z2
    integer, private :: m1, n1
    integer, private :: m2, n2
    character(len=1), private  :: mode1, mode2

    real :: nz_quad ! nz index of the first integration loop

  contains




    !
    !
    ! Coupling admittance integral evaluation
    !
    !
    function eval_K(aa1, bb1, yy1, zz1, mod1, mm1, nn1, aa2, bb2, yy2, zz2, mod2, mm2, nn2) result(K)
      use aloha2d_constants

      implicit none

      real, intent(in) :: aa1, bb1, aa2, bb2, yy1, yy2, zz1, zz2
      integer, intent(in) :: mm1, nn1, mm2, nn2
      character(len=1), intent(in) :: mod1, mod2
      complex :: K

      real :: Kr=0.0,Ki=0.0

      ! Copy the input parameters into the module,
      ! in order to use these parameters from the module instead of
      ! having to pass them by arguments...
      a1 = aa1
      a2 = aa2
      b1 = bb1
      b2 = bb2
      y1 = yy1
      y2 = yy2
      z1 = zz1
      z2 = zz2
      mode1 = mod1
      mode2 = mod2
      m1 = mm1
      m2 = mm2
      n1 = nn1
      n2 = nn2

      ! Real an imaginary parts
      Call integral_real_f(Kr)
      Call integral_imag_f(Ki)
      K = one*Kr+j*Ki


    end function eval_K


    subroutine spectralFunction(ny,nz,f_ny_nz)
      use aloha2d_globalParameters, only : k0
      use aloha2d_constants, only : Y0
      use aloha2d_plasma_admittance
      use aloha2d_rectangularWaveguides, only : spectralEigenVector

      implicit none

      real, intent(in) :: ny,nz
      complex, intent(out):: f_ny_nz

      complex :: g_S_interp,g_F_interp
      complex :: Ys_yy_interp, Ys_yz_interp, Ys_zy_interp, Ys_zz_interp
      complex :: Ey1,Ez1,Hy1,Hz1,Ey2,Ez2,Hy2,Hz2

      ! calculates the modal function in the (ny,nz) space
      call spectralEigenVector(ny,nz,a1,b1,m1,n1,y1,z1,mode1,Ey1,Ez1,Hy1,Hz1)
      call spectralEigenVector(ny,nz,a2,b2,m2,n2,y2,z2,mode2,Ey2,Ez2,Hy2,Hz2)


      ! evaluate the plasma admitance tensor from interpolation
      ! TODO : correct q bug ??
      call interpolate_gS_gF(ny,nz, g_S_interp,g_F_interp)
      !g_S_interp = eval_gS_fem(ny,nz)
      !g_F_interp = eval_gF_fem(ny,nz)


      ! Part to be replaced by OLGA evaluation of the admittance
      call eval_admittanceTensor(ny,nz, &
                g_S_interp,g_F_interp, &
                Ys_yy_interp,Ys_yz_interp,Ys_zy_interp,Ys_zz_interp)


      ! Integrande
      f_ny_nz = Y0*k0**2/(4.*pi**2)*&
                (conjg(Hy1)*(Ys_yy_interp*Ey2 + Ys_yz_interp*Ez2)+&
                 conjg(Hz1)*(Ys_zy_interp*Ey2 + Ys_zz_interp*Ez2))


    end subroutine spectralFunction



Subroutine integral_real_f(Kr)
      use aloha2d_constants, only : epsrel, epsabs
      use aloha2d_plasma

      implicit none

      real, intent(out) :: Kr

      integer :: npts2=2,neval,ier!,limit=400,leniw=802,lenw=2000,last,iwork(802)
      real :: a_nz,b_nz, K_tamp=0.0
      real :: points(2), abserr !,work(2000)
      real :: demi_intervalle=0.2

      real, external :: f_real_nz

      ! TODO
      ! This section could be easily distributed, each all the integral being calculated by one core

      !write(*,*) 'real part'
      !write(*,*) 'nz domain 1/5'

      a_nz=nz_min
      b_nz=-SQRT(ABS(S))-demi_intervalle
      !Call qag(f_real_nz,a_nz,b_nz,epsabs,epsrel,quad_key,K_tamp, &
      !      abserr,neval,ier,limit,lenw,last,iwork,work)
      Call qags(f_real_nz,a_nz,b_nz,epsabs,epsrel,&
                K_tamp, abserr,neval,ier)

      Kr=K_tamp

      !write(*,*) Kr
      !write(*,*) 'nz domain 2/5'

      a_nz=-SQRT(ABS(S))-demi_intervalle
      b_nz=-SQRT(ABS(S))+demi_intervalle
      points(1)=-SQRT(ABS(S))
      Call qagp(f_real_nz,a_nz,b_nz,npts2,points,epsabs,epsrel,&
          K_tamp ,abserr,neval,ier)
      Kr=Kr+K_tamp

      !write(*,*) Kr
      !write(*,*) 'nz domain 3/5'

      a_nz=-SQRT(ABS(S))+demi_intervalle
      b_nz=SQRT(ABS(S))-demi_intervalle
      Call qags(f_real_nz,a_nz,b_nz,epsabs,epsrel,&
            K_tamp,abserr, neval,ier)
      Kr=Kr+K_tamp

      !write(*,*) Kr
      !write(*,*) 'nz domain 4/5'

      a_nz=+SQRT(ABS(S))-demi_intervalle
      b_nz=+SQRT(ABS(S))+demi_intervalle
      points(1)=SQRT(ABS(S))
      Call qagp(f_real_nz,a_nz,b_nz,npts2,points,epsabs,epsrel,&
                K_tamp,abserr,neval,ier)
      Kr=Kr+K_tamp

      !write(*,*) Kr
      !write(*,*) 'nz domain 5/5'

      a_nz=+SQRT(ABS(S))+demi_intervalle
      b_nz=nz_max
      !Call qag(f_real_nz,a_nz,b_nz,epsabs,epsrel,quad_key,K_tamp,abserr, &
      !       neval,ier,limit,lenw,last,iwork,work)

      Call qags(f_real_nz,a_nz,b_nz,epsabs,epsrel,&
              K_tamp,abserr, neval,ier)
      Kr=Kr+K_tamp

      !write(*,*) Kr

End subroutine integral_real_f


Subroutine integral_imag_f(Ki)
      use aloha2d_constants, only : quad_key, epsrel, epsabs
      use aloha2d_plasma

      implicit none

      real, intent(out) :: Ki

      integer :: npts2=2,neval,ier!,limit=400,leniw=802,lenw=2000,last,iwork(802)
      real :: a_nz,b_nz, Ki_tamp=0.0
      real :: points(2),abserr !,work(2000),
      real :: demi_intervalle=0.2

      real, external :: f_imag_nz

      ! TODO
      ! This section could be easily distributed, each all the integral being calculated by one core

      !write(*,*) 'imag part'
      !write(*,*) 'nz domain 1/5'

      a_nz=nz_min
      b_nz=-SQRT(ABS(S))-demi_intervalle
      Call qags(f_imag_nz,a_nz,b_nz,epsabs,epsrel,&
                Ki_tamp, abserr,neval,ier)

      Ki=Ki_tamp

      !write(*,*) Ki
      !write(*,*) 'nz domain 2/5'

      a_nz=-SQRT(ABS(S))-demi_intervalle
      b_nz=-SQRT(ABS(S))+demi_intervalle
      points(1)=-SQRT(ABS(S))
      Call qagp(f_imag_nz,a_nz,b_nz,npts2,points,epsabs,epsrel,&
          Ki_tamp ,abserr,neval,ier)
      Ki=Ki+Ki_tamp

      !write(*,*) Ki
      !write(*,*) 'nz domain 3/5'

      a_nz=-SQRT(ABS(S))+demi_intervalle
      b_nz=SQRT(ABS(S))-demi_intervalle
      Call qags(f_imag_nz,a_nz,b_nz,epsabs,epsrel,&
            Ki_tamp,abserr, neval,ier)
      Ki=Ki+Ki_tamp

      !write(*,*) Ki
      !write(*,*) 'nz domain 4/5'

      a_nz=+SQRT(ABS(S))-demi_intervalle
      b_nz=+SQRT(ABS(S))+demi_intervalle
      points(1)=SQRT(ABS(S))
      Call qagp(f_imag_nz,a_nz,b_nz,npts2,points,epsabs,epsrel,&
                Ki_tamp,abserr,neval,ier)
      Ki=Ki+Ki_tamp

      !write(*,*) Ki
      !write(*,*) 'nz domain 5/5'

      a_nz=+SQRT(ABS(S))+demi_intervalle
      b_nz=nz_max
      Call qags(f_imag_nz,a_nz,b_nz,epsabs,epsrel,&
              Ki_tamp,abserr, neval,ier)
      Ki=Ki+Ki_tamp

      !write(*,*) Ki

End subroutine integral_imag_f


! Calculates arrays of spectral electric and magnetic field components
! Returned array of size nb_ny x nb_nz
subroutine spect_discr(aa, bb, yy, zz, mode, mm, nn, &
                       eyt_ny_nz, ezt_ny_nz, hyt_ny_nz, hzt_ny_nz)
    use aloha2d_plasma_admittance
    use aloha2d_rectangularwaveguides

    implicit none

    real, intent(in) :: aa, bb, yy, zz
    character(len=1), intent(in) :: mode
    integer, intent(in) :: mm, nn
    complex, dimension(:), allocatable, intent(out) :: eyt_ny_nz, ezt_ny_nz, hyt_ny_nz, hzt_ny_nz

    integer :: incr ! arrays index
    real :: ny, nz
    complex :: ey_ny_nz, ez_ny_nz, hy_ny_nz, hz_ny_nz ! spectral components of the field

    allocate(eyt_ny_nz(GRID_NY_NB*GRID_NZ_NB), ezt_ny_nz(GRID_NY_NB*GRID_NZ_NB))
    allocate(hyt_ny_nz(GRID_NY_NB*GRID_NZ_NB), hzt_ny_nz(GRID_NY_NB*GRID_NZ_NB))

    incr=1

    ! nz first then ny
    Do n1=0,GRID_NZ_NB-1
      !write(*,*) n1+1,'/',GRID_NZ_NB
      nz = GRID_NZ_MIN + n1*GRID_DNZ

      Do n2=0,GRID_NY_NB-1
            ny = GRID_NY_MIN + n2*GRID_DNY

            call spectralEigenVector(ny, nz, aa, bb, mm, nn, yy, zz, mode, &
                                    ey_ny_nz, ez_ny_nz, hy_ny_nz, hz_ny_nz)

            ! save into the arrays
            eyt_ny_nz(incr)=ey_ny_nz
            ezt_ny_nz(incr)=ez_ny_nz
            hyt_ny_nz(incr)=hy_ny_nz
            hzt_ny_nz(incr)=hz_ny_nz

            incr=incr+1
        end do
    end do
end subroutine spect_discr



end module aloha2d_coupling


      !
      ! integrand of the spectral function
      !
      !
      function f_real(ny) result(f_r)
        use aloha2d_constants, only : wp
        use aloha2d_coupling, only : nz_quad, spectralFunction

        implicit none

        real, intent(in) :: ny ! quadpack uses default precision real
        real :: f_r

        complex :: f_ny_nz

        call spectralFunction(ny,nz_quad,f_ny_nz)

        f_r = real(f_ny_nz)

      End function f_real

      !
      ! integrand of the spectral function
      !
      !
      function f_imag(ny) result(f_i)
        use aloha2d_constants, only : wp
        use aloha2d_coupling, only : nz_quad, spectralFunction

        implicit none

        real, intent(in) :: ny ! quadpack uses default precision real
        real :: f_i

        complex :: f_ny_nz

        call spectralFunction(ny,nz_quad,f_ny_nz)

        f_i = aimag(f_ny_nz)

      End function f_imag


  !
  ! Integrand F(ny,nz), integrated along ny
  !
  function f_real_nz(nz_dp) result(f_nz_r)
    use aloha2d_globalParameters, only : ny_min, ny_max !
    use aloha2d_constants ! epsrel, epsabs, quadkey
    use aloha2d_plasma
    use aloha2d_coupling, only : nz_quad

    implicit none

    real, intent(in) :: nz_dp ! quadpack use default real precision
    real :: f_nz_r

    real :: nz

    real :: tamp, demi_intervalle=0.2
    real :: a_ny, b_ny, points(2), abserr
    integer :: npts2=1+2 ! one break point +2 (cf doc)
    integer :: ier,neval!,limit=400,leniw=802,lenw=2000,last,work,iwork(802)

    real, external :: f_real

      !write(*,*) '    f_real_nz'
    ! copy the current nz value to global memory
    ! casting the values to default real for quadpack
    nz = nz_dp
    nz_quad = nz

      If (SQRT(ABS(S-nz**2.)).LE.ny_max) Then

         If (SQRT(ABS(S-nz**2.)).GE.(ny_max-demi_intervalle)) Then
            !write(*,*) '    1st ny domain'

            a_ny=real(ny_min) ! cast to default real precision for compatibility with quadpack
            b_ny=real(-SQRT(ABS(S-nz**2.))+demi_intervalle)
            points(1)=-SQRT(ABS(S-nz**2.))
            Call qagp(f_real,a_ny,b_ny,npts2,points,epsabs, &
              epsrel,tamp,abserr,neval,ier)
            f_nz_r=tamp
            !write(*,*) '    2nd ny domain'
            a_ny=real(-SQRT(ABS(S-nz**2.))+demi_intervalle)
            b_ny=real(SQRT(ABS(S-nz**2.))-demi_intervalle)
            Call qags(f_real,a_ny,b_ny,epsabs,epsrel,tamp,abserr, &
                     neval,ier)
            f_nz_r=f_nz_r+tamp
            !write(*,*) '    3rd and last ny domain'
            a_ny=real(+SQRT(ABS(S-nz**2.))-demi_intervalle)
            b_ny=real(ny_max)
            points(1)=SQRT(ABS(S-nz**2.))
            Call qagp(f_real,a_ny,b_ny,npts2,points,epsabs, &
              epsrel,tamp,abserr,neval,ier)
            f_nz_r=f_nz_r+tamp

         Else
            !write(*,*) '    1st ny domain'
            a_ny=real(ny_min)
            b_ny=real(-SQRT(ABS(S-nz**2.))-demi_intervalle)
            Call qags(f_real,a_ny,b_ny,epsabs,epsrel,tamp,abserr, &
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_r=tamp
            !write(*,*) '    2nd ny domain'
            a_ny=real(-SQRT(ABS(S-nz**2.))-demi_intervalle)
            b_ny=real(-SQRT(ABS(S-nz**2.))+demi_intervalle)
            points(1)=-SQRT(ABS(S-nz**2.))
            Call qagp(f_real,a_ny,b_ny,npts2,points,epsabs, &
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_r=f_nz_r+tamp
            !write(*,*) '    3rd ny domain'
            a_ny=real(-SQRT(ABS(S-nz**2.))+demi_intervalle)
            b_ny=real(SQRT(ABS(S-nz**2.))-demi_intervalle)
            Call qags(f_real,a_ny,b_ny,epsabs,epsrel,tamp,abserr, &
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_r=f_nz_r+tamp
            !write(*,*) '    4th ny domain'
            a_ny=real(+SQRT(ABS(S-nz**2.))-demi_intervalle)
            b_ny=real(+SQRT(ABS(S-nz**2.))+demi_intervalle)
            points(1)=SQRT(ABS(S-nz**2.))
            Call qagp(f_real,a_ny,b_ny,npts2,points,epsabs, &
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_r=f_nz_r+tamp
            !write(*,*) '    5th and last ny domain'
            a_ny=real(SQRT(ABS(S-nz**2.))+demi_intervalle)
            b_ny=real(ny_max)
            Call qags(f_real,a_ny,b_ny,epsabs,epsrel,tamp,abserr, &
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_r=f_nz_r+tamp

         Endif

      Else
        !write(*,*) '    1st and last ny domain'
         a_ny=real(ny_min)
         b_ny=real(ny_max)
         Call qags(f_real,a_ny,b_ny,epsabs,epsrel,f_nz_r,abserr, &
             neval,ier)

      Endif
  End function f_real_nz

  !
  ! Integrand F(ny,nz), integrated along ny
  !
  function f_imag_nz(nz) result(f_nz_i)
    use aloha2d_globalParameters
    use aloha2d_constants ! epsrel, epsabs, quadkey
    use aloha2d_plasma
    use aloha2d_coupling, only : nz_quad

    implicit none

    real, intent(in) :: nz
    real :: f_nz_i

    real :: tamp, demi_intervalle=0.2
    real :: a_ny, b_ny, points(2), abserr
    integer :: npts2=2
    integer :: ier,neval!,limit=400,leniw=802,lenw=2000,last,work,iwork(802)

    real, external :: f_imag

    ! copy the current nz value to global memory
    nz_quad = nz

      If (SQRT(ABS(S-nz**2.)).LE.ny_max) Then

         If (SQRT(ABS(S-nz**2.)).GE.(ny_max-demi_intervalle)) Then

            a_ny=ny_min
            !TODO : est ce qu il n'y a pas un probleme avec les bornes ci-dessous ??
            b_ny=-SQRT(ABS(S-nz**2.))+demi_intervalle
            points(1)=-SQRT(ABS(S-nz**2.))
            Call qagp(f_imag,a_ny,b_ny,npts2,points,epsabs,&
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=tamp

            a_ny=-SQRT(ABS(S-nz**2.))+demi_intervalle
            b_ny=SQRT(ABS(S-nz**2.))-demi_intervalle
            Call qags(f_imag,a_ny,b_ny,epsabs,epsrel,tamp,abserr,&
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=f_nz_i+tamp

            a_ny=+SQRT(ABS(S-nz**2.))-demi_intervalle
            b_ny=ny_max
            points(1)=SQRT(ABS(S-nz**2.))
            Call qagp(f_imag,a_ny,b_ny,npts2,points,epsabs,&
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=f_nz_i+tamp

          Else
            a_ny=real(ny_min)
            b_ny=real(-SQRT(ABS(S-nz**2.))-demi_intervalle)
            Call qags(f_imag,a_ny,b_ny,epsabs,epsrel,tamp,abserr,&
                     neval,ier)
            ! JH 24/03/2013
            ! The follwing isNaN test below intend to avoid
            ! a NaN value in tamp for some high order mode (TE0n?)
            if (isNaN(tamp)) tamp = 0.0;


            f_nz_i=tamp

            a_ny=real(-SQRT(ABS(S-nz**2.))-demi_intervalle)
            b_ny=real(-SQRT(ABS(S-nz**2.))+demi_intervalle)
            points(1)=-SQRT(ABS(S-nz**2.))
            Call qagp(f_imag,a_ny,b_ny,npts2,points,epsabs,&
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=f_nz_i+tamp

            a_ny=-SQRT(ABS(S-nz**2.))+demi_intervalle
            b_ny=SQRT(ABS(S-nz**2.))-demi_intervalle
            Call qags(f_imag,a_ny,b_ny,epsabs,epsrel,tamp,abserr,&
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;

            f_nz_i=f_nz_i+tamp

            a_ny=+SQRT(ABS(S-nz**2.))-demi_intervalle
            b_ny=+SQRT(ABS(S-nz**2.))+demi_intervalle
            points(1)=SQRT(ABS(S-nz**2.))
            Call qagp(f_imag,a_ny,b_ny,npts2,points,epsabs,&
              epsrel,tamp,abserr,neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=f_nz_i+tamp

            a_ny=SQRT(ABS(S-nz**2.))+demi_intervalle
            b_ny=ny_max
            Call qags(f_imag,a_ny,b_ny,epsabs,epsrel,tamp,abserr,&
                     neval,ier)
            if (isNaN(tamp)) tamp = 0.0;
            f_nz_i=f_nz_i+tamp

         Endif

      Else ! Pas de pole dans le domaine d' integration : un seul interval d' integration

         a_ny=ny_min
         b_ny=ny_max
         Call qags(f_imag,a_ny,b_ny,epsabs,epsrel,f_nz_i,abserr,&
             neval,ier)

      Endif

End function f_imag_nz
