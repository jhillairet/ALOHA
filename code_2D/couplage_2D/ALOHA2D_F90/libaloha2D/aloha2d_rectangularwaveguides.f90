! Rectangular waveguide module
! This module contains all the specific calculations
! related to rectangular waveguide
module aloha2d_rectangularwaveguides
  use aloha2d_constants

  implicit none

  contains

  !
  ! Cut-off wavenumber kc of a rectangular waveguide
  !
  function rectwg_kc(a,b,m,n) result(kc)
      real, intent(in)     :: a, b
      integer,intent(in)   :: m, n
      real                 :: kc

      kc = sqrt((m*pi/a)**2 + (n*pi/b)**2);

  end function rectwg_kc

  !
  ! Guided wavenumber beta (kg) of a rectangular waveguide
  !
  function rectwg_beta(f,a,b,m,n) result(beta)
      real, intent(in)     :: f, a, b
      integer,intent(in)   :: m, n
      complex              :: beta

      real :: kc, k0

      ! vacuum wavenumber (evaluated here, in order to avoid using a global parameter)
      k0 = 2*pi*f/c0
      ! Cutoff wavenumber
      kc = rectwg_kc(a,b,m,n)
      ! guided wavenumber
      ! TODO : shall we treat differently the case of complex valued (negative sqrt) ? Or do fortran do the thing correctly ?
      beta = sqrt(one*(k0**2 - kc**2))

  end function rectwg_beta

  !
  ! Complex propagation constant gamma=j*kg
  !
  function rectwg_gamma(f,a,b,m,n)
      real, intent(in)     :: f, a, b
      integer,intent(in)   :: m, n
      complex              :: rectwg_gamma

      real :: kc, k0

      ! vacuum wavenumber (evaluated here, in order to avoid using a global parameter)
      k0 = 2*pi*f/c0
      ! Cutoff wavenumber
      kc = rectwg_kc(a,b,m,n)


      ! guided wavenumber
      ! TODO : shall we treat differently the case of complex valued (negative sqrt) ? Or do fortran do the thing correctly ?
      rectwg_gamma = sqrt(one*(kc**2 - k0**2))

  end function rectwg_gamma

  !
  ! Characteristic Impedance Zc of a rectangular waveguide
  !
  ! NB 22/8/12 : the sign of the real part of the complex square root is important for the .
  !              correct value determination. It seems that the calculation using gamma is better
  !               that the one using beta
  ! TODO ???
  !
  function rectwg_Zc(f,a,b,m,n,mode) result(Zc)
    real, intent(in)     :: f, a, b
    integer,intent(in)   :: m, n
    character(len=1), intent(in)  :: mode
    complex              :: Zc

    complex  :: beta, gamma_mn
    real     :: k0

    ! vacuum wavenumber (evaluated here, in order to avoid using a global parameter)
    k0 = 2*pi*f/c0
    ! guided wavenumber
    beta = rectwg_beta(f,a,b,m,n)

    gamma_mn = rectwg_gamma(f,a,b,m,n)

    select case (mode)

      case ('H')
        print*,'H mode (TE)'
        !Zc = k0*Z0/beta
        Zc = j*k0*Z0/gamma_mn

      case ('E')
        print*,'E mode (TM)'
        !print*,'beta=',rectwg_beta(f,a,b,m,n)
        !print*,'gamma=', rectwg_gamma(f,a,b,m,n)
        !Zc = beta*Z0/k0
        Zc = gamma_mn*Z0/(j*k0)

    end select

  end function rectwg_Zc


  !
  !
  ! Spectral eigenvector function
  !
  subroutine spectralEigenVector(ny,nz,aa,bb,m,n,yy,zz,mode,Ey,Ez,Hy,Hz)
    use aloha2d_globalParameters, only : k0

    implicit none

    real, intent(in) :: nz,ny,aa,bb,yy,zz
    integer, intent(in) :: m,n
    character(len=1), intent(in)  :: mode
    complex,intent(out) :: Ey,Ez,Hy,Hz

    real :: coeff_mode_y, coeff_mode_z
    complex :: TF_cos_m, TF_cos_n, TF_sin_m, TF_sin_n

    ! Normalization factor
    select case (mode)
      case ('H')
        coeff_mode_y=(n/bb)/SQRT(m**2*bb/aa+n**2*aa/bb)
        coeff_mode_z=-(m/aa)/SQRT(m**2*bb/aa+n**2*aa/bb)

        If (m.GT.0) Then
          coeff_mode_y=coeff_mode_y*SQRT2
          coeff_mode_z=coeff_mode_z*SQRT2
        Endif

        If (n.GT.0) Then
          coeff_mode_y=coeff_mode_y*SQRT2
          coeff_mode_z=coeff_mode_z*SQRT2
        Endif

      case ('E')
        coeff_mode_y=-2.*(m/aa)/SQRT(m**2*bb/aa+n**2*aa/bb)
        coeff_mode_z=-2.*(n/bb)/SQRT(m**2*bb/aa+n**2*aa/bb)

    end select

    ! Fourier transform of the spatial eingenfunction cos and sin
     TF_sin_m = (-m*pi/(aa*k0**2))* &
         (1.-(-1.)**m*exp(-j*k0*ny*aa))* &
         1.0/(ny**2-(m*pi/(k0*aa))**2)* &
         exp(-j*k0*ny*yy)

     TF_cos_m=(-j*ny/k0)*&
         (1.-(-1.)**m*exp(-j*k0*ny*aa))*&
         1/(ny**2-(m*pi/(k0*aa))**2)*&
         exp(-j*k0*ny*yy)

     TF_sin_n=(-n*pi/(bb*k0**2))*&
         (1.-(-1.)**n*exp(-j*k0*nz*bb))*&
         1/(nz**2-(n*pi/(k0*bb))**2)*&
         exp(-j*k0*nz*zz)

     TF_cos_n=(-j*nz/k0)*&
         (1.-(-1.)**n*exp(-j*k0*nz*bb))*&
         1/(nz**2-(n*pi/(k0*bb))**2)*&
         exp(-j*k0*nz*zz)

    ! Spectral compoments
    Hy = -coeff_mode_z*TF_sin_m*TF_cos_n
    Hz = coeff_mode_y*TF_cos_m*TF_sin_n

    Ez = coeff_mode_z*TF_sin_m*TF_cos_n
    Ey = coeff_mode_y*TF_cos_m*TF_sin_n
  end subroutine spectralEigenVector


end module aloha2D_rectangularwaveguides
