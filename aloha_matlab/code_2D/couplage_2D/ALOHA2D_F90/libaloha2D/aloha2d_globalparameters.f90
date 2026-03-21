! Global parameters for ALOHA2D
!
!
module aloha2d_globalParameters
  use aloha2d_constants

  implicit none

  ! debug
  logical, parameter :: DEBUG=.TRUE.

  ! Input parameters
  integer :: wg_nb
  integer :: wg_modes_nb
  real :: B0, f
  real :: ne, dne ! Electron density and gradient
  real :: nD, dnD, nT, dnT, nH, dnH, nHe, dnHE ! ion species density and gradient
  real, allocatable, dimension(:)  :: a, b, y, z
  real :: ny_min, nz_min
  real :: ny_max, nz_max
  real :: nz_nb, ny_nb

  ! Output parameters
  complex, allocatable, dimension(:) :: Zc_he
  complex, allocatable, dimension(:,:) :: K

  ! Usual global parameters
  real  :: k0 ! vacuum wavenumber

  ! Plasma admittance matrix elements
  complex, dimension(:), allocatable :: Ys_yy, YS_yz, Ys_zy, Ys_zz
  ! Boolean
  logical :: isAdmittanceCalculated = .FALSE.

!contains


end module aloha2d_globalParameters
