! mathematical and physical constants module
module aloha2d_constants

  implicit none

  !--! specific precisions, usually same as real and double precision
  integer, parameter :: r6 = selected_real_kind(6)
  integer, parameter :: r15 = selected_real_kind(15)
  ! Code working precision
  integer, parameter :: wp = 8!kind(1.0D0)
  ! 10 decimals numbers and a dynamic of 10^60
  !     integer, parameter :: wp = selected_real_kind(10,60)

  ! mathematical constants
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510
  complex, parameter :: j = (0.0,1.0)   ! Imaginary Unit shortcut
  complex, parameter       :: one = (1.0, 0.0)! Real Unit shortcut
  real, parameter          :: SQRT2 = sqrt(2.0)

  ! physics constants
  real, parameter :: c0 = 299792458.0     ! Vacuum light velocity
  real, parameter :: mu0= 4.0*pi*1.E-7    ! Vacuum permeability
  real, parameter :: Eps0=8.854187E-12    ! Vacuum permittivity
  real, parameter :: Z0 = mu0*c0          ! Vacuum impedance
  real, parameter :: Y0 = 1.0/Z0          ! Vacuum inductance

  ! Electron properties
  real, parameter :: me = 9.10938215E-31  ! Electron mass
  real, parameter :: qe = 1.602176487E-19 ! Electron Electric charge (modulus)
  ! misc.
  real, parameter          :: pc = 1.*1.E-3

  ! Ion species masses
  real, parameter :: mD=3.345E-27     ! Deuton mass
  real, parameter :: mT=5.018E-2     ! Tritium mass
  real, parameter ::  mH=1.673E-27  ! Hydrogen mass
  real, parameter :: mHe=5.018E-27        ! Helium mass

  ! Had-Hoc losses
  real, parameter :: pertes = -10.E-2 ! converge from values > 5E-2
  ! Integration error
  ! relative accuracy requested
  real, parameter :: epsrel = 0.001
  ! absolute accuracy requested
  real, parameter :: epsabs = 0.01
  integer, parameter :: quad_key = 6

end module aloha2d_constants
