module aloha_constants
    ! Code working precision
    integer, parameter :: wp = kind(1.0D0)
    ! 10 decimals numbers and a dynamic of 10^60 
!     integer, parameter :: wp = selected_real_kind(10,60)
    ! mathematical constants
    real(kind=wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510
    complex, parameter       :: ImU = (0.0,1.0) 	! Imaginary Unit
    ! physics constants
    real(kind=wp), parameter :: cl = 299792458.0     ! Vacuum light velocity
    real(kind=wp), parameter :: mu0= 4.0*pi*1.E-7    ! Vacuum permeability
    real(kind=wp), parameter :: Eps0=8.854187E-12    ! Vacuum permittivity
    real(kind=wp), parameter :: Z0 = mu0*cl          ! Vacuum impedance
    real(kind=wp), parameter :: Y0 = 1.0/Z0          ! Vacuum inductance
    ! Electron properties
    real(kind=wp), parameter :: me = 9.10938215E-31  ! Electron mass
    real(kind=wp), parameter :: qe = 1.602176487E-19 ! Electron Electric charge (modulus)
    ! misc.
    real, parameter          :: pc = 1.*1.E-3

end module aloha_constants
