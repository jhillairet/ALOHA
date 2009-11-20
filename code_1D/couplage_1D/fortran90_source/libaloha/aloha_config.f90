! ALOHA global configuration parameters
module aloha_config
    USE aloha_constants, ONLY : wp

    ! Vacuum wavenumber
    real(kind=wp)       :: k0

    ! waveguide parameters
    character, dimension(1) :: bcte
    real(kind=wp)       :: a
    real(kind=wp)       :: b(60), z(60)

    ! ? 
    real(kind=wp)       :: epsabs
    ! plasma ??
    real(kind=wp)       :: X0, D0
    real(kind=wp)       :: X1, D1, d_couche, pertes !V6
    !modmn/
    integer             :: m, n
    ! guid
    integer             :: i,j


    ! ########## numerical integration configuration
    integer :: nlimit
    real(kind=wp) ::  max_nz

    ! Mode number
    integer :: Nmhm, Nmem
    
    integer :: Nmax=10 ! nombre de mode max
    integer :: Gmax=60 ! nombre de guide max

    ! Airy
    integer :: app

    ! telnum  
    integer :: knout

end module aloha_config