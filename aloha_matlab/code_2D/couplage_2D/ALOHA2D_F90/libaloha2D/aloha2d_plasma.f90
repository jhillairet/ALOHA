! ALOHA2D Plasma module
!
!
module aloha2d_plasma
  use aloha2d_constants
  use aloha2d_globalparameters

  real :: wce, wcD, wcT, wcH, wcHe ! cyclotron pulsations
  real :: wpe, wpD, wpT, wpH, wpHe ! plasma pulsations
  real :: dwpe, dwpD, dwpT, dwpH, dwpHe ! plasma pulsations derivative with respect to ne
  real :: S, D, P, dD, dP ! Cold plasma dielectric tensor elements

contains

  !
  ! Calculate the cyclotron pulsations :
  ! wce, wcD, wcT, wcH, wcHe
  !
  subroutine eval_cyclotron_pulsations()
    wce = omega_c(qe, B0, me)
    wcD = omega_c(qe, B0, mD)
    wcT = omega_c(qe, B0, mT)
    wcH = omega_c(qe, B0, mH)
    wcHe= omega_c(2*qe, B0, mHe)
  end subroutine eval_cyclotron_pulsations

  !
  ! Calculate the plasma pulsations
  ! wpe, wpD, wpT, wpH, wpHe
  !
  subroutine eval_plasma_pulsations()
    wpe = omega_p(ne, qe, me)
    wpD = omega_p(nD, qe, mD)
    wpT = omega_p(nT, qe, mT)
    wpH = omega_p(nH, qe, mH)
    wpHe= omega_p(nHe, 2*qe, mHe)
  end subroutine eval_plasma_pulsations

  !
  ! Calculate the plasma pulsations derivatives
  ! with respect to species's density
  ! dwpe, dwpD, dwpT, dwpH, dwpHe
  !
  subroutine eval_plasma_pulsations_derivatives()
    dwpe = domega_p(dne, qe, me)
    dwpD = domega_p(dnD, qe, mD)
    dwpT = domega_p(dnT, qe, mT)
    dwpH = domega_p(dnH, qe, mH)
    dwpHe = domega_p(2*nD, 2*qe, mHe)
  end subroutine eval_plasma_pulsations_derivatives

  !
  ! Calculate the cold plasma dielectric tensor element (Stix elements)
  ! S, D, P
  ! and derivatives with respect to the electron density
  ! dD, dP
  !
  subroutine eval_cold_tensor_elements()
    S = cold_plasma_dielectric_S()
    D = cold_plasma_dielectric_D()
    P = cold_plasma_dielectric_P()
    dD = cold_plasma_dielectric_dD()
    dP = cold_plasma_dielectric_dP()
  end subroutine eval_cold_tensor_elements

  !
  ! Cyclotron pulsation
  !
  function omega_c(q,B,m)
    real, intent(in) :: q,B,M ! electric charge, DC magnetic field, mass
    real  :: omega_c

    omega_c = q*B/m
  end function omega_c

  !
  ! Plasma pulsation
  !
  function omega_p(n, q, m)
    real, intent(in) :: n, q, m
    real :: omega_p

    omega_p = sqrt(n*q**2/(Eps0*m))
  end function omega_p

  !
  ! Plasma pulsation derivative with respect to the species density
  !
  function domega_p(dn, q, m)
    real, intent(in) :: dn, q, m
    real :: domega_p

    domega_p = sqrt(dn*q**2/(Eps0*m))
  end function domega_p

  !
  ! Cold plasma dielectric tensor elements S, D, P
  !
  function cold_plasma_dielectric_S()
    real :: cold_plasma_dielectric_S

    cold_plasma_dielectric_S = 1-( wpe**2/((2*pi*f)**2-wce**2) + &
      wpD**2/((2*pi*f)**2-wcD**2) + &
      wpT**2/((2*pi*f)**2-wcT**2) + &
      wpH**2/((2*pi*f)**2-wcH**2) + &
      wpHe**2/((2*pi*f)**2-wcHe**2) )

  end function cold_plasma_dielectric_S

  function cold_plasma_dielectric_D()
    real :: cold_plasma_dielectric_D
    cold_plasma_dielectric_D= &
      -wpe**2*wce/(2*pi*f*((2*pi*f)**2 - wce**2)) + &
      wpD**2*wcD/(2*pi*f*((2*pi*f)**2 - wcD**2)) + &
      wpT**2*wcT/(2*pi*f*((2*pi*f)**2 - wcT**2)) + &
      wpH**2*wcH/(2*pi*f*((2*pi*f)**2 - wcH**2)) + &
      wpHe**2*wcHe/(2*pi*f*((2*pi*f)**2 - wcHe**2))
  end function cold_plasma_dielectric_D

  function cold_plasma_dielectric_P()
    real :: cold_plasma_dielectric_P

    cold_plasma_dielectric_P=&
      1-( (wpe/(2*pi*f))**2 + &
      (wpD/(2*pi*f))**2 + &
      (wpT/(2*pi*f))**2 + &
      (wpH/(2*pi*f))**2 + &
      (wpHe/(2*pi*f))**2 )

  end function cold_plasma_dielectric_P

  !
  ! Cold plasma tensor element derivative with respect to electron density dD, dP
  !
  function cold_plasma_dielectric_dD()
    real :: cold_plasma_dielectric_dD

    cold_plasma_dielectric_dD=&
      (-dwpe**2*wce/(2*pi*f*((2*pi*f)**2 - wce**2)) + &
      dwpD**2*wcD/(2*pi*f*((2*pi*f)**2 - wcD**2)) +   &
      dwpT**2*wcT/(2*pi*f*((2*pi*f)**2 - wcT**2)) +   &
      dwpH**2*wcH/(2*pi*f*((2*pi*f)**2 - wcH**2)) +   &
      dwpHe**2*wcHe/(2*pi*f*((2*pi*f)**2 - wcHe**2))) &
      /k0
  end function cold_plasma_dielectric_dD

  function cold_plasma_dielectric_dP()
    real :: cold_plasma_dielectric_dP

    cold_plasma_dielectric_dP=&
    -((dwpe/(2*pi*f))**2 + &
      (dwpD/(2*pi*f))**2 +   &
      (dwpT/(2*pi*f))**2 +   &
      (dwpH/(2*pi*f))**2 +   &
      (dwpHe/(2*pi*f))**2)/k0

  end function cold_plasma_dielectric_dP



end module aloha2d_plasma
