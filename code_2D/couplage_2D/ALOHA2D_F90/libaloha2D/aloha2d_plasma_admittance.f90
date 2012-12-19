!
! ALOHA 2D module
!
! Contains all the routine used to evaluate the plasma admitance
! for the slow wave and the fast wave
!
module aloha2d_plasma_admittance
  use aloha2d_plasma

  implicit none

  ! g_S=dE_S/E_S , g_F=dE_F/E_F
  complex, dimension(:), allocatable :: g_S, g_F

  ! TODO : currently the mesh grid is fixed to
  ! ny = [-3,+3]
  ! nz = [0 ,24]
  ! 361501 ??= 6/0.02 * 24/0.02 + 1 + 1500??
  real, parameter :: GRID_NZ_MIN=0.0,  GRID_NY_MIN=-3.0  ! was : 0 and -3
  real, parameter :: GRID_NZ_MAX=24,   GRID_NY_MAX=+3.0  ! was : 24 and 3
  real, parameter :: GRID_DNZ=0.02,    GRID_DNY=0.02     ! was : 0.02
  ! calculate the number of point for the grid for the array allocations
  integer :: GRID_NY_NB=ceiling( (GRID_NY_MAX-GRID_NY_MIN)/GRID_DNY) + 1
  integer :: GRID_NZ_NB=ceiling( (GRID_NZ_MAX-GRID_NZ_MIN)/GRID_DNZ) + 1

contains

  !
  ! Forms the plasma admittance tensor from the g_S and g_F functions
  !
  subroutine eval_admittanceTensor(ny,nz, gg_S, gg_F, YYs_yy,YYs_yz,YYs_zy,YYs_zz)
    use aloha2d_globalParameters
    use aloha2d_constants, only : one,j,pertes

    implicit none

    real, intent(in)     :: ny, nz
    complex, intent(in)  :: gg_S, gg_F
    complex, intent(out) :: YYs_yy,YYs_yz,YYs_zy,YYs_zz



    ! Plasma admittance matrix elements
    ! For the slow wave only (ny=0)
    !YYs_yy= (0.,0.)
    !YYs_yz=-j*S*gg_S/(S-nz**2.)
    !YYs_yz=-j*(one*S+j*pertes)*gg_S/(one*(S-nz**2.) + j*pertes)
    !YYs_zy= (0.,0.)
    !YYs_zz= (0.,0.)

    !        ! TODO : check the equations

    YYs_yy= -j*ny*nz*gg_F/(one*(S-ny**2.-nz**2.)+j*pertes)&
      +j*nz*D/(one*(S-ny**2.-nz**2.)+j*pertes)

    YYs_yz=-j*(one*S+j*pertes)*gg_S/(one*(S-nz**2.)+j*pertes)&
      -j*ny**2.*nz**2.*gg_F/((one*(S-nz**2.)+j*pertes)*&
      (one*(S-ny**2.-nz**2.)+j*pertes))&
      +j*ny*nz**2.*D/&
      ((one*(S-nz**2.)+one*pertes)*&
      (one*(S-ny**2.-nz**2.)+j*pertes))

    YYs_zy=+j*(one*(S-nz**2.)+j*pertes)*gg_F/(one*(S-nz**2.)+j*pertes)&
      -j*ny*D/(one*(S-nz**2.)+j*pertes)

    YYs_zz=+j*ny*nz*gg_F/(one*(S-nz**2.)+j*pertes)&
      -j*ny**2.*nz*D/&
      ((one*(S-nz**2.)+j*pertes)*&
      (one*(S-ny**2.-nz**2.)+j*pertes))
  end subroutine eval_admittanceTensor

  !
  ! Evaluate the surface admittance Ys(ny,nz) on a grid mesh
  ! Future needs of the value Ys(ny,nz) are deduced from an interpolation of this grid.
  ! The plasma admittance functions g_S and g_F are stored in the global parameters module
  !
  ! TODO : rewrite the array as 2D array and change the order of the element
  ! accessing in the do loop ('first index vary fastest' in fortran)
  !
  subroutine eval_plasma_admittance_ongrid()
    implicit none

    integer  :: n1, n2, incr=1
    real     :: ny, nz
    complex  :: Ys_yy_temp, Ys_yz_temp, Ys_zy_temp, Ys_zz_temp, temp
    allocate(g_S(GRID_NY_NB*GRID_NZ_NB), g_F(GRID_NY_NB*GRID_NZ_NB))
    allocate(Ys_yy(GRID_NY_NB*GRID_NZ_NB), Ys_yz(GRID_NY_NB*GRID_NZ_NB), &
      Ys_zy(GRID_NY_NB*GRID_NZ_NB), Ys_zz(GRID_NY_NB*GRID_NZ_NB))

    ! OLGA comparizon
    !open(70, file='ALOHA2D.out.admittance', form='formatted', status='replace')
    !write(70,*) 'nz, ny, Y_yy, Y_yz, Y_zy, Ys_zz'

    !!$OMP do
    ! For all elements of the ny,nz grid
    Do n1=0,GRID_NZ_NB-1
      write(*,*) n1,'/',GRID_NZ_NB-1
      nz = GRID_NZ_MIN + n1*GRID_DNZ

      Do n2=0,GRID_NY_NB-1
        ny = GRID_NY_MIN + n2*GRID_DNY

        ! TODO : calculate correctly the slow and fast admittance functions from the following routines :
        g_S(incr) = eval_gS_fem(ny,nz)
        g_F(incr) = eval_gF_fem(ny,nz)

        ! Plasma admittance tensor elements
        call eval_admittanceTensor(ny,nz,g_S(incr),g_F(incr), &
          Ys_yy_temp, Ys_yz_temp, Ys_zy_temp, Ys_zz_temp)

        Ys_yy(incr) = Ys_yy_temp
        Ys_yz(incr) = Ys_yz_temp
        Ys_zy(incr) = Ys_zy_temp
        Ys_zz(incr) = Ys_zz_temp

        !write(70,*) nz, ny, Ys_yy(incr), Ys_yz(incr), Ys_zy(incr), Ys_zz(incr)
        !print*,'incr=',incr,'/',GRID_NY_NB*GRID_NZ_NB, '(ny,nz)=',ny,nz,'g_S=',g_S(incr)
        !print*,'incr=',incr,'/',GRID_NY_NB*GRID_NZ_NB, '(ny,nz)=',ny,nz,'g_F=',g_F(incr)
        !print*,'Ys_yy(incr)=',Ys_yy(incr)
        !print*,'Ys_yz(incr)=',Ys_yz(incr)
        !print*,'Ys_zy(incr)=',Ys_zy(incr)
        !print*,'Ys_zz(incr)=',Ys_zz(incr)

        incr=incr+1

      end do ! n2
    end do ! n1
    !close(70)

    write(*,*) 'The admittance has been evaluated on the (ny,nz) grid'

  end subroutine eval_plasma_admittance_ongrid



  ! Interpolate the gS and gF value from the grid on which they have been evaluated previously
  !
  ! The interpolation is done using a 2D barycentre law
  !
  ! TODO : remove completly all the magic numbers (eg. 301)
  ! TODO : make a more robust interpolation routine ? (ny,nz) ?
  Subroutine interpolate_gS_gF(ny,nz, g_S_interp, g_F_interp)

    implicit none

    real, intent(in) :: ny,nz
    complex, intent(out)  :: g_S_interp, g_F_interp

    integer :: n1,n2,nb_ny
    real  :: abs_nz
    real  :: d1,d2,d3,d4
    complex  :: f1,f2,f3,f4

    abs_nz=abs(nz)

    n1=INT(abs_nz/GRID_DNZ)+1
    n2=INT((ny-GRID_NY_MIN)/GRID_DNY)+1

    nb_ny = INT(abs(GRID_NY_MAX - GRID_NY_MIN)/GRID_DNY)+1

    f1=g_S(n2+(n1-1)*nb_ny)
    f2=g_S(n2+1+(n1-1)*nb_ny)
    f3=g_S(n2+1+n1*nb_ny)
    f4=g_S(n2+n1*nb_ny)

    d1=SQRT((abs_nz-(n1-1)*GRID_DNZ)**2. + (ny-(n2-1)*GRID_DNY-GRID_NY_MIN)**2.)
    d2=SQRT((abs_nz-(n1-1)*GRID_DNZ)**2. + (ny-n2*GRID_DNY-GRID_NY_MIN)**2.)
    d3=SQRT((abs_nz-n1*GRID_DNZ)**2. + (ny-n2*GRID_DNY-GRID_NY_MIN)**2.)
    d4=SQRT((abs_nz-n1*GRID_DNZ)**2. + (ny-(n2-1)*GRID_DNY-GRID_NY_MIN)**2.)


    g_S_interp = (d2*d3*d4*f1 + d1*d3*d4*f2 + d1*d2*d4*f3 + d1*d2*d3*f4)/ &
      (d2*d3*d4 + d1*d3*d4 + d1*d2*d4 + d1*d2*d3)

    f1=g_F(n2+(n1-1)*nb_ny)
    f2=g_F(n2+1+(n1-1)*nb_ny)
    f3=g_F(n2+1+n1*nb_ny)
    f4=g_F(n2+n1*nb_ny)


    g_F_interp = (d2*d3*d4*f1 + d1*d3*d4*f2 + d1*d2*d4*f3 + d1*d2*d3*f4)/ &
      (d2*d3*d4 + d1*d3*d4 + d1*d2*d4 + d1*d2*d3)


  end subroutine interpolate_gS_gF


  !
  ! Evaluate the gS function using the FEM method
  !
  function eval_gS_fem(ny,nz) result(gS)
    implicit none

    real, intent(in) :: ny, nz
    complex  :: gS

    complex :: KK
    Integer :: nbre_intervalles,info,ind_piv(2000)
    Integer :: incr,incr_l,incr_c


    real :: P_tab(1001),M1(2001,2001),M2(2001,2001),L(2001,2001)
    real :: M1_ref_l1(3,3),M1_ref_l2(3,3),M2_ref(3,3),L_ref(3,3)
    real :: Long,Long_0,Long_1,Long_norm,h,h_norm
    real :: Eps_eq,dEps_eq,nbre_per,nbre_eva
    Complex*16 :: AA(7,2000),BB(2000,1)

    nbre_per=25
    nbre_eva=10

    Eps_eq = (P/S)*(S-nz**2.)-ny**2.


    ! Détermination du nombre d'intervalles
    ! pour la resolution par elements finis
    !
    If ((dP/S)*(S-nz**2.).GT.0.) then ! (dP/S * (S-nz^2) > 0)

      If (Eps_eq.GT.0.) Then ! Eps_eq > 0

        Long=nbre_per*2*pi/(k0*SQRT(Eps_eq))
        Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
          (dP/S)*(S-nz**2.)*k0*Long

        ! Deprecated version
        !
        !        10 If (Long.GT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) Then
        !        Long=Long/1.5
        !        Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
          !          (dP/S)*(S-nz**2.)*k0*Long
        !        GOTO 10
        !      EndIf

        do
          If (Long.LT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) exit
          Long=Long/1.5
          Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
            (dP/S)*(S-nz**2.)*k0*Long

        end do

        nbre_intervalles=500
        dEps_eq = (dP/S)*(S-nz**2.)
        KK=(1.,0.)*dEps_eq/Eps_eq+(0.,1.)*SQRT(Eps_eq)

      Else! Eps_eq < 0

        Long_0=(ny**2.-(P/S)*(S-nz**2.))/((dP/S)*(S-nz**2.)*k0)
        If (SQRT(-Eps_eq).GT.nbre_eva/(k0*Long_0)) Then
          Long=nbre_eva/(k0*SQRT(-Eps_eq))
          nbre_intervalles=100
          KK = (0.,0.)
        Else
          Long_1=((nbre_per*2.*pi/SQRT((dP/S)*(S-nz**2.)))**(2./3.))/k0
          Long=Long_0+Long_1
          Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
            (dP/S)*(S-nz**2.)*k0*Long
          dEps_eq = (dP/S)*(S-nz**2.)
          KK=(1.,0.)*dEps_eq/Eps_eq+(0.,1.)*SQRT(Eps_eq)
          nbre_intervalles=MIN(1000, &
            MAX(600,100*(1+INT(Long_1/Long_0))))

        End If
      End if

    Else !

      If (Eps_eq.LT.0.) Then ! Eps_eq

        Long=nbre_eva/(k0*SQRT(-Eps_eq))
        Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
          (dP/S)*(S-nz**2.)*k0*Long

        ! Deprecated version
        !
        !      20         If (Long.GT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) Then
        !      Long=Long/1.5
        !      Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
          !        (dP/S)*(S-nz**2.)*k0*Long
        !      GOTO 20
        !    End If

        do
          If (Long.LT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) exit
          Long=Long/1.5
          Eps_eq = (P/S)*(S-nz**2.)-ny**2.+ &
            (dP/S)*(S-nz**2.)*k0*Long
        end do

        nbre_intervalles=100
        KK=(0.,0.)

      Else ! Eps_eq

        Long_0=(ny**2.-(P/S)*(S-nz**2.))/ &
          ((dP/S)*(S-nz**2.)*k0)

        Long_1=((nbre_eva/SQRT(-(dP/S)*(S-nz**2.)))**(2./3.))/k0
        Long=Long_0+Long_1
        KK=(0.,0.)
        nbre_intervalles=MIN(1000, &
          MAX(600,100*(1+INT(Long_0/Long_1))))

      End If ! Eps_eq
    End if !


    h=Long/nbre_intervalles
    Long_norm=k0*Long
    h_norm=k0*h

    Do incr=1,nbre_intervalles+1
      P_tab(incr)=P+dP*h_norm*(incr-1)
    end do


    Do incr_c=1,2*nbre_intervalles+1
      Do incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles+1,incr_c+2)
        M1(incr_l,incr_c)=0.
        M2(incr_l,incr_c)=0.
        L(incr_l,incr_c)=0.
      end do
    end do


    M1_ref_l1(1,1)=(h_norm/60.)*7.
    M1_ref_l1(1,2)=(h_norm/60.)*4.
    M1_ref_l1(1,3)=(h_norm/60.)*(-1)
    M1_ref_l1(2,1)=(h_norm/60.)*4.
    M1_ref_l1(2,2)=(h_norm/60.)*16.
    M1_ref_l1(2,3)=(h_norm/60.)*0.
    M1_ref_l1(3,1)=(h_norm/60.)*(-1.)
    M1_ref_l1(3,2)=(h_norm/60.)*0.
    M1_ref_l1(3,3)=(h_norm/60.)*1.

    M1_ref_l2(1,1)=(h_norm/60.)*1.
    M1_ref_l2(1,2)=(h_norm/60.)*0.
    M1_ref_l2(1,3)=(h_norm/60.)*(-1.)
    M1_ref_l2(2,1)=(h_norm/60.)*0.
    M1_ref_l2(2,2)=(h_norm/60.)*16.
    M1_ref_l2(2,3)=(h_norm/60.)*4.
    M1_ref_l2(3,1)=(h_norm/60.)*(-1.)
    M1_ref_l2(3,2)=(h_norm/60.)*4.
    M1_ref_l2(3,3)=(h_norm/60.)*7.

    M2_ref(1,1)=(h_norm/15.)*2.
    M2_ref(1,2)=(h_norm/15.)*1.
    M2_ref(1,3)=(h_norm/15.)*(-1./2.)
    M2_ref(2,1)=(h_norm/15.)*1.
    M2_ref(2,2)=(h_norm/15.)*8.
    M2_ref(2,3)=(h_norm/15.)*1.
    M2_ref(3,1)=(h_norm/15.)*(-1./2.)
    M2_ref(3,2)=(h_norm/15.)*1.
    M2_ref(3,3)=(h_norm/15.)*2.

    L_ref(1,1)=(1./(3.*h_norm))*7.
    L_ref(1,2)=(1./(3.*h_norm))*(-8.)
    L_ref(1,3)=(1./(3.*h_norm))*1.
    L_ref(2,1)=(1./(3.*h_norm))*(-8.)
    L_ref(2,2)=(1./(3.*h_norm))*16.
    L_ref(2,3)=(1./(3.*h_norm))*(-8.)
    L_ref(3,1)=(1./(3.*h_norm))*1.
    L_ref(3,2)=(1./(3.*h_norm))*(-8.)
    L_ref(3,3)=(1./(3.*h_norm))*7.


    Do incr=1,nbre_intervalles

      M1(2*(incr-1)+1,2*(incr-1)+1)=M1(2*(incr-1)+1,2*(incr-1)+1) &
        +M1_ref_l1(1,1)*P_tab(incr)+M1_ref_l2(1,1)*P_tab(incr+1)
      M1(2*(incr-1)+1,2*(incr-1)+2)=M1(2*(incr-1)+1,2*(incr-1)+2) &
        +M1_ref_l1(1,2)*P_tab(incr)+M1_ref_l2(1,2)*P_tab(incr+1)
      M1(2*(incr-1)+1,2*(incr-1)+3)=M1(2*(incr-1)+1,2*(incr-1)+3) &
        +M1_ref_l1(1,3)*P_tab(incr)+M1_ref_l2(1,3)*P_tab(incr+1)
      M1(2*(incr-1)+2,2*(incr-1)+1)=M1(2*(incr-1)+2,2*(incr-1)+1) &
        +M1_ref_l1(2,1)*P_tab(incr)+M1_ref_l2(2,1)*P_tab(incr+1)
      M1(2*(incr-1)+2,2*(incr-1)+2)=M1(2*(incr-1)+2,2*(incr-1)+2) &
        +M1_ref_l1(2,2)*P_tab(incr)+M1_ref_l2(2,2)*P_tab(incr+1)
      M1(2*(incr-1)+2,2*(incr-1)+3)=M1(2*(incr-1)+2,2*(incr-1)+3) &
        +M1_ref_l1(2,3)*P_tab(incr)+M1_ref_l2(2,3)*P_tab(incr+1)
      M1(2*(incr-1)+3,2*(incr-1)+1)=M1(2*(incr-1)+3,2*(incr-1)+1) &
        +M1_ref_l1(3,1)*P_tab(incr)+M1_ref_l2(3,1)*P_tab(incr+1)
      M1(2*(incr-1)+3,2*(incr-1)+2)=M1(2*(incr-1)+3,2*(incr-1)+2) &
        +M1_ref_l1(3,2)*P_tab(incr)+M1_ref_l2(3,2)*P_tab(incr+1)
      M1(2*(incr-1)+3,2*(incr-1)+3)=M1(2*(incr-1)+3,2*(incr-1)+3) &
        +M1_ref_l1(3,3)*P_tab(incr)+M1_ref_l2(3,3)*P_tab(incr+1)


      M2(2*(incr-1)+1,2*(incr-1)+1)=M2(2*(incr-1)+1,2*(incr-1)+1) &
        +M2_ref(1,1)
      M2(2*(incr-1)+1,2*(incr-1)+2)=M2(2*(incr-1)+1,2*(incr-1)+2) &
        +M2_ref(1,2)
      M2(2*(incr-1)+1,2*(incr-1)+3)=M2(2*(incr-1)+1,2*(incr-1)+3) &
        +M2_ref(1,3)
      M2(2*(incr-1)+2,2*(incr-1)+1)=M2(2*(incr-1)+2,2*(incr-1)+1) &
        +M2_ref(2,1)
      M2(2*(incr-1)+2,2*(incr-1)+2)=M2(2*(incr-1)+2,2*(incr-1)+2) &
        +M2_ref(2,2)
      M2(2*(incr-1)+2,2*(incr-1)+3)=M2(2*(incr-1)+2,2*(incr-1)+3) &
        +M2_ref(2,3)
      M2(2*(incr-1)+3,2*(incr-1)+1)=M2(2*(incr-1)+3,2*(incr-1)+1) &
        +M2_ref(3,1)
      M2(2*(incr-1)+3,2*(incr-1)+2)=M2(2*(incr-1)+3,2*(incr-1)+2) &
        +M2_ref(3,2)
      M2(2*(incr-1)+3,2*(incr-1)+3)=M2(2*(incr-1)+3,2*(incr-1)+3) &
        +M2_ref(3,3)

      L(2*(incr-1)+1,2*(incr-1)+1)=L(2*(incr-1)+1,2*(incr-1)+1) &
        +L_ref(1,1)
      L(2*(incr-1)+1,2*(incr-1)+2)=L(2*(incr-1)+1,2*(incr-1)+2) &
        +L_ref(1,2)
      L(2*(incr-1)+1,2*(incr-1)+3)=L(2*(incr-1)+1,2*(incr-1)+3) &
        +L_ref(1,3)
      L(2*(incr-1)+2,2*(incr-1)+1)=L(2*(incr-1)+2,2*(incr-1)+1) &
        +L_ref(2,1)
      L(2*(incr-1)+2,2*(incr-1)+2)=L(2*(incr-1)+2,2*(incr-1)+2) &
        +L_ref(2,2)
      L(2*(incr-1)+2,2*(incr-1)+3)=L(2*(incr-1)+2,2*(incr-1)+3) &
        +L_ref(2,3)
      L(2*(incr-1)+3,2*(incr-1)+1)=L(2*(incr-1)+3,2*(incr-1)+1) &
        +L_ref(3,1)
      L(2*(incr-1)+3,2*(incr-1)+2)=L(2*(incr-1)+3,2*(incr-1)+2) &
        +L_ref(3,2)
      L(2*(incr-1)+3,2*(incr-1)+3)=L(2*(incr-1)+3,2*(incr-1)+3) &
        +L_ref(3,3)



    end do



    Do incr_c=1,2*nbre_intervalles

      Do incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles,incr_c+2)

        AA(5+incr_l-incr_c,incr_c) = -(1.,0.)*L(incr_l+1,incr_c+1)        &
          + ((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/              &
          ((1.,0.)*S+(0.,1.)*pertes)*                     &
          ((1.,0.)*M1(incr_l+1,incr_c+1)+            &
          (0.,1.)*pertes*M2(incr_l+1,incr_c+1))  &
          - (1.,0.)*(ny**2.)*M2(incr_l+1,incr_c+1)

      end do ! incr_l

      BB(incr_c,1)=(0.,0.)

    end do ! incr_c


    BB(1,1)=(1.,0.)*L(2,1)                                           &
      -((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/             &
      ((1.,0.)*S+(0.,1.)*pertes)*                  &
      ((1.,0.)*M1(2,1)+(0.,1.)*pertes*M2(2,1))  &
      + (1.,0.)*(ny**2.)*M2(2,1)
    BB(2,1)=(1.,0.)*L(3,1)                                           &
      -((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)/            &
      ((1.,0.)*S+(0.,1.)*pertes)*                 &
      ((1.,0.)*M1(3,1)+(0.,1.)*pertes*M2(3,1))  &
      + (1.,0.)*(ny**2.)*M2(3,1)


    AA(5,2*nbre_intervalles)=AA(5,2*nbre_intervalles)+KK

    Call ZGBSV(2*nbre_intervalles,2,2,1,AA,7,ind_piv,BB,  &
      2*nbre_intervalles,info)

    gS=(1/h_norm)*((-3.,0)+4.*BB(1,1)-BB(2,1))

  end function eval_gS_fem


  !
  ! Evaluate the fast wave plasma admittance using the FEM method
  !
  function eval_gF_fem(ny,nz) result(gF)

    implicit none

    real, intent(in) :: ny,nz
    complex :: gF

    complex :: KK
    Integer :: nbre_intervalles,info,ind_piv(2000)
    Integer :: incr,incr_l,incr_c

    real :: D_carre_tab(1001)
    real :: M1_bis(2001,2001),M2(2001,2001),L(2001,2001)
    real :: M1_ref_l1(3,3),M1_ref_l2(3,3),M2_ref(3,3),L_ref(3,3)
    real :: Long,Long_0,Long_1,Long_norm,h,h_norm
    real :: Eps_eq,dEps_eq,nbre_per,nbre_eva
    complex*16 :: A_bis(7,2000),B_bis(2000,1)


    nbre_per=25
    nbre_eva=10

    If ((D/dD).LT.0.) then
      write(*,*) 'warning D/dD<0 non traite dans g_f'
    Endif


    Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
      (dD**2./(nz**2.-S))*((D/dD)**2.)


    If ((nz**2.).GT.S) Then ! (nz^2 > S)

      If (Eps_eq.GT.0.) Then ! Eps_eq > 0

        Long=nbre_per*2*pi/(k0*SQRT(Eps_eq))
        Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        ! Deprecated version

        !        10 If (Long.GT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) Then
        !          Long=Long/1.5
        !          Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          !            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        !        GOTO 10
        !        end if

        do
          If (Long.LT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) exit
          Long=Long/1.5
          Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        end do


        nbre_intervalles=500
        dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
        KK=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*SQRT(Eps_eq)

      Else ! Eps_eq < 0

        Long_0=(SQRT((ny*dD/(S-nz**2.)-(S-ny**2.-nz**2.))/ &
          (dD**2./(nz**2.-S)))-D/dD)/k0
        If (SQRT(-Eps_eq).GT.nbre_eva/(k0*Long_0)) Then
          Long=nbre_eva/(k0*SQRT(-Eps_eq))
          nbre_intervalles=100
          KK = (0.,0.)
        Else
          Long_1=((nbre_per*2*pi/ &
            SQRT((2*(dD**2.)/(nz**2.-S))* &
            (k0*Long_0+D/dD)))**(2./3.))/k0
          Long=Long_0+Long_1
          Eps_eq=S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)

          ! deprecated version

          ! 15            If (Long_1.GT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) Then
          !                  Long_1=Long_1/1.5
          !                  Long=Long_0+Long_1
          !                  Eps_eq=S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
            !                          (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
          !                  GOTO 15
          !               EndIf

          do
            If (Long_1.LT.1.5*nbre_per*2*pi/(k0*SQRT(Eps_eq))) exit
            Long_1=Long_1/1.5
            Long=Long_0+Long_1
            Eps_eq=S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
              (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
          end do

          dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
          KK=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*SQRT(Eps_eq)
          nbre_intervalles=MIN(1000, MAX(600,100*(1+INT(Long_1/Long_0))))
        End If

      end if ! eps

    else ! (nz^2 < S)

      If (Eps_eq.LT.0.) Then ! Eps_eq < 0

        Long=nbre_eva/(k0*SQRT(-Eps_eq))
        Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)

        ! deprecated version
        ! 20         If (Long.GT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) Then
        !               Long=Long/1.5
        !               Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          !                 (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        !               GOTO 20
        !            EndIf

        do
          If (Long.LT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) exit
          Long=Long/1.5
          Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        end do

        nbre_intervalles=100
        KK=(0.,0.)

      else! Eps_eq > 0

        Long_0=(SQRT((ny*dD/(S-nz**2.)-(S-ny**2.-nz**2.))/ &
          (dD**2./(nz**2.-S)))-D/dD)/k0
        Long_1=((nbre_eva/ &
          SQRT((-2*(dD**2.)/(nz**2.-S))* &
          (k0*Long_0+D/dD)))**(2./3.))/k0
        Long=Long_0+Long_1
        Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)

        ! deprecated version
        ! 25         If (Long_1.GT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) Then
        !               Long_1=Long_1/1.5
        !               Long=Long_0+Long_1
        !               Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
          !                           (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        !               GOTO 25
        !            EndIf

        do
          If (Long_1.LT.1.5*nbre_eva/(k0*SQRT(-Eps_eq))) exit
          Long_1=Long_1/1.5
          Long=Long_0+Long_1
          Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+ &
            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
        end do


        KK=(0.,0.)
        nbre_intervalles=MIN(1000,MAX(600,100*(1+INT(Long_0/Long_1))))

      end if ! Eps_eq
    end if ! (nz^2 <> S)

    !c$$$      Long=0.2
    !c$$$      nbre_intervalles=1000
    !c$$$
    !c$$$      If ((nz**2.).GT.S) Then
    !c$$$         Eps_eq =S-ny**2.-nz**2.-ny*dD/(S-nz**2.)+
    !c$$$     &            (dD**2./(nz**2.-S))*((k0*Long+D/dD)**2.)
    !c$$$         dEps_eq=2*(dD**2./(nz**2.-S))*(k0*Long+D/dD)
    !c$$$         If (Eps_eq.GT.0.) Then
    !c$$$             K=(1.,0.)*dEps_eq/Eps_eq-(0.,1.)*SQRT(Eps_eq)
    !c$$$         Else
    !c$$$             K=(0.,0.)
    !c$$$         Endif
    !c$$$      Else
    !c$$$         K = (0.,0.)
    !c$$$      Endif

    h=Long/nbre_intervalles
    Long_norm=k0*Long
    h_norm=k0*h

    Do incr=1,nbre_intervalles+1
      D_carre_tab(incr)=(D+dD*h_norm*(incr-1))**2.
    end do

    Do incr_c=1,2*nbre_intervalles+1
      Do incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles+1,incr_c+2)
        M1_bis(incr_l,incr_c)=0.
        M2(incr_l,incr_c)=0.
        L(incr_l,incr_c)=0.
      end do
    end do

    M1_ref_l1(1,1)=(h_norm/60.)*7.
    M1_ref_l1(1,2)=(h_norm/60.)*4.
    M1_ref_l1(1,3)=(h_norm/60.)*(-1)
    M1_ref_l1(2,1)=(h_norm/60.)*4.
    M1_ref_l1(2,2)=(h_norm/60.)*16.
    M1_ref_l1(2,3)=(h_norm/60.)*0.
    M1_ref_l1(3,1)=(h_norm/60.)*(-1.)
    M1_ref_l1(3,2)=(h_norm/60.)*0.
    M1_ref_l1(3,3)=(h_norm/60.)*1.

    M1_ref_l2(1,1)=(h_norm/60.)*1.
    M1_ref_l2(1,2)=(h_norm/60.)*0.
    M1_ref_l2(1,3)=(h_norm/60.)*(-1.)
    M1_ref_l2(2,1)=(h_norm/60.)*0.
    M1_ref_l2(2,2)=(h_norm/60.)*16.
    M1_ref_l2(2,3)=(h_norm/60.)*4.
    M1_ref_l2(3,1)=(h_norm/60.)*(-1.)
    M1_ref_l2(3,2)=(h_norm/60.)*4.
    M1_ref_l2(3,3)=(h_norm/60.)*7.

    M2_ref(1,1)=(h_norm/15.)*2.
    M2_ref(1,2)=(h_norm/15.)*1.
    M2_ref(1,3)=(h_norm/15.)*(-1./2.)
    M2_ref(2,1)=(h_norm/15.)*1.
    M2_ref(2,2)=(h_norm/15.)*8.
    M2_ref(2,3)=(h_norm/15.)*1.
    M2_ref(3,1)=(h_norm/15.)*(-1./2.)
    M2_ref(3,2)=(h_norm/15.)*1.
    M2_ref(3,3)=(h_norm/15.)*2.

    L_ref(1,1)=(1./(3.*h_norm))*7.
    L_ref(1,2)=(1./(3.*h_norm))*(-8.)
    L_ref(1,3)=(1./(3.*h_norm))*1.
    L_ref(2,1)=(1./(3.*h_norm))*(-8.)
    L_ref(2,2)=(1./(3.*h_norm))*16.
    L_ref(2,3)=(1./(3.*h_norm))*(-8.)
    L_ref(3,1)=(1./(3.*h_norm))*1.
    L_ref(3,2)=(1./(3.*h_norm))*(-8.)
    L_ref(3,3)=(1./(3.*h_norm))*7.


    Do incr=1,nbre_intervalles

      M1_bis(2*(incr-1)+1,2*(incr-1)+1)= &
        M1_bis(2*(incr-1)+1,2*(incr-1)+1) &
        +M1_ref_l1(1,1)*D_carre_tab(incr) &
        +M1_ref_l2(1,1)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+1,2*(incr-1)+2)= &
        M1_bis(2*(incr-1)+1,2*(incr-1)+2) &
        +M1_ref_l1(1,2)*D_carre_tab(incr) &
        +M1_ref_l2(1,2)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+1,2*(incr-1)+3)= &
        M1_bis(2*(incr-1)+1,2*(incr-1)+3) &
        +M1_ref_l1(1,3)*D_carre_tab(incr) &
        +M1_ref_l2(1,3)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+2,2*(incr-1)+1)= &
        M1_bis(2*(incr-1)+2,2*(incr-1)+1) &
        +M1_ref_l1(2,1)*D_carre_tab(incr) &
        +M1_ref_l2(2,1)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+2,2*(incr-1)+2)= &
        M1_bis(2*(incr-1)+2,2*(incr-1)+2) &
        +M1_ref_l1(2,2)*D_carre_tab(incr) &
        +M1_ref_l2(2,2)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+2,2*(incr-1)+3)= &
        M1_bis(2*(incr-1)+2,2*(incr-1)+3) &
        +M1_ref_l1(2,3)*D_carre_tab(incr) &
        +M1_ref_l2(2,3)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+3,2*(incr-1)+1)= &
        M1_bis(2*(incr-1)+3,2*(incr-1)+1) &
        +M1_ref_l1(3,1)*D_carre_tab(incr) &
        +M1_ref_l2(3,1)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+3,2*(incr-1)+2)= &
        M1_bis(2*(incr-1)+3,2*(incr-1)+2) &
        +M1_ref_l1(3,2)*D_carre_tab(incr) &
        +M1_ref_l2(3,2)*D_carre_tab(incr+1)
      M1_bis(2*(incr-1)+3,2*(incr-1)+3)= &
        M1_bis(2*(incr-1)+3,2*(incr-1)+3) &
        +M1_ref_l1(3,3)*D_carre_tab(incr) &
        +M1_ref_l2(3,3)*D_carre_tab(incr+1)

      M2(2*(incr-1)+1,2*(incr-1)+1)=M2(2*(incr-1)+1,2*(incr-1)+1)+M2_ref(1,1)
      M2(2*(incr-1)+1,2*(incr-1)+2)=M2(2*(incr-1)+1,2*(incr-1)+2)+M2_ref(1,2)
      M2(2*(incr-1)+1,2*(incr-1)+3)=M2(2*(incr-1)+1,2*(incr-1)+3)+M2_ref(1,3)
      M2(2*(incr-1)+2,2*(incr-1)+1)=M2(2*(incr-1)+2,2*(incr-1)+1)+M2_ref(2,1)
      M2(2*(incr-1)+2,2*(incr-1)+2)=M2(2*(incr-1)+2,2*(incr-1)+2)+M2_ref(2,2)
      M2(2*(incr-1)+2,2*(incr-1)+3)=M2(2*(incr-1)+2,2*(incr-1)+3)+M2_ref(2,3)
      M2(2*(incr-1)+3,2*(incr-1)+1)=M2(2*(incr-1)+3,2*(incr-1)+1)+M2_ref(3,1)
      M2(2*(incr-1)+3,2*(incr-1)+2)=M2(2*(incr-1)+3,2*(incr-1)+2)+M2_ref(3,2)
      M2(2*(incr-1)+3,2*(incr-1)+3)=M2(2*(incr-1)+3,2*(incr-1)+3)+M2_ref(3,3)

      L(2*(incr-1)+1,2*(incr-1)+1)=L(2*(incr-1)+1,2*(incr-1)+1)+L_ref(1,1)
      L(2*(incr-1)+1,2*(incr-1)+2)=L(2*(incr-1)+1,2*(incr-1)+2)+L_ref(1,2)
      L(2*(incr-1)+1,2*(incr-1)+3)=L(2*(incr-1)+1,2*(incr-1)+3)+L_ref(1,3)
      L(2*(incr-1)+2,2*(incr-1)+1)=L(2*(incr-1)+2,2*(incr-1)+1)+L_ref(2,1)
      L(2*(incr-1)+2,2*(incr-1)+2)=L(2*(incr-1)+2,2*(incr-1)+2)+L_ref(2,2)
      L(2*(incr-1)+2,2*(incr-1)+3)=L(2*(incr-1)+2,2*(incr-1)+3)+L_ref(2,3)
      L(2*(incr-1)+3,2*(incr-1)+1)=L(2*(incr-1)+3,2*(incr-1)+1)+L_ref(3,1)
      L(2*(incr-1)+3,2*(incr-1)+2)=L(2*(incr-1)+3,2*(incr-1)+2)+L_ref(3,2)
      L(2*(incr-1)+3,2*(incr-1)+3)=L(2*(incr-1)+3,2*(incr-1)+3)+L_ref(3,3)



    end do


    Do incr_c=1,2*nbre_intervalles
      Do incr_l=MAX(1,incr_c-2),MIN(2*nbre_intervalles,incr_c+2)

        A_bis(5+incr_l-incr_c,incr_c) = -(1.,0.)*L(incr_l+1,incr_c+1) &
          + ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes)) &
          *M1_bis(incr_l+1,incr_c+1) &
          + ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes- &
          (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes)) &
          *M2(incr_l+1,incr_c+1)

      end do
      B_bis(incr_c,1)=(0.,0.)

    end do


    B_bis(1,1)=(1.,0.)*L(2,1) &
      - ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes))*M1_bis(2,1) &
      - ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes- &
      (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes))*M2(2,1)

    B_bis(2,1)=(1.,0.)*L(3,1) &
      - ((1.,0.)/((1.,0.)*(nz**2.-S)-(0.,1.)*pertes))*M1_bis(3,1) &
      - ((1.,0.)*(S-ny**2.-nz**2.)+(0.,1.)*pertes- &
      (1.,0.)*ny*dD/((1.,0.)*(S-nz**2.)+(0.,1.)*pertes))*M2(3,1)

    A_bis(5,2*nbre_intervalles)=A_bis(5,2*nbre_intervalles)+KK


    Call ZGBSV(2*nbre_intervalles,2,2,1,A_bis,7,ind_piv,B_bis, &
      2*nbre_intervalles,info)

    gF=(1/h_norm)*((-3.,0)+4.*B_bis(1,1)-B_bis(2,1))


  End function eval_gF_fem



end module aloha2d_plasma_admittance
