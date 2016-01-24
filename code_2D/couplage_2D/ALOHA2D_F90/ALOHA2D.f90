! MAIN FILE
!
! ALOHA 2D F90 version
! Author: J.Hillairet
! August 2012
!
PROGRAM ALOHA_2D
  use aloha2d_constants         ! math and physical constants
  use aloha2d_globalparameters  ! Input parameters to be read from the ascii file and global variables
  use aloha2d_plasma            ! general plasma routines
  use aloha2d_plasma_admittance ! plasma admittance routines
  use aloha2d_rectangularwaveguides ! retangular waveguide routines
  use aloha2d_coupling          ! Coupling element evaluation

  implicit none

  ! Constants
  character(len=*), parameter :: input_file  ='ALOHA2D.in'
  character(len=*), parameter :: output_file ='ALOHA2D.out'



  ! namelists
  !namelist /general/ wg_nb, wg_modes_nb, f, nz_min,nz_max, ny_min, ny_max
  !namelist /plasma/ ne, dne, B0
  !namelist /output_param/ K, Zc_he


  ! Get the input parameters from the input ascii file
  ! The variable memory allocation is also done inside this subroutine
  call get_input_parameters()

  ! Ion species density and gradient
  ! TODO : make this cleaner (as input parameter instead)
  nD=0.9*ne
  dnD=0.9*dne
  nT=0.
  dnT=0.
  nH=0.1*ne
  dnH=0.1*dne
  nHe=0.
  dnHe=0.

  ! Calculate the cyclotron pulsations :
  ! wce, wcD, wcT, wcH, wcHe
  call eval_cyclotron_pulsations()
  write(*,*) 'wce, wcD, wcT, wcH, wcHe =', wce, wcD, wcT, wcH, wcHe

  ! Calculate the plasma pulsations
  ! wpe, wpD, wpT, wpH, wpHe
  call eval_plasma_pulsations()
  write(*,*) 'wpe, wpD, wpT, wpH, wpHe =', wpe, wpD, wpT, wpH, wpHe

  ! Calculate the plasma pulsations derivatives
  ! with respect to species's density
  ! dwpe, dwpD, dwpT, dwpH, dwpHe
  call eval_plasma_pulsations_derivatives()
  write(*,*) 'dwpe, dwpD, dwpT, dwpH, dwpHe =', dwpe, dwpD, dwpT, dwpH, dwpHe

  ! Calculate the cold plasma dielectric tensor element (Stix elements)
  ! S, D, P
  ! and derivatives with respect to the electron density
  ! dD, dP
  call eval_cold_tensor_elements()
  write(*,*) 'S, D, P, dD, dP=', S, D, P, dD, dP

  ! Main part
  call eval_coupling()

  ! ! Quick test to verify the way arrays are stored in Namelist
  !K(1,1) = one
  !K(2,1) = 2*one
  !K(1,2) = j
  !K(2,2) = 2*j
  !
  ! ! NB: The array in fortran is written inline with the following way!
  ! ! Array content
  ! ! K = 1,2
  ! !   3,4
  ! ! Array in the namelist ==>
  ! ! K = 1,3,2,4

  !
  ! saving output parameters
  !
  call set_output_parameters()
  write(*,*) size(K)

  !
  !
  !
  !
  ! Main file Subroutines
  !
  !
  !
  !
  contains
    !
    ! Get the input parameters from the input ascii file
    !
    subroutine get_input_parameters()
      character(len=255) :: path
      logical :: exist_file

      integer :: fu, ios = 0, line = 0
      ! ATTENTION : the following char size is the maximum char to be read on a line
      ! Thus it determines the maximum number of guides the code can support !!
      ! At that part of the programm, the number of waveguides has not yet be read. So we must infer the maximum number of char we could read !
      ! Assuming that we need 9 char for one  value (8 char + 1 space, eg. " 5.88E-03"), and assuming than 300 waveguides
      ! is the maximum number of waveguides we could encouter, then 9*300=2700 should be sufficient...
      character(len=2700) :: buffer, param_value, param_label

      ! Display the current path
      call getcwd(path)
      write(*,*) 'Current path : ', path

      ! Test if the input file exists
      inquire(File=input_file, Exist=exist_file)
      if (exist_file .eqv. .FALSE.) then
        write(*,*) 'ERROR: the input file does not exist'
        stop
      endif

      ! process config input file
      !
      open(unit=fu, file=input_file, form='formatted', status='old', action='read')

      ! ios is:
      ! - negative if an end of record or endfile condition is encoutered
      ! - positive if an error is detected
      ! - zero otherwise
      do while(ios == 0)
        ! read a line as a character
        read(fu,'(A)',iostat=ios) buffer

        if (ios == 0) then
          line = line + 1

          ! Parameter name in the control file :
          ! find first instance of '=' and split data label and data value
          param_value = adjustl(buffer(index(buffer,'=')+1:))
          param_label = trim(buffer(:index(buffer,'=')-1))

          ! convert value string into the appropriate type
          ! depending of the parameter name
          select case (param_label)
            case ('wg_nb')
            ! convert substring into integer
            read(param_value,*) wg_nb
            write(*,*) '+ wg_nb = ', wg_nb


          case('wg_modes_nb')
            read(param_value,*) wg_modes_nb
            write(*,*) '+ wg_mode_nb = ', wg_modes_nb

            ! at this point we can allocate the arrays of waveguide dimensions,
            ! assuming this has not been made before ; this imply that the number
            ! of waveguide and number of mode is read BEFORE the reading of the array values...

            ! allocate the array to the number of waveguides
            allocate(a(wg_nb), b(wg_nb), y(wg_nb), z(wg_nb))
            ! Characteristic impedances of each port
            allocate(Zc_he(wg_nb*wg_modes_nb))
            ! Coupling matrix
            allocate(K(wg_nb*wg_modes_nb, wg_nb*wg_modes_nb))

          case('f')
            ! convert substring into float
            read(param_value, *) f
            write(*,*) '+ f = ', f
            ! calculate the free spacewavenumber
            k0 = 2*pi*f/c0

          case('nz_min')
            read(param_value,*) nz_min
            write(*,*) '+ nz_min = ', nz_min

          case('nz_max')
            read(param_value,*) nz_max
            write(*,*) '+ nz_max = ', nz_max

          case('ny_min')
            read(param_value,*) ny_min
            write(*,*) '+ ny_min = ', ny_min

           case('ny_max')
            read(param_value,*) ny_max
            write(*,*) '+ ny_max = ', ny_max

          case('ny_nb')
            read(param_value,*) ny_nb
            write(*,*) '+ ny_nb = ', ny_nb

          case('nz_nb')
            read(param_value,*) nz_nb
            write(*,*) '+ nz_nb = ', nz_nb

          case('b')
            read(param_value, *) b
            write(*,'((A),(f7.2))') '+ Waveguides width in [mm], b=', b*1e3

          case('a')
            read(param_value, *) a
            write(*,'((A),(f7.2))') '+ Waveguides height in [mm], a=', a*1e3

          case('y')
            read(param_value, *) y
            write(*,'((A),(f7.2))') '+ Waveguides polidal location in [mm], y=', y*1e3

          case('z')
            read(param_value, *) z
            write(*,'((A),(f7.2))') '+ Waveguides toroidal location in [mm], z=', z*1e3

          case('ne')
            read(param_value, *) ne
            write(*,*) '+ Edge density [m^-3], ne0=', ne

          case('dne')
            read(param_value, *) dne
            write(*,*) '+ Edge density gradient [m^-4], dne=', dne

          case('B0')
            read(param_value, *) B0
            write(*,*) '+ B0 = ', B0

          end select ! buffer parsing
        end if ! ios = 0
      end do ! parsing file

      close(fu)

      print*,'+============ End of input parameters ============+'

    end subroutine get_input_parameters

    !
    ! Set the output parameters into the output ascii file
    !
    subroutine set_output_parameters()
        integer :: p,q, fu
        ! newunit is a fortran 2008 Feature
        !open(newunit=fu, file='ALOHA2D.out.K.dat', form='formatted', status='replace')
        fu=1
        open(fu, file='ALOHA2D.out.K.dat', form='formatted', status='replace')
        !write(fu,*) ((real(K(p,q)),p=1,wg_nb*wg_modes_nb),q=1,wg_nb*wg_modes_nb)
        !write(fu,*) ((imag(K(p,q)),p=1,wg_nb*wg_modes_nb),q=1,wg_nb*wg_modes_nb)
        write(fu,*) 'K(real)   K(imag)'
        write(fu,'(2g20.10)') ((K(p,q),p=1,wg_nb*wg_modes_nb),q=1,wg_nb*wg_modes_nb)
        close(fu)

        fu=1
        !open(newunit=fu, file='ALOHA2D.out.Zc.dat', form='formatted', status='replace')
        open(fu, file='ALOHA2D.out.Zc.dat', form='formatted', status='replace')
        !write(fu,*) (real(Zc_he(p)),p=1,wg_nb*wg_modes_nb)
        !write(fu,*) (imag(Zc_he(p)),p=1,wg_nb*wg_modes_nb)
        write(*,*) 'Zc(real)  Zc(imag)'
        write(fu,'(2g20.10)') (Zc_he(p),p=1,wg_nb*wg_modes_nb)
        close(fu)
    end subroutine set_output_parameters

    !
    ! Main calculation procedure
    !
    subroutine eval_coupling()
      use aloha2d_plasma_admittance

      implicit none

      integer :: id_wg, id_port, id_port1, id_port2, id_mode
      integer :: fu, fu2 ! file descriptors
      integer :: q ! array index
      integer, dimension(wg_nb*wg_modes_nb)  :: m, n ! port modal indexes
      real, dimension(wg_nb*wg_modes_nb) :: a_port, b_port, y_port, z_port ! port waveguide dimensions
      character(len=1), dimension(wg_nb*wg_modes_nb) :: mode_port ! mode type index ('E' (TM) or 'H' (TE) modes)

      complex, dimension(:), allocatable :: eyt_ny_nz, ezt_ny_nz, hyt_ny_nz, hzt_ny_nz
      allocate(eyt_ny_nz(GRID_NY_NB*GRID_NZ_NB), ezt_ny_nz(GRID_NY_NB*GRID_NZ_NB))
      allocate(hyt_ny_nz(GRID_NY_NB*GRID_NZ_NB), hzt_ny_nz(GRID_NY_NB*GRID_NZ_NB))

!
!          ! Julien
!          ! 24/03/2013
!          !
!          ! Temporary piece of code for AOLGA comparizon
!          ! For the 2013 RF Conference
!          !
!          ! Assuming the number of mode is 7:
!          character(len=1), dimension(7) :: AOLGA_mode_port
!          integer, dimension(7) :: AOLGA_m, AOLGA_n
!          ! We want here to compare the following modes
!          !        ( TE10,TE01,TE11,TE02,TE12,TM11,TM12)
!
!          AOLGA_mode_port = (/ 'H', 'H', 'H', 'H', 'H', 'E', 'E' /)
!                  AOLGA_m = (/  1 ,  0 ,  1 ,  0 ,  1 ,  1 ,  1  /)
!                  AOLGA_n = (/  0 ,  1 ,  1 ,  2 ,  2 ,  1 ,  2  /)


      print*,'Preparing the calculation...'
      ! For all waveguides
      !
      ! Set up the parameters before the coupling calculations
      do id_wg=1,wg_nb
        print*,'[ Waveguide #', id_wg, '/', wg_nb, ']'

        ! for all the modes
        !
        ! Evaluate the charactertic impedances of all the port
        ! This do loop also set up the port indexing array id_port and
        ! the mode index (m,n) for all the ports
        !
        ! TODO : Evaluating the additional modes by sorting the cut-off wavelength
        !
        do id_mode=1,wg_modes_nb
            print*,'   mode #', id_mode
            ! index of the port corresponding to the current waveguide
            id_port=(id_wg-1)*wg_modes_nb + id_mode

          ! set the port's waveguide dimensions and location
          a_port(id_port) = a(id_wg)
          b_port(id_port) = b(id_wg)
          y_port(id_port) = y(id_wg)
          z_port(id_port) = z(id_wg)

          ! TODO : here we may select the modes to use depending of
          ! their cut-off wavelengths

!          ! 24/03/2013
!          ! Temporary piece of code for the comparizon to AOLGA
!          mode_port(id_port) = AOLGA_mode_port(id_mode)
!          m(id_port) = AOLGA_m(id_mode)
!          n(id_port) = AOLGA_n(id_mode)

          ! Here the first mode is the TE10 mode, while other modes are TM_1n
          if (id_mode .EQ. 1) then
            ! TE_10 (H_10) characteristic impedance
            mode_port(id_port) = 'H'
            m(id_port)=1
            n(id_port)=0
          else
            ! TM_1n (E_1n) characteristic impedance
            mode_port(id_port) = 'E'
            m(id_port)=1
            n(id_port)=id_mode-1
          end if

          ! define the port characteristic impedance depending the kind of mode is used
          if (mode_port(id_port) .EQ. 'H') then ! H (TE) modes
            Zc_he(id_port) = rectwg_Zc(real(f),a_port(id_port),b_port(id_port),m(id_port),n(id_port),'H')
          else ! E (TM) modes
            Zc_he(id_port) = rectwg_Zc(real(f),a_port(id_port),b_port(id_port),m(id_port),n(id_port),'E')
          endif

          print*, 'Zc_he(id_port)=', Zc_he(id_port)

        end do ! id_mode
      end do ! id_wg

      !
      ! Evaluate the plasma admittance on a definite fixed grid.
      ! The future needed values of the plasma admittance will thus be
      ! calculated from an interpolation of this grid
      !

      print*,'Plasma admittance grid evaluation...'
      call eval_plasma_admittance_ongrid()
      isAdmittanceCalculated = .TRUE.
      print*,'Plasma admittance grid evaluation : Done.'

      ! write the results in a text file
      ! which contains also the kind of mode
      fu=1
      open(fu, file='ALOHA2D.out.K_details.dat', form='formatted', status='replace')
      ! newunit is a fortran 2008 feature
      !open(newunit=fu, file='ALOHA2D.out.K_details.dat', form='formatted', status='replace')
      write(fu,*) 'port1  mode1  m1 n1  port2  mode2 m2  n2   K(port1,port2)'

      print*,'Launching coupling calculation...'
      ! For all ports, double loop
      ! Array storage in Fortran is column-major ; in order to optimize the loop,
      ! we access element column-wise ("the first index vary fastest").
      ! In this particular case, this play little importance, since each loop takes time to compute.
!      !$OMP parallel do
!      do id_port2=1,wg_nb*wg_modes_nb
!        !print*,' Port#1 : ',id_port1,'/',wg_nb*wg_modes_nb
!        do id_port1=1,wg_nb*wg_modes_nb
!          !print*,'++++++++ Port#2 : ',id_port2,'/',wg_nb*wg_modes_nb
!          K(id_port1,id_port2) = eval_K( &
!                a_port(id_port1),b_port(id_port1), &
!                y_port(id_port1), z_port(id_port1), &
!                mode_port(id_port1), m(id_port1),n(id_port1), &
!                a_port(id_port2),b_port(id_port2), &
!                y_port(id_port2), z_port(id_port2), &
!                mode_port(id_port2), m(id_port2),n(id_port2))
!          write(*,'(A,I3,A,I3,A,2g15.5)') 'K(',id_port1,',',id_port2,')=',K(id_port1,id_port2)
!
!          write(fu,*) id_port1, mode_port(id_port1), m(id_port1), n(id_port1), &
!                      id_port2, mode_port(id_port2), m(id_port2), n(id_port2), &
!                      K(id_port1,id_port2)
!
!        end do ! id_port2
!      end do ! id_port1
!      !$OMP end parallel do
      print*,'Coupling calculation : Done.'
      close(fu)

      print*,'size(K)=',size(K)
      print*,'size(Zc_he)=',size(Zc_he)

      print*,'Starting Spectrum calculation...'
      ! Spectral field for all ny,nz
      fu2=1
      open(fu2, file='ALOHA2D.out.spectralFields.dat', form='formatted', status='replace')
      write(fu2,*) '==================== ALOHA 2D result file - 2D Spectral Ey,Ez,Hy,Hz components ========'
      write(fu2,*) '    nbre_modes, nbre_guides'
      write(fu2,'(2g20.10)') wg_modes_nb, wg_nb
      write(fu2,*) '    nbre_ny, nbre_nz'
      write(fu2,'(2g20.10)') GRID_NY_NB, GRID_NZ_NB
      write(fu2,*) '   Ey, Ez, Hy, Hz (in spectral domain)'
      do id_port1=1,wg_nb*wg_modes_nb
            ! calculates the spectral components of the E-H field and write them into the output file
            call spect_discr(a_port(id_port1),b_port(id_port1), &
                             y_port(id_port1), z_port(id_port1), &
                             mode_port(id_port1), m(id_port1),n(id_port1), &
                             eyt_ny_nz, ezt_ny_nz, hyt_ny_nz, hzt_ny_nz)
            do q=1,GRID_NY_NB*GRID_NZ_NB
                write(fu2,*) real(eyt_ny_nz(q)), imag(eyt_ny_nz(q)), &
                             real(ezt_ny_nz(q)), imag(ezt_ny_nz(q)), &
                             real(hyt_ny_nz(q)), imag(hyt_ny_nz(q)), &
                             real(hzt_ny_nz(q)), imag(hzt_ny_nz(q))
            end do ! q

      end do !id_port


      print*,'Spectrum calculation : Done.'
      close(fu2)
    end subroutine eval_coupling

END PROGRAM ALOHA_2D
