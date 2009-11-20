SUBROUTINE DSPSVX1_F95(A, B, X, UPLO, AF, IPIV, FACT, &
                      FERR, BERR, RCOND, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: SPSVX_F77 => LA_SPSVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, FACT
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND, FERR, BERR
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN) :: A(:), B(:)
   REAL(WP), INTENT(OUT) :: X(:)
   INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: IPIV(:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: AF(:)
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SPSVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LFACT, LUPLO
   INTEGER :: LINFO, N, NN, ISTAT, ISTAT1, SIPIV, SAF
   REAL(WP) :: LRCOND, LFERR, LBERR
   COMPLEX(WP) :: WW
!  .. LOCAL POINTERS ..
   INTEGER, POINTER :: IWORK(:), LPIV(:)
   REAL(WP),  POINTER :: WORK(:), LAF(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; NN = SIZE(A)
   WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
   IF( PRESENT(RCOND) ) RCOND = 1.0_WP
   IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
   IF( PRESENT(AF) )THEN; SAF = SIZE(AF); ELSE; SAF = NN; END IF
!  .. TEST THE ARGUMENTS
    IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
   ELSE IF( SIZE(B) /= N )THEN; LINFO = -2
   ELSE IF( SIZE(X) /= N )THEN; LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
   ELSE IF( SAF /= NN ) THEN; LINFO = -5
   ELSE IF( SIPIV /= N )THEN; LINFO = -6
   ELSE IF( ( .NOT. LSAME(LFACT,'F') .AND. .NOT. LSAME(LFACT,'N') ) .OR. &
     ( LSAME(LFACT,'F') .AND. .NOT.( PRESENT(AF) .AND. PRESENT(IPIV) ) ) )THEN; LINFO = -7
   ELSE IF ( N > 0 )THEN
      IF( .NOT.PRESENT(AF) ) THEN; ALLOCATE( LAF(NN), STAT=ISTAT )
      ELSE; LAF => AF; END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(IPIV) )THEN; ALLOCATE( LPIV(N), STAT=ISTAT )
         ELSE; LPIV => IPIV; END IF
      END IF
      IF( ISTAT == 0 )THEN
         ALLOCATE(WORK(MAX(1,3*N)), IWORK(N), STAT=ISTAT)
      END IF
      IF( ISTAT == 0 )THEN
!        .. CALL LAPACK77 ROUTINE
         CALL SPSVX_F77( LFACT, LUPLO, N, 1, A, LAF, LPIV, B, N, X, N, &
                         LRCOND, LFERR, LBERR, WORK, IWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(AF) ) DEALLOCATE( LAF, STAT=ISTAT1 )
      IF( .NOT.PRESENT(IPIV) ) DEALLOCATE( LPIV, STAT=ISTAT1 )
      IF( PRESENT(FERR) ) FERR = LFERR
      IF( PRESENT(BERR) ) BERR = LBERR
      IF( PRESENT(RCOND) ) RCOND=LRCOND
      DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DSPSVX1_F95
