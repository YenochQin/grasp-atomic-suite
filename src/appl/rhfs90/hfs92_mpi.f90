!***********************************************************************
!                                                                      *
      PROGRAM HFS92_MPI
!                                                                      *
!   Entry routine for the MPI hyperfine-structure program.             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE default_C
      USE iounit_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcsla_I
      USE gethfd_I
      USE getmixblock_I
      USE strsum_I
      USE factt_I
      USE hfsgg_mpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ARGC, IARG, K, NCI, NCORE_NOT_USED, NCOUNT1
      INTEGER :: NSTAGE1, NSTAGE2, NCOUNT_RATE, NCOUNT_MAX
      LOGICAL :: YES, USE_ARGS
      REAL(DOUBLE) :: TIME_SETCSLA, TIME_GETHFD, TIME_GETMIXBLOCK, TIME_HFSGG
      REAL(DOUBLE), DIMENSION(1) :: TIMBUF
      CHARACTER(LEN=24) :: ARG
      CHARACTER :: NAME*24
!-----------------------------------------------
!
      CALL MPIX_STARTUP(MYID, NPROCS, HOST, LENHOST, NCOUNT1, &
         'RHFS_MPI', 'This is the MPI hyperfine structure program', &
         'isodata, name.c, name.(c)m, name.w', &
         'name.(c)h, name.(c)hoffd', .TRUE.)
!
      NAME = ' '
      NCI = 1
      NDEF = 0
      USE_ARGS = .FALSE.
!
      IF (MYID == 0) THEN
         ARGC = COMMAND_ARGUMENT_COUNT()
         IF (ARGC > 0) THEN
            USE_ARGS = .TRUE.
            DO IARG = 1, ARGC
               ARG = ' '
               CALL GET_COMMAND_ARGUMENT(IARG, ARG)
               IF (IARG == 1 .AND. ARG(1:2) /= '--') THEN
                  NAME = ARG
               ELSE IF (TRIM(ARG) == '--ci') THEN
                  NCI = 0
               ELSE IF (TRIM(ARG) == '--nonci') THEN
                  NCI = 1
               ELSE IF (TRIM(ARG) == '--help' .OR. TRIM(ARG) == '-h') THEN
                  WRITE (ISTDE, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
                  WRITE (ISTDE, *) 'If no arguments are given, interactive mode is used.'
                  NDEF = -1
               ELSE IF (IARG /= 1 .OR. ARG(1:2) == '--') THEN
                  WRITE (ISTDE, *) 'Unrecognized argument: ', TRIM(ARG)
                  WRITE (ISTDE, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
                  NDEF = -1
               ENDIF
            END DO
            IF (NDEF == 0 .AND. LEN_TRIM(NAME) == 0) THEN
               WRITE (ISTDE, *) 'Usage: rhfs_mpi [name] [--ci|--nonci]'
               NDEF = -1
            ENDIF
         ENDIF
!
         IF (.NOT.USE_ARGS .AND. NDEF == 0) THEN
            WRITE (ISTDE, *)
            WRITE (ISTDE, *) 'Default settings?'
            YES = GETYN()
            WRITE (ISTDE, *)
            IF (YES) THEN
               NDEF = 0
            ELSE
               NDEF = 1
               WRITE (ISTDE, *) 'RHFS_MPI currently supports default settings only.'
            ENDIF
         ENDIF
!
         IF (NDEF == 0 .AND. .NOT.USE_ARGS) THEN
   10       CONTINUE
            WRITE (ISTDE, *) 'Name of state'
            READ (*, '(A)') NAME
            K = INDEX(NAME,' ')
            IF (K == 1) THEN
               WRITE (ISTDE, *) 'Names may not start with a blank'
               GO TO 10
            ENDIF
!
            WRITE (ISTDE, *)
            WRITE (ISTDE, *) 'Mixing coefficients from a CI calc.?'
            YES = GETYN()
            IF (YES) THEN
               NCI = 0
            ELSE
               NCI = 1
            ENDIF
         ENDIF
      ENDIF
!
      CALL MPI_BCAST(NDEF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      IF (NDEF /= 0) THEN
         CALL MPI_FINALIZE(IERR)
         STOP
      ENDIF
!
      CALL MPI_BCAST(NAME, LEN(NAME), MPI_CHARACTER, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST(NCI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
!
      CALL SETDBG
      CALL SETMC
      CALL SETCON
      IF (MYID == 0) CALL SETSUM(NAME, NCI)
      CALL SYSTEM_CLOCK(NSTAGE1, NCOUNT_RATE, NCOUNT_MAX)
      CALL SETCSLA(NAME, NCORE_NOT_USED)
      CALL SYSTEM_CLOCK(NSTAGE2, NCOUNT_RATE, NCOUNT_MAX)
      TIME_SETCSLA = DBLE(NSTAGE2 - NSTAGE1) / DBLE(NCOUNT_RATE)
      TIMBUF(1) = TIME_SETCSLA
      CALL GDMAXMPI_ROOT(TIMBUF, 1)
      IF (MYID == 0) TIME_SETCSLA = TIMBUF(1)
!
      CALL SYSTEM_CLOCK(NSTAGE1, NCOUNT_RATE, NCOUNT_MAX)
      CALL GETHFD(NAME)
      CALL SYSTEM_CLOCK(NSTAGE2, NCOUNT_RATE, NCOUNT_MAX)
      TIME_GETHFD = DBLE(NSTAGE2 - NSTAGE1) / DBLE(NCOUNT_RATE)
      TIMBUF(1) = TIME_GETHFD
      CALL GDMAXMPI_ROOT(TIMBUF, 1)
      IF (MYID == 0) TIME_GETHFD = TIMBUF(1)
!
      CALL SYSTEM_CLOCK(NSTAGE1, NCOUNT_RATE, NCOUNT_MAX)
      CALL GETMIXBLOCK(NAME, NCI)
      CALL SYSTEM_CLOCK(NSTAGE2, NCOUNT_RATE, NCOUNT_MAX)
      TIME_GETMIXBLOCK = DBLE(NSTAGE2 - NSTAGE1) / DBLE(NCOUNT_RATE)
      TIMBUF(1) = TIME_GETMIXBLOCK
      CALL GDMAXMPI_ROOT(TIMBUF, 1)
      IF (MYID == 0) TIME_GETMIXBLOCK = TIMBUF(1)
      IF (MYID == 0) CALL STRSUM
      CALL FACTT
!
      CALL SYSTEM_CLOCK(NSTAGE1, NCOUNT_RATE, NCOUNT_MAX)
      CALL HFSGG_MPI
      CALL SYSTEM_CLOCK(NSTAGE2, NCOUNT_RATE, NCOUNT_MAX)
      TIME_HFSGG = DBLE(NSTAGE2 - NSTAGE1) / DBLE(NCOUNT_RATE)
      TIMBUF(1) = TIME_HFSGG
      CALL GDMAXMPI_ROOT(TIMBUF, 1)
      IF (MYID == 0) TIME_HFSGG = TIMBUF(1)
!
      IF (MYID == 0) THEN
         WRITE (6, '(A, F10.3, A)') 'RHFS_MPI setcsla wall time  (max rank): ', &
            TIME_SETCSLA, ' s'
         WRITE (6, '(A, F10.3, A)') 'RHFS_MPI gethfd wall time   (max rank): ', &
            TIME_GETHFD, ' s'
         WRITE (6, '(A, F10.3, A)') 'RHFS_MPI getmix wall time   (max rank): ', &
            TIME_GETMIXBLOCK, ' s'
         WRITE (6, '(A, F10.3, A)') 'RHFS_MPI hfsgg wall time    (max rank): ', &
            TIME_HFSGG, ' s'
      ENDIF
!
      CALL MPIX_SHUTDOWN(MYID, NCOUNT1, 'RHFS_MPI')
      STOP
      END PROGRAM HFS92_MPI
