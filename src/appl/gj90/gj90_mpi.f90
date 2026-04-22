!***********************************************************************
!                                                                      *
      PROGRAM GJ90_MPI
!                                                                      *
!   Entry routine for the MPI Landé-factor program.                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE default_C
      USE iounit_C
      USE memory_man
      USE mpi_C
      USE prnt_C,          ONLY: NVEC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE compute_gj_only_mpi_I
      USE gethfd_I
      USE getmixblock_I
      USE getyn_I
      USE setcon_I
      USE setcsla_I
      USE setdbg_I
      USE setmc_I
      USE setout_gj_I
      USE factt_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ARGC, IARG, K, NCI, NCORE_NOT_USED
      LOGICAL :: YES, USE_ARGS
      REAL(DOUBLE), DIMENSION(:), POINTER :: GJC_DIAG, DGJC_DIAG
      CHARACTER(LEN=24) :: ARG
      CHARACTER :: NAME*24
!-----------------------------------------------
!
      CALL STARTMPI(MYID, NPROCS, HOST, LENHOST)
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
                  WRITE (ISTDE, *) 'Usage: gj90_mpi [name] [--ci|--nonci]'
                  WRITE (ISTDE, *) 'If no arguments are given, interactive mode is used.'
                  NDEF = -1
               ELSE IF (IARG /= 1 .OR. ARG(1:2) == '--') THEN
                  WRITE (ISTDE, *) 'Unrecognized argument: ', TRIM(ARG)
                  WRITE (ISTDE, *) 'Usage: gj90_mpi [name] [--ci|--nonci]'
                  NDEF = -1
               ENDIF
            END DO
            IF (NDEF == 0 .AND. LEN_TRIM(NAME) == 0) THEN
               WRITE (ISTDE, *) 'Usage: gj90_mpi [name] [--ci|--nonci]'
               NDEF = -1
            ENDIF
         ENDIF
!
         WRITE (ISTDE, *)
         WRITE (ISTDE, *) 'GJ90_MPI'
         WRITE (ISTDE, *) 'This is the MPI Lande-factor program'
         WRITE (ISTDE, *) 'Input files:  isodata, name.c, name.(c)m, name.w'
         WRITE (ISTDE, *) 'Output files: name.gj or name.cgj'
         WRITE (ISTDE, *) 'Running on ', NPROCS, ' MPI ranks'
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
               WRITE (ISTDE, *) 'GJ90_MPI currently supports default settings only.'
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
      CALL SETCSLA(NAME, NCORE_NOT_USED)
      CALL GETHFD(NAME)
      CALL GETMIXBLOCK(NAME, NCI)
      CALL FACTT
!
      CALL ALLOC(GJC_DIAG, NVEC, 'GJC_DIAG', 'GJ90_MPI')
      CALL ALLOC(DGJC_DIAG, NVEC, 'DGJC_DIAG', 'GJ90_MPI')
!
      CALL COMPUTE_GJ_ONLY_MPI(GJC_DIAG, DGJC_DIAG)
      IF (MYID == 0) CALL SETOUT_GJ(NAME, NCI, GJC_DIAG, DGJC_DIAG)
!
      CALL DALLOC(GJC_DIAG, 'GJC_DIAG', 'GJ90_MPI')
      CALL DALLOC(DGJC_DIAG, 'DGJC_DIAG', 'GJ90_MPI')
!
      IF (MYID == 0) THEN
         WRITE (ISTDE, *)
         WRITE (ISTDE, *) 'GJ90_MPI: Execution complete.'
      ENDIF
!
      CALL MPI_FINALIZE(IERR)
      STOP
      END PROGRAM GJ90_MPI
