!***********************************************************************
!                                                                      *
      PROGRAM RIS_MPI
!                                                                      *
!   Entry routine for the MPI RIS program.                             *
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY: DOUBLE
      USE default_C
      USE iounit_C
      USE debug_C,         ONLY: CUTOFF
      USE mpi_C
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcsla_I
      USE getmixblock_I
      USE getsmd_I
      USE factt_I
      USE ris_cal_mpi_I
      IMPLICIT NONE
      INTEGER :: ARGC, IARG, K, NCI, NCORE_NOT_USED, NCOUNT1
      LOGICAL :: YES, USE_ARGS
      CHARACTER(LEN=24) :: ARG
      CHARACTER :: NAME*24

      CALL MPIX_STARTUP(MYID, NPROCS, HOST, LENHOST, NCOUNT1, &
         'RIS4_MPI', 'This is the MPI RIS program', &
         'isodata, name.c, name.(c)m, name.w', &
         'name.i or name.ci', .TRUE.)

      NAME = ' '
      NCI = 1
      NDEF = 0
      USE_ARGS = .FALSE.
      CUTOFF = 1.0D-10

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
                  WRITE (ISTDE, *) 'Usage: ris4_mpi [name] [--ci|--nonci]'
                  WRITE (ISTDE, *) 'If no arguments are given, interactive mode is used.'
                  NDEF = -1
               ELSE
                  WRITE (ISTDE, *) 'Unrecognized argument: ', TRIM(ARG)
                  WRITE (ISTDE, *) 'Usage: ris4_mpi [name] [--ci|--nonci]'
                  NDEF = -1
               ENDIF
            END DO
            IF (NDEF == 0 .AND. LEN_TRIM(NAME) == 0) THEN
               WRITE (ISTDE, *) 'Usage: ris4_mpi [name] [--ci|--nonci]'
               NDEF = -1
            ENDIF
         ENDIF

         IF (.NOT.USE_ARGS .AND. NDEF == 0) THEN
            WRITE (ISTDE, *)
            WRITE (ISTDE, *) 'Default settings?'
            YES = GETYN()
            WRITE (ISTDE, *)
            IF (YES) THEN
               NDEF = 0
            ELSE
               NDEF = 1
               WRITE (ISTDE, *) 'RIS4_MPI currently supports default settings only.'
            ENDIF
         ENDIF

         IF (NDEF == 0 .AND. .NOT.USE_ARGS) THEN
    9       WRITE (ISTDE, *) 'Name of state'
            READ(*,'(A)') NAME
            K=INDEX(NAME,' ')
            IF (K.EQ.1) THEN
               WRITE (ISTDE, *) 'Names may not start with a blank'
               GOTO 9
            ENDIF
            WRITE (ISTDE, *)
            WRITE (ISTDE, *) 'Mixing coefficients from a CI calc.?'
            YES = GETYN ()
            IF (YES) THEN
               NCI = 0
            ELSE
               NCI = 1
            ENDIF
            WRITE (ISTDE, *)
         ENDIF
      ENDIF

      CALL MPI_BCAST(NDEF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      IF (NDEF /= 0) THEN
         CALL MPI_FINALIZE(IERR)
         STOP
      ENDIF
      CALL MPI_BCAST(NAME, LEN(NAME), MPI_CHARACTER, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST(NCI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)

      CALL SETDBG
      CALL SETMC
      CALL SETCON
      IF (MYID == 0) CALL SETSUM(NAME,NCI)
      CALL SETCSLA(NAME,NCORE_NOT_USED)
      CALL GETSMD(NAME)
      CALL GETMIXBLOCK(NAME,NCI)
      CALL FACTT
      CALL RIS_CAL_MPI(NAME)

      CALL MPIX_SHUTDOWN(MYID, NCOUNT1, 'RIS4_MPI')
      STOP
      END PROGRAM RIS_MPI
