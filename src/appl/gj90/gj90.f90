!***********************************************************************
!                                                                      *
      PROGRAM GJ90
!                                                                      *
!   Entry routine for the serial Landé-factor program.                 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE default_C
      USE iounit_C
      USE memory_man
      USE prnt_C,          ONLY: NVEC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE compute_gj_only_I
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
      INTEGER :: K, NCI, NCORE_NOT_USED
      LOGICAL :: YES
      REAL(DOUBLE), DIMENSION(:), POINTER :: GJC_DIAG, DGJC_DIAG
      CHARACTER :: NAME*24
!-----------------------------------------------
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'GJ90'
      WRITE (ISTDE, *) 'This is the serial Lande-factor program'
      WRITE (ISTDE, *) 'Input files:  isodata, name.c, name.(c)m, name.w'
      WRITE (ISTDE, *) 'Output files: name.gj or name.cgj'
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'Default settings?'
      YES = GETYN()
      WRITE (ISTDE, *)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF
!
   10 CONTINUE
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
!
      CALL SETDBG
      CALL SETMC
      CALL SETCON
      CALL SETCSLA(NAME, NCORE_NOT_USED)
      CALL GETHFD(NAME)
      CALL GETMIXBLOCK(NAME, NCI)
      CALL FACTT
!
      CALL ALLOC(GJC_DIAG, NVEC, 'GJC_DIAG', 'GJ90')
      CALL ALLOC(DGJC_DIAG, NVEC, 'DGJC_DIAG', 'GJ90')
!
      CALL COMPUTE_GJ_ONLY(GJC_DIAG, DGJC_DIAG)
      CALL SETOUT_GJ(NAME, NCI, GJC_DIAG, DGJC_DIAG)
!
      CALL DALLOC(GJC_DIAG, 'GJC_DIAG', 'GJ90')
      CALL DALLOC(DGJC_DIAG, 'DGJC_DIAG', 'GJ90')
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'GJ90: Execution complete.'
      STOP
      END PROGRAM GJ90
