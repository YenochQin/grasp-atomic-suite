!***********************************************************************
!                                                                      *
      SUBROUTINE COMPUTE_GJ_ONLY(GJC_DIAG, DGJC_DIAG)
!                                                                      *
!   This routine evaluates only the diagonal ASF contributions needed  *
!   for g_J, delta g_J, and total g_J.                                 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE decide_C
      USE EIGV_C
      USE foparm_C,        ONLY: ICCUT
      USE iounit_C,        ONLY: ISTDE
      USE orb_C,           ONLY: NCF, NW
      USE prnt_C,          ONLY: NVEC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      USE itjpo_I
      USE matelt_I
      USE oneparticlejj_I
      USE rint_I
      USE rinthf_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(OUT) :: GJC_DIAG(NVEC), DGJC_DIAG(NVEC)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, IPT, IC, LCNUM, IR, ITJPOC, ITJPOR, IDIFF, IA, IB, &
         K, LOC
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: RINTGJ, RINTDGJ, GJMELT, DGJMELT
      REAL(DOUBLE) :: APART, GJPART, DGJPART, ELEMNTGJ, ELEMNTDGJ, CONTRGJ, &
         CONTRDGJ
      CHARACTER :: CNUM*11
!-----------------------------------------------
!
      GJC_DIAG(:NVEC) = 0.0D00
      DGJC_DIAG(:NVEC) = 0.0D00
!
!   Precompute the radial integrals and angular factors required by
!   the two one-body operators contributing to the Landé factor.
!
      DO I = 1, NW
         DO J = 1, NW
            RINTGJ(I,J) = RINTHF(I,J,1)
            RINTDGJ(I,J) = RINT(I,J,0)
            CALL MATELT(I, 1, J, APART, GJPART, DGJPART)
            GJMELT(I,J) = GJPART
            DGJMELT(I,J) = DGJPART
         END DO
      END DO
!
      IPT = 1
!
!   Sweep over the CSF pairs. Only the diagonal-in-J branch from HFSGG
!   is required, and only ASF diagonal matrix elements are accumulated.
!
      DO IC = 1, NCF
         IF (MOD(IC,100) == 0) THEN
            CALL CONVRT(IC, CNUM, LCNUM)
            WRITE (ISTDE, *) 'Column '//CNUM(1:LCNUM)//' complete;'
         ENDIF
!
         DO IR = 1, NCF
            IF (LFORDR .AND. IC>ICCUT .AND. IC/=IR) CYCLE
!
            ITJPOC = ITJPO(IC)
            ITJPOR = ITJPO(IR)
            IDIFF = ITJPOC - ITJPOR
            IF (.NOT.(IDIFF == 0 .AND. IR >= IC)) CYCLE
!
            ELEMNTGJ = 0.0D00
            ELEMNTDGJ = 0.0D00
!
            CALL ONEPARTICLEJJ(1, IPT, IC, IR, IA, IB, TSHELL)
            IF (IA == 0) CYCLE
!
            IF (IA == IB) THEN
               DO I = 1, NW
                  IF (ABS(TSHELL(I)) <= CUTOFF) CYCLE
                  ELEMNTGJ = ELEMNTGJ + GJMELT(I,I)*RINTGJ(I,I)*TSHELL(I)
                  ELEMNTDGJ = ELEMNTDGJ + DGJMELT(I,I)*RINTDGJ(I,I)*TSHELL(I)
               END DO
            ELSE
               IF (ABS(TSHELL(1)) <= CUTOFF) CYCLE
               ELEMNTGJ = ELEMNTGJ + GJMELT(IA,IB)*RINTGJ(IA,IB)*TSHELL(1)
               ELEMNTDGJ = ELEMNTDGJ + DGJMELT(IA,IB)*RINTDGJ(IA,IB)*TSHELL(1)
            ENDIF
!
            IF (ABS(ELEMNTGJ) <= CUTOFF .AND. ABS(ELEMNTDGJ) <= CUTOFF) CYCLE
!
            DO K = 1, NVEC
               LOC = (K - 1)*NCF
               IF (IR /= IC) THEN
                  CONTRGJ = ELEMNTGJ*(EVEC(IC + LOC)*EVEC(IR + LOC) + EVEC(IR &
                     + LOC)*EVEC(IC + LOC))
                  CONTRDGJ = ELEMNTDGJ*(EVEC(IC + LOC)*EVEC(IR + LOC) + EVEC(IR&
                      + LOC)*EVEC(IC + LOC))
               ELSE
                  CONTRGJ = ELEMNTGJ*EVEC(IC + LOC)*EVEC(IR + LOC)
                  CONTRDGJ = ELEMNTDGJ*EVEC(IC + LOC)*EVEC(IR + LOC)
               ENDIF
               GJC_DIAG(K) = GJC_DIAG(K) + CONTRGJ
               DGJC_DIAG(K) = DGJC_DIAG(K) + CONTRDGJ
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE COMPUTE_GJ_ONLY
