! Stopping policy for the independent rmcdhf_mem_strict_mpi executable.
! SCNSTY is the legacy occupation-weighted update measure, not an
! Euler-equation residual or a difference between saved normalized orbitals.
      MODULE STRICT_SCF
      USE vast_kind_param, ONLY: DOUBLE
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: UPDATE_STRICT_SCF

      CONTAINS

      SUBROUTINE UPDATE_STRICT_SCF(EOL, NIT, ACCY, SCNSTY, FIXED, &
                                  ENERGY, PREVIOUS_ENERGY, STREAK, &
                                  MAX_SCNSTY, RELATIVE_CHANGE, &
                                  ORBITAL_OK, ENERGY_OK, CONVERGED)
      LOGICAL, INTENT(IN) :: EOL, FIXED(:)
      INTEGER, INTENT(IN) :: NIT
      REAL(DOUBLE), INTENT(IN) :: ACCY, SCNSTY(:)
      REAL(DOUBLE), INTENT(IN) :: ENERGY, PREVIOUS_ENERGY
      INTEGER, INTENT(INOUT) :: STREAK
      REAL(DOUBLE), INTENT(OUT) :: MAX_SCNSTY, RELATIVE_CHANGE
      LOGICAL, INTENT(OUT) :: ORBITAL_OK, ENERGY_OK, CONVERGED
      INTEGER :: J

      ORBITAL_OK = IEEE_IS_FINITE(ACCY) .AND. ACCY > 0.D0
      MAX_SCNSTY = 0.D0
      DO J = 1, SIZE(SCNSTY)
         IF (FIXED(J)) CYCLE
         IF (.NOT.IEEE_IS_FINITE(SCNSTY(J))) THEN
            ORBITAL_OK = .FALSE.
            MAX_SCNSTY = HUGE(1.D0)
         ELSE
            MAX_SCNSTY = MAX(MAX_SCNSTY, ABS(SCNSTY(J)))
            IF (SCNSTY(J) < 0.D0 .OR. SCNSTY(J) > ACCY) &
               ORBITAL_OK = .FALSE.
         ENDIF
      END DO

! EAL has no diagonalized weighted energy. Require two orbital passes.
      ENERGY_OK = .NOT.EOL
      RELATIVE_CHANGE = HUGE(1.D0)
      IF (EOL .AND. NIT > 1) THEN
         IF (IEEE_IS_FINITE(ENERGY) .AND. &
             IEEE_IS_FINITE(PREVIOUS_ENERGY)) THEN
            IF (ABS(ENERGY) > TINY(1.D0)) THEN
               RELATIVE_CHANGE = ABS((ENERGY-PREVIOUS_ENERGY)/ENERGY)
               ENERGY_OK = IEEE_IS_FINITE(RELATIVE_CHANGE) .AND. &
                           RELATIVE_CHANGE < 0.001D0*ACCY
            ENDIF
         ENDIF
      ENDIF

! A fresh restart does not inherit a convergence streak from its .w file.
      IF (NIT <= 1) STREAK = 0
      IF (ORBITAL_OK .AND. ENERGY_OK) THEN
         STREAK = MIN(STREAK + 1, 2)
      ELSE
         STREAK = 0
      ENDIF
      CONVERGED = STREAK >= 2
      END SUBROUTINE UPDATE_STRICT_SCF
      END MODULE STRICT_SCF
