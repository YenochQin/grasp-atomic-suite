      MODULE setout_gj_I
      INTERFACE
      SUBROUTINE setout_gj (name, nci, gjc_diag, dgjc_diag)
      USE vast_kind_param, ONLY: DOUBLE
      USE prnt_C, ONLY: NVEC
      CHARACTER(LEN=24), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: nci
      REAL(DOUBLE), INTENT(IN) :: gjc_diag(NVEC), dgjc_diag(NVEC)
      END SUBROUTINE
      END INTERFACE
      END MODULE
