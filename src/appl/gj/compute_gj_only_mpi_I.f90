      MODULE compute_gj_only_mpi_I
      INTERFACE
      SUBROUTINE compute_gj_only_mpi (gjc_diag, dgjc_diag)
      USE vast_kind_param, ONLY: DOUBLE
      USE prnt_C, ONLY: NVEC
      REAL(DOUBLE), INTENT(OUT) :: gjc_diag(NVEC), dgjc_diag(NVEC)
      END SUBROUTINE
      END INTERFACE
      END MODULE
