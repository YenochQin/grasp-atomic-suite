      MODULE smsnew_mpi_I
      INTERFACE
      SUBROUTINE smsnew_mpi (DOIT,VINT,VINT2)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def, ONLY: NNNW
      INTEGER, INTENT(IN) :: DOIT
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      END SUBROUTINE
      END INTERFACE
      END MODULE
