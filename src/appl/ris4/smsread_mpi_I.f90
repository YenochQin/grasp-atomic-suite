      MODULE smsread_mpi_I
      INTERFACE
      SUBROUTINE smsread_mpi (VINT,VINT2)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def, ONLY: NNNW
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      END SUBROUTINE
      END INTERFACE
      END MODULE
