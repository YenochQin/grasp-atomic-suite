      MODULE densnew_mpi_I
      INTERFACE
      SUBROUTINE densnew_mpi (DOIT,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def, ONLY: NNNW
      INTEGER, INTENT(IN) :: DOIT
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2, &
                                                        DINT3, DINT4, &
                                                        DINT5, DINT6, &
                                                        DINT7
      END SUBROUTINE
      END INTERFACE
      END MODULE
