      MODULE densread_seltz_mpi_I
      INTERFACE
      SUBROUTINE densread_seltz_mpi (DINT1,DINT2,DINT3,DINT4,DINT5,DINT6, &
         DINT7,DINT1VEC,DENS1VEC,NRNUC)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def, ONLY: NNNW
      USE prnt_C, ONLY: NVEC
      INTEGER, INTENT(IN) :: NRNUC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      REAL(DOUBLE), DIMENSION(NVEC,NRNUC), INTENT(OUT) :: DENS1VEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NRNUC), INTENT(IN) :: DINT1VEC
      END SUBROUTINE
      END INTERFACE
      END MODULE
