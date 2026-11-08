      MODULE rkco_GG_I
      INTERFACE
!
      SUBROUTINE RKCO_GG (JA,JB,CORD,INCOR,ICOLBREI)
!      USE vast_kind_param, ONLY:  DOUBLE
      USE cord_I, ONLY: cord_callback => cord
      PROCEDURE(cord_callback) :: CORD
      INTEGER, INTENT(IN) :: JA,JB,INCOR,ICOLBREI
      END SUBROUTINE
      END INTERFACE
      END MODULE
