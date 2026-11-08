      MODULE mcpin_I
      INTERFACE
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      SUBROUTINE mcpin (NAME, startdir, IK, NTESTG, INPCI)
      CHARACTER (LEN = 24), INTENT(IN) :: NAME
      CHARACTER(LEN=128), INTENT(IN) :: startdir
      INTEGER, INTENT(IN) :: IK
      INTEGER, INTENT(IN) :: NTESTG
      INTEGER :: INPCI
      END SUBROUTINE
      END INTERFACE
      END MODULE
