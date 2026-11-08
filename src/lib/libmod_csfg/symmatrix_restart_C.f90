      module symmatrix_restart_C
! CYC 2023/01/11
! Collection all variables 

      USE csfg_parameters, ONLY: CUTOFF => CSFG_RESTART_CUTOFF
      USE vast_kind_param, ONLY:  DOUBLE, LONG
! CYC 2023/01/11

      implicit none

! Restart mode (and, reuse the existing old.cm)
      INTEGER       :: IFNCOL= 27 
      INTEGER       :: ncsfDF1_core
      INTEGER(LONG) :: NELMNTres
      CHARACTER*256 :: FOLDCM
      LOGICAL       :: LOLDCM, RESTRT
      INTEGER       :: NEVBLK_OLD, NCFBLK_OLD
      INTEGER       :: IATJPO_OLD, IASPAR_OLD
      INTEGER, DIMENSION(:), ALLOCATABLE          :: IVEC_OLD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EVAL_OLD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EVEC_OLD


      end module symmatrix_restart_C
