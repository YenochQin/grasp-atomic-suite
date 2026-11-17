!
!***********************************************************************
!                                                                      *
      module stat_c
!
!***********************************************************************
!                                                                      *
      use, intrinsic :: iso_fortran_env, only: double => real64, byte => int8
      use, intrinsic :: iso_fortran_env, only: int8
      implicit none (type, external)
      private :: int8
      integer, pointer :: jqsar(:, :, :) => null()
      integer(int8), pointer :: csf_subshell_quantum_codes(:, :, :) => null()
      integer, pointer :: jcupar(:, :) => null()
      integer(int8), pointer :: csf_coupling_quantum_codes(:, :) => null()
      end module stat_c
