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
      integer, pointer :: jqsar(:, :, :)
      integer(int8), pointer :: csf_subshell_quantum_codes(:, :, :)
      integer, pointer :: jcupar(:, :)
      integer(int8), pointer :: csf_coupling_quantum_codes(:, :)
      end module stat_c
