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
!...Created by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, pointer :: jqsar(:, :, :)
      integer(int8), pointer :: csf_subshell_quantum_codes(:, :, :)
      integer, pointer :: jcupar(:, :)
      integer(int8), pointer :: csf_coupling_quantum_codes(:, :)
      end module stat_c
