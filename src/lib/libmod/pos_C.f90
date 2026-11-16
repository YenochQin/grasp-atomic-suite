      module pos_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: total_csf_count, mixing_coefficient_storage_size
      integer, pointer :: csf_offsets_by_block(:)
      integer, pointer :: eigenvalue_offsets_by_block(:), mixing_offsets_by_block(:)

      end module pos_c
