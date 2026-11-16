      module pos_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer :: total_csf_count, mixing_coefficient_storage_size
      integer, pointer :: csf_offsets_by_block(:)
      integer, pointer :: eigenvalue_offsets_by_block(:), mixing_offsets_by_block(:)

      end module pos_c
