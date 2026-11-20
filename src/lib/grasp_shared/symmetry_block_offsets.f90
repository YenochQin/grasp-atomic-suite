      module symmetry_block_offsets
      implicit none (type, external)
      integer :: total_csf_count, mixing_coefficient_storage_size
      integer, pointer :: csf_offsets_by_block(:)
      integer, pointer :: eigenvalue_offsets_by_block(:), mixing_offsets_by_block(:)

      end module symmetry_block_offsets
