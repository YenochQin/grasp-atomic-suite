module symmetry_block_offsets
    implicit none (type, external)

    ! Module state
    integer :: total_csf_count, mixing_coefficient_storage_size
    ! Each offset counts entries in preceding blocks. Add a one-based block-local index to
    ! obtain the position in the combined array.
    integer, pointer :: csf_offsets_by_block(:)
    integer, pointer :: eigenvalue_offsets_by_block(:), mixing_offsets_by_block(:)

end module symmetry_block_offsets
