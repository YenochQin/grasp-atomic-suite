module csf_block_cutoffs
    implicit none (type, external)

    ! Module state
    integer :: active_csf_cutoffs(100)
    integer, pointer :: csf_cutoff_by_block(:)
end module csf_block_cutoffs
