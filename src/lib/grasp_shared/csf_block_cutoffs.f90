module csf_block_cutoffs
    implicit none (type, external)

    ! Module state
    integer :: iccut(100)
    integer, pointer :: iccutblk(:)
end module csf_block_cutoffs
