subroutine setcsl(name, ncore, idblk)
   use rmcdhf_descriptor_input, only: load_rmcdhf_descriptors
   use hblock_c, only: symmetry_block_count
   implicit none (type, external)
   character(len=*), intent(in) :: name
   integer, intent(out) :: ncore
   character(len=8), intent(in) :: idblk(*)
   integer :: status
   call load_rmcdhf_descriptors(name,ncore,idblk(:symmetry_block_count),status)
   if (status /= 0) error stop "RMCDHF descriptor input failed"
end subroutine setcsl
