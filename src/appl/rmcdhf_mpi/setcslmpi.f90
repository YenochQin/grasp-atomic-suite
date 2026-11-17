subroutine setcslmpi(name, ncore, idblk)
   use rmcdhf_descriptor_mpi, only: load_rmcdhf_descriptors_mpi
   use hblock_c, only: symmetry_block_count
   implicit none (type, external)
   character(len=*), intent(in) :: name
   integer, intent(out) :: ncore
   character(len=8), intent(in) :: idblk(*)
   call load_rmcdhf_descriptors_mpi(name,ncore,idblk(:symmetry_block_count))
end subroutine setcslmpi
