subroutine setcslmpi(name, ncore, idblk)
   use mpi_rmcdhf_csf_input, only: load_rmcdhf_descriptors_mpi
   use hamiltonian_symmetry_blocks, only: symmetry_block_count
   implicit none (type, external)
   character(len=*), intent(in) :: name
   integer, intent(out) :: ncore
   character(len=8), intent(in) :: idblk(*)
   call load_rmcdhf_descriptors_mpi(name,ncore,idblk(:symmetry_block_count))
end subroutine setcslmpi
