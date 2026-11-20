subroutine setcsl(name, core_subshell_count, block_labels)
   use rmcdhf_csf_input, only: load_rmcdhf_descriptors
   use hamiltonian_symmetry_blocks, only: symmetry_block_count
   implicit none (type, external)
   character(len=*), intent(in) :: name
   integer, intent(out) :: core_subshell_count
   character(len=8), intent(in) :: block_labels(*)
   integer :: status
   call load_rmcdhf_descriptors(name,core_subshell_count,block_labels(:symmetry_block_count),status)
   if (status /= 0) error stop "RMCDHF descriptor input failed"
end subroutine setcsl
