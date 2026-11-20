subroutine merg12(state_name,initial_core_subshell_count,final_core_subshell_count)
   use csf_basis_merging, only: merge_csf_bases
   use mpi_csf_collectives, only: check_descriptor_status
   implicit none (type, external)
   character(len=128), intent(in) :: state_name(2)
   integer, intent(in) :: initial_core_subshell_count,final_core_subshell_count
   integer :: status
   character(len=256) :: message
   call merge_csf_bases(initial_core_subshell_count,final_core_subshell_count,status,message)
   if (status /= 0) write(*,*) "CSF descriptor merge ",trim(state_name(1))," / ",trim(state_name(2)),": ",trim(message)
   call check_descriptor_status(status)
end subroutine merg12
