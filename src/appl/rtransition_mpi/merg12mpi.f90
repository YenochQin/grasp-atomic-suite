subroutine merg12(name,ncorer,ncore)
   use csf_basis_merging, only: merge_csf_bases
   use mpi_csf_collectives, only: check_descriptor_status
   implicit none (type, external)
   character(len=128), intent(in) :: name(2)
   integer, intent(in) :: ncorer,ncore
   integer :: status
   character(len=256) :: message
   call merge_csf_bases(ncorer,ncore,status,message)
   if (status /= 0) write(*,*) "CSF descriptor merge ",trim(name(1))," / ",trim(name(2)),": ",trim(message)
   call check_descriptor_status(status)
end subroutine merg12
