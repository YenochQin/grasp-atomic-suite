subroutine merg12(name,ncorer,ncore)
   use csf_basis_merging, only: merge_csf_bases
   implicit none (type, external)
   character(len=24), intent(in) :: name(2)
   integer, intent(in) :: ncorer,ncore
   integer :: status
   character(len=256) :: message
   call merge_csf_bases(ncorer,ncore,status,message)
   if (status /= 0) write(*,*) "CSF descriptor merge ",trim(name(1))," / ",trim(name(2)),": ",trim(message)
   if (status /= 0) error stop "CSF descriptor merge failed"
end subroutine merg12
