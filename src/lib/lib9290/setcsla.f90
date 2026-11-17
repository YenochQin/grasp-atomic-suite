module setcsla_i
   implicit none (type, external)
   private
   public :: setcsla
contains
   subroutine setcsla(name,ncore)
      use csf_basis_input, only: load_csf_basis
      implicit none (type, external)
      character(len=24), intent(in) :: name
      integer, intent(out) :: ncore
      integer :: status
      call load_csf_basis(trim(name),ncore,status)
      if (status /= 0) error stop "Descriptor basis load failed"
   end subroutine setcsla
end module setcsla_i
