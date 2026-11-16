      module setcslmpi_i
      implicit none (type, external)
      interface
      subroutine setcslmpi (name, ncore, idblk)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      integer :: ncore
      character (len = 8) :: idblk(*)
!VAST...Dummy argument IDBLK is not referenced in this routine.
      end subroutine setcslmpi
      end interface
      end module setcslmpi_i
