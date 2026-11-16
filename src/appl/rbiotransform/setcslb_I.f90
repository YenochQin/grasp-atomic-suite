      module setcslb_i
      implicit none (type, external)
      interface
      subroutine setcslb (name, ncore,igg)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer :: ncore,igg
      end subroutine
      end interface
      end module
