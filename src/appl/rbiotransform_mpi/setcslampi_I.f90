      module setcsla_i
      implicit none (type, external)
      interface
      subroutine setcsla (name, ncore,igg)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer :: ncore,igg
      end subroutine setcsla
      end interface
      end module setcsla_i
