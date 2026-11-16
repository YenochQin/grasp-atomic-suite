      module oscl_i
      implicit none (type, external)
      interface
      subroutine oscl (name)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      end subroutine
      end interface
      end module
